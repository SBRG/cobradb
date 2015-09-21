# -*- coding: utf-8 -*-

from ome import settings, timing, base
from ome.base import *
from ome.components import Gene, Protein
from ome.util import scrub_gene_id, get_or_create_data_source

import sys, os, math, re
from warnings import warn
from sqlalchemy import text, or_, and_, func
import logging


class BadGenomeError(Exception):
    pass


def get_bioproject_id(genbank_filepath, fast=False):
    """Load the file and return the Bioproject ID and the loaded file.

    Returns a tuples of (BioProject ID, SeqIO object).

    Arguments
    ---------

    genbank_filepath: The path to the genbank file.

    fast: If True, then only look in the first 100 lines for a BioProject ID,
    and return a tuple of (BioProject ID, None)

    """
    # imports
    from Bio import SeqIO

    if fast:
        # try to find the BioProject ID in the first 100 lines. Otherwise, use
        # the full SeqIO.read
        line_limit = 100
        with open(genbank_filepath, 'r') as f:
            for i, line in enumerate(f.readlines()):
                s = line.split('BioProject: ')
                if len(s) > 1:
                    return s[1].strip(), None
                if i > line_limit:
                    break
        logging.debug('Could not use get_bioproject_id in fast mode. Falling back.')

    # load the genbank file
    logging.debug('Loading file: %s' % genbank_filepath)
    try:
        gb_file = SeqIO.read(genbank_filepath, 'gb')
    except IOError:
        raise BadGenomeError("File '%s' not found" % genbank_filepath)
    except Exception as e:
        raise BadGenomeError('BioPython failed to parse %s with error "%s"' %
                             (genbank_filepath, e.message))

    bioproject_id = None
    for value in gb_file.dbxrefs[0].split():
        if 'BioProject' in value:
            bioproject_id = value.split(':')[1]

    if bioproject_id is None:
        raise BadGenomeError('Invalid genbank file %s: Does not contain a BioProject ID' % genbank_filepath)

    return bioproject_id, gb_file


def load_gene_synonym(session, gene_db, synonym, data_source_name):
    """Load the synonym for this gene from the given genome."""

    data_source_id = get_or_create_data_source(session, data_source_name)
    synonym_db = (session
                  .query(Synonym)
                  .filter(Synonym.ome_id == gene_db.id)
                  .filter(Synonym.synonym == synonym)
                  .filter(Synonym.data_source_id == data_source_id)
                  .first())

    if synonym_db is None:
        synonym_db = Synonym(type='gene',
                             ome_id=gene_db.id,
                             synonym=synonym,
                             data_source_id=data_source_id)
        session.add(synonym_db)
        session.flush()

    return synonym_db.id


def _get_qual(feat, name, get_first=False):
    """Get a non-null attribute from the feature."""
    try:
        qual = feat.qualifiers[name]
    except KeyError:
        if get_first:
            return None
        else:
            return []

    def nonempty_str(s):
        s = s.strip()
        return None if s == '' else s

    if get_first:
        return nonempty_str(qual[0])
    else:
        return [y for y in (nonempty_str(x) for x in qual)
                if y is not None]


@timing
def load_genome(genbank_filepath, session):
    bioproject_id, gb_file = get_bioproject_id(genbank_filepath)
    old_locus_tag_id = get_or_create_data_source(session, 'refseq_old_locus_tag')
    organism = gb_file.annotations['organism']
    genome = (session
              .query(base.Genome)
              .filter(base.Genome.bioproject_id == bioproject_id)
              .filter(base.Genome.organism == organism)
              .first())

    if not genome:
        logging.debug('Adding new genome: %s' % bioproject_id)
        ome_genome = { 'bioproject_id': bioproject_id,
                       'organism': organism,
                       'taxon_id':'None' }
        genome = base.Genome(**ome_genome)
        session.add(genome)
        session.flush()
    else:
        logging.debug('Genome already loaded for bioproject_id %s' % bioproject_id)

    chromosome = (session
                  .query(base.Chromosome)
                  .filter(base.Chromosome.genbank_id == gb_file.annotations['gi'])
                  .filter(base.Chromosome.genome_id == genome.id)
                  .first())

    if not chromosome:
        logging.debug('Adding new chromosome: %s' % gb_file.annotations['gi'])
        chromosome = base.Chromosome(genome_id=genome.id,
                                     genbank_id=gb_file.annotations['gi'],
                                     ncbi_id=gb_file.id)
        session.add(chromosome)
        session.flush()
    else:
        logging.debug('Chromosome already loaded: %s' % gb_file.annotations['gi'])

    bigg_id_warnings = 0
    duplicate_genes_warnings = 0
    warning_num = 5
    for i, feature in enumerate(gb_file.features):
        # get the source information
        if feature.type == 'source':
            for ref in _get_qual(feature, 'db_xref'):
                if 'taxon' == ref.split(':')[0]:
                    genome.taxon_id = ref.split(':')[1]
                    break
            continue

        # only read in CDSs
        if feature.type != 'CDS':
            continue

        # bigg_id required
        bigg_id = None
        gene_name = None
        refseq_name = None
        locus_tag = None

        t = _get_qual(feature, 'locus_tag', True)
        if t is not None:
            locus_tag = t
            bigg_id = scrub_gene_id(t)

        t = _get_qual(feature, 'gene', True)
        if t is not None:
            gene_name = t
            refseq_name = t

        if gene_name is not None and bigg_id is None:
            if bigg_id_warnings <= warning_num:
                msg = 'No locus_tag for gene. Using Gene name as bigg_id: %s' % gene_name
                if bigg_id_warnings == warning_num:
                    msg += ' (Warnings limited to %d)' % warning_num
                logging.warn(msg)
                bigg_id_warnings += 1
            bigg_id = scrub_gene_id(gene_name)
            gene_name = bigg_id
        elif bigg_id is None:
            logging.error(('No locus_tag or gene name for gene %d in chromosome '
                           '%s' % (i, chromosome.genbank_id)))
            continue

        gene_db = (session
                   .query(Gene)
                   .filter(Gene.bigg_id == bigg_id)
                   .filter(Gene.chromosome_id == chromosome.id)
                   .first())
        if gene_db is None:
            # get the strand and positions
            strand = None
            if feature.strand == 1:
                strand = '+'
            elif feature.strand == -1:
                strand = '-'
            leftpos = int(feature.location.start)
            rightpos = int(feature.location.end)

            # finally, create the gene
            gene_db = Gene(bigg_id=bigg_id,
                           locus_tag=locus_tag,
                           chromosome_id=chromosome.id,
                           name=gene_name,
                           leftpos=leftpos,
                           rightpos=rightpos,
                           strand=strand,
                           mapped_to_genbank=True)
            session.add(gene_db)
            session.commit()
        else:
            # warn about duplicate genes.
            #
            # TODO The only downside to loading CDS's this way is that the
            # leftpos and rightpos correspond to a CDS, not the whole gene. So
            # these need to be fixed eventually.
            if duplicate_genes_warnings <= warning_num:
                msg = 'Duplicate genes %s on chromosome %s' % (bigg_id, chromosome.id)
                if duplicate_genes_warnings == warning_num:
                    msg += ' (Warnings limited to %d)' % warning_num
                logging.warn(msg)
                duplicate_genes_warnings += 1

        # load the synonyms for the gene
        if locus_tag is not None:
            load_gene_synonym(session, gene_db, locus_tag, 'locus_tag')

        if refseq_name is not None:
            load_gene_synonym(session, gene_db, refseq_name, 'refseq_name')

        for ref in _get_qual(feature, 'gene_synonym'):
            synonyms = [x.strip() for x in ref.split(';')]
            for syn in synonyms:
                load_gene_synonym(session, gene_db, syn, 'refseq_synonym')

        for ref in _get_qual(feature, 'db_xref'):
            splitrefs = [x.strip() for x in ref.split(':')]
            if len(splitrefs) == 2:
                load_gene_synonym(session, gene_db, splitrefs[1], splitrefs[0])

        for ref in _get_qual(feature, 'old_locus_tag'):
            for syn in [x.strip() for x in ref.split(';')]:
                load_gene_synonym(session, gene_db, syn, 'refseq_old_locus_tag')


        for ref in _get_qual(feature, 'note'):
            for value in [x.strip() for x in ref.split(';')]:
                sp = value.split(':')
                if len(sp) == 2 and sp[0] == 'ORF_ID':
                    load_gene_synonym(session, gene_db, sp[1], 'refseq_orf_id')


    session.commit()

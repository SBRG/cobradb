# -*- coding: utf-8 -*-

from ome import settings, timing, base
from ome.components import Gene, Protein

import sys, os, math, re
from warnings import warn
from sqlalchemy import text, or_, and_, func
import logging

def getAttributes(file_name):
    file = open(settings.data_directory + '/annotation/MetaCyc/17.1/data/'+file_name,'r')
    atts = []
    attflag = 0
    for line in file.readlines():
        if line[0] != '#': return atts
        if line[0:13] == '# Attributes:':
            attflag = 1
            continue

        if attflag == 1 and line != '#\n': atts.append(line.lstrip('#    ').rstrip('\n'))

def get_metacyc_data():
    metacyc_genes = parse_metacyc_dat('genes.dat')
    metacyc_promoters = parse_metacyc_dat('promoters.dat')
    metacyc_proteins = parse_metacyc_dat('proteins.dat')
    metacyc_ligands = parse_metacyc_dat('compounds.dat')
    metacyc_protein_cplxs = parse_metacyc_dat('protligandcplxes.dat')
    metacyc_tus = parse_metacyc_dat('transunits.dat')
    return (metacyc_genes, metacyc_promoters, metacyc_proteins, metacyc_ligands,
            metacyc_protein_cplxs, metacyc_tus)

def parse_metacyc_dat(file_name):
    file = open(settings.data_directory + '/annotation/MetaCyc/17.1/data/'+file_name,'r')
    master_dict = {}
    attributes = {}
    unique_id = ''

    for line in file.readlines():
        if line[0] == '#': continue

        vals = line.rstrip('\n').split(" - ")

        if vals[0] == 'UNIQUE-ID':
            master_dict[unique_id] = attributes
            unique_id = vals[1]
            attributes = {}
        try:
            tmp = attributes[vals[0]]
            tmp.append(vals[1])
            attributes[vals[0]] = tmp
        except:
            try: attributes[vals[0]] = [vals[1]]
            except: None

    return master_dict


def parse_metacyc_col(file_name):
    file = open(settings.data_directory + '/annotation/MetaCyc/17.1/data/'+file_name,'r')

    master_dict = {}
    attributes = {}
    unique_id = ''
    header = []

    for line in file.readlines():
        if line[0] == '#': continue

        if line[0:9] == 'UNIQUE-ID':
            header = line.rstrip('\n').split('\t')
            continue

        vals = line.rstrip('\n').split('\t')

        for i in range(len(vals)):
            if vals[i] == '': continue

            try:
                tmp = attributes[header[i]]
                tmp.append(vals[i])
                attributes[header[i]] = tmp
            except:
                try: attributes[header[i]] = [vals[i]]
                except: None

        master_dict[vals[0]] = attributes
        attributes = {}

    return master_dict


def parse_metacyc_fsa(file_name):
    file = open(settings.data_directory + '/annotation/MetaCyc/17.1/data/'+file_name,'r')
    id,seq = '',''
    seq_dict = {}
    for line in file.readlines():
        if line[0] == '>':
            seq_dict[id] = seq
            seq = ''
            id = line.split()[0].lstrip('>')
        else: seq+=line.rstrip('\n')
    return seq_dict


def make_citations(ome, metacyc_entry, id):
    if "CITATIONS" not in metacyc_entry:
        return
    citations = metacyc_entry['CITATIONS']
    citation_pmid_list = []
    for citation_entry in citations:
        citation_pmid = citation_entry.strip().strip("[").strip("]").split(":")[0]
        try:
            citation_pmid_list.append(int(citation_pmid))
        except:
            pass
    for citation_pmid in set(citation_pmid_list):
        citation = get_or_create(session, Citation, pmid=citation_pmid)
        ome.execute("INSERT INTO citation_list (table_id, citation_id) " + \
                        "VALUES (%d, %d);" % (id, citation.id))


def checkEqual(list): return list[1:] == list[:-1]


def uniquify(seq, idfun=None):
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result


def parse_gpr(gpr):
    gprs = []
    elems = gpr.split(' or ')
    for e in elems: gprs.append(e.replace(')','').replace('(','').split(' and '))
    return gprs


def scrub_metacyc_entry(entry, args=['UNIQUE-ID','TYPES'],extra_args=[]):
    try:
        return {arg: entry[arg] for arg in args+extra_args}
    except:
        return None


def get_gene_with_metacyc(session, base, components, genome, gene_entry):
    vals = scrub_metacyc_entry(gene_entry,
                               args=['UNIQUE-ID','ACCESSION-1','LEFT-END-POSITION',
                                     'RIGHT-END-POSITION','TRANSCRIPTION-DIRECTION',
                                     'COMMON-NAME'])
    if vals is None: return None

    gene = (session
            .query(components.Gene)
            .filter(and_(components.Gene.chromosome_id == genome.id, 
                         or_(components.Gene.bigg_id == vals['ACCESSION-1'][0], 
                             components.Gene.name == vals['COMMON-NAME'][0])))
            .first())
    if gene is None:
        print 'Exception, MetaCyc gene:'+vals['ACCESSION-1'][0]+' not found in genbank'

    return gene


def update_gene_with_metacyc(session, base, components, gene_entry):
    vals = scrub_metacyc_entry(gene_entry,
                               args=['UNIQUE-ID','ACCESSION-1','LEFT-END-POSITION',
                                     'RIGHT-END-POSITION','TRANSCRIPTION-DIRECTION', 'COMMON-NAME'])
    if vals is None: return None

    gene = get_gene_with_metacyc(session, base, components, gene_entry)
    if not gene: return None

    gene.name = vals['COMMON-NAME'][0]

    session.add(gene)
    session.flush()


def get_protein_with_metacyc(session, base, components, genome, protein_entry, metacyc_genes):
    vals = scrub_metacyc_entry(protein_entry, extra_args=['GENE'])
    if vals is None:
        return None

    # First see if the gene exists from genbank, if not log an error message and return

    try:
        gene_entry = metacyc_genes[vals['GENE'][0]]
    except:
        print 'MetaCyc issue, gene: '+vals['GENE'][0]+' does not exist'
        return None

    gene = get_gene_with_metacyc(session, base, components, genome, gene_entry)

    # If the gene exists in genbank, then the protein should exist in genbank as
    # well, if not log an error message and return

    try:
        return session.query(components.Protein).filter(components.Protein.gene_id == gene.id).one()
    except:
        return None


def update_protein_with_metacyc(session, base, components, genome, protein_entry):
    vals = scrub_metacyc_entry(protein_entry,extra_args=['GENE', 'COMMON-NAME'])
    if vals is None: return None

    protein = get_protein_with_metacyc(session, base, components, genome,
                                       protein_entry, metacyc_genes)
    if not protein: return None

    protein.name = vals['COMMON-NAME'][0]

    session.add(protein)
    session.flush()


def get_or_create_metacyc_ligand(session, base, components, ligand_entry):
    vals = scrub_metacyc_entry(ligand_entry,extra_args=['COMMON-NAME', 'SMILES'])
    if vals is None:
        return None

    name = vals['COMMON-NAME'][0].replace('\'','[prime]')
    markup = re.compile("<.+?>")
    name = markup.sub("", name)

    return session.get_or_create(components.Metabolite,
                                 bigg_id=vals['UNIQUE-ID'][0], kegg_id=None,
                                 cas_number=None, formula=None, name=name,
                                 flag=True, smiles=vals['SMILES'][0])


def get_or_create_metacyc_protein_complex(session, base, components, genome,
                                          protein_complex_entry,
                                          metacyc_proteins, metacyc_ligands,
                                          metacyc_protein_cplxs):
    #if protein_complex_entry['UNIQUE-ID'][0] == 'CPLX0-226': print protein_complex_entry
    #print protein_complex_entry['UNIQUE-ID'][0]
    vals = scrub_metacyc_entry(protein_complex_entry, extra_args=['COMMON-NAME', 'COMPONENTS'])
    if vals is None:
        return None

    protein_complex = session.get_or_create(components.Complex,
                                            bigg_id=vals['UNIQUE-ID'][0],
                                            name=vals['COMMON-NAME'][0])

    for component in vals['COMPONENTS']:

        component_vals = None
        complex_component = None

        if component in metacyc_proteins:
            component_vals = scrub_metacyc_entry(metacyc_proteins[component])

        elif component in metacyc_ligands:
            component_vals = scrub_metacyc_entry(metacyc_ligands[component])

        elif component in metacyc_protein_cplxs:
            component_vals = scrub_metacyc_entry(metacyc_protein_cplxs[component])

        if component_vals is None: continue

        if 'Protein-Complexes' in component_vals['TYPES']:
            complex_component = get_or_create_metacyc_protein_complex(session, base, components, genome, metacyc_proteins[component])

        elif 'Protein-Small-Molecule-Complexes' in component_vals['TYPES']:
            complex_component = get_or_create_metacyc_protein_complex(session, base, components, genome, metacyc_protein_cplxs[component])

        elif 'Polypeptides' in component_vals['TYPES']:
            complex_component = get_protein_with_metacyc(session, base,
                                                         components, genome,
                                                         metacyc_proteins[component],
                                                         metacyc_genes)

        elif 'Compounds' in component_vals['TYPES']:
            complex_component = get_or_create_metacyc_ligand(session, base, components, metacyc_ligands[component])

        if complex_component is None: continue

        session.get_or_create(components.ComplexComposition, complex_id=protein_complex.id,\
                                                                 component_id=complex_component.id,\
                                                                 stoichiometry=1.)
    return protein_complex


def get_or_create_metacyc_transcription_unit(session, base, components, genome,
                                             tu_entry, metacyc_genes, metacyc_promoters):

    vals = scrub_metacyc_entry(tu_entry, extra_args=['COMPONENTS', 'COMMON-NAME'])
    if vals is None:
        vals = scrub_metacyc_entry(tu_entry, extra_args=['COMPONENTS'])
    if vals is None: return None

    genes = []
    tss = None
    for metacyc_component in vals['COMPONENTS']:
        try:
            gene_entry = metacyc_genes[metacyc_component]
            genes.append(get_gene_with_metacyc(session, base, components, genome, gene_entry))
        except: None
        try:
            promoter_entry = metacyc_promoters[metacyc_component]
            tss_vals = scrub_metacyc_entry(promoter_entry,extra_args=['ABSOLUTE-PLUS-1-POS'])
            tss = tss_vals['ABSOLUTE-PLUS-1-POS'][0]
        except: None

    if len(genes) == 0: return

    try:
        leftpos = min([x.leftpos for x in genes])
        rightpos = max([x.rightpos for x in genes])
        strand = genes[0].strand
    except: return

    try:
        bigg_id = vals['UNIQUE-ID'][0]
    except:
        bigg_id = vals['COMMON-NAME'][0]
    try:
        name = vals['COMMON-NAME'][0]
    except:
        name = ''

    if not tss and strand == '+':
        tss = leftpos
    elif not tss and strand =='-':
        tss = rightpos

    if strand == '+':
        tu = session.get_or_create(components.TU, bigg_id=bigg_id, name=name,
                                   leftpos=tss, rightpos=rightpos,
                                   strand=strand, genome_id=genome.id)
    else:
        tu = session.get_or_create(components.TU, bigg_id=bigg_id, name=name,
                                   leftpos=leftpos, rightpos=tss, strand=strand,
                                   genome_id=genome.id)

    for gene in genes:
        session.get_or_create(components.TUGenes, tu_id=tu.id, gene_id=gene.id)

    return tu


@timing
def load_genome(genbank_file, session):
    # imports 
    from Bio import SeqIO

    # load the genbank file
    logging.debug('Loading file: %s' % genbank_file)
    try:
        gb_file = SeqIO.read(genbank_file, 'gb')
    except IOError:
        raise Exception("File '%s' not found" % genbank_file)
    except Exception as e:
        raise Exception('biopython failed to parse %s with error "%s"' %
                        (genbank_file, e.message))

    bioproject_id = None
    for value in gb_file.dbxrefs[0].split():
        if 'BioProject' in value:
            bioproject_id = value.split(':')[1]

    if bioproject_id is None:
        raise Exception('Invalid genbank file %s: Does not contain a BioProject ID' % genbank_file)

    genome = (session
              .query(base.Genome)
              .filter(base.Genome.bioproject_id == bioproject_id)
              .first())

    if not genome:
        logging.debug('Adding new genome: %s' % bioproject_id)
        ome_genome = { 'bioproject_id': bioproject_id,
                       'organism': gb_file.annotations['organism'] }
        genome = base.Genome(**ome_genome)
        session.add(genome)
        session.flush()

    chromosome = (session
                  .query(base.Chromosome)
                  .filter(base.Chromosome.genbank_id == gb_file.annotations['gi'])
                  .filter(base.Chromosome.genome_id == genome.id)
                  .first())

    if not chromosome:
        logging.debug('Adding new chromosome: %s' % gb_file.annotations['gi'])
        ome_chromosome = { 'genome_id': genome.id,
                           'genbank_id': gb_file.annotations['gi'],
                           'ncbi_id': gb_file.id }
        chromosome = base.Chromosome(**ome_chromosome)
        session.add(chromosome)
        session.flush()
    else:
        logging.debug('Chromosome already loaded: %s' % gb_file.annotations['gi'])

    db_xref_data_source_id = { data_source.name: data_source.id for data_source
                               in session.query(base.DataSource).all() }

    bigg_id_warnings = 0
    for i, feature in enumerate(gb_file.features):
        # only read in CDSs
        if feature.type != 'CDS':
            continue

        # bigg_id required
        bigg_id = None
        gene_name = ''

        if 'locus_tag' in feature.qualifiers:
            bigg_id = feature.qualifiers['locus_tag'][0]

        if 'gene' in feature.qualifiers:
            gene_name = feature.qualifiers['gene'][0]

        if gene_name != '' and bigg_id is None:
            warning_num = 5
            if bigg_id_warnings <= warning_num:
                msg = 'No locus_tag for gene. Using Gene name as bigg_id: %s' % gene_name
                if bigg_id_warnings == warning_num:
                    msg += ' (Warnings limited to %d)' % warning_num
                logging.warn(msg)
                bigg_id_warnings += 1
            bigg_id = gene_name
            gene_name = ''
        elif bigg_id is None:
            logging.error(('No locus_tag or gene name for gene %d in chromosome '
                           '%s' % (i, chromosome.genbank_id)))
            continue

        ome_gene = {}
        ome_gene['bigg_id'] = bigg_id
        ome_gene['name'] = gene_name
        ome_gene['leftpos'] = int(feature.location.start)
        ome_gene['rightpos'] = int(feature.location.end)
        ome_gene['chromosome_id'] = chromosome.id
        ome_gene['mapped_to_genbank'] = True
        ome_gene['info'] = ''
        if feature.strand == 1: ome_gene['strand'] = '+'
        elif feature.strand == -1: ome_gene['strand'] = '-'

        # record the notes in info
        elif 'note' in feature.qualifiers:
            ome_gene['info'] = feature.qualifiers['note'][0][0:300]
            
        # record the function in info
        if 'function' in feature.qualifiers:
            ome_gene['info'] = ome_gene['info'] + ',' + feature.qualifiers['function'][0]

        # finally, create the gene
        gene = session.get_or_create(Gene, **ome_gene)
        session.commit()

        # get the protein
        if 'protein_id' in feature.qualifiers and len(feature.qualifiers['protein_id']) > 0:

            ome_protein = {}
            ome_protein['bigg_id'] = feature.qualifiers['protein_id'][0]
            ome_protein['gene_id'] = gene.id

            if 'product' in feature.qualifiers and len(feature.qualifiers['product']) > 0:
                ome_protein['name'] = feature.qualifiers['product'][0]
            else:
                ome_protein['name'] = ''

            session.get_or_create(Protein, **ome_protein)
        
        # add the synonyms
        if 'db_xref' in feature.qualifiers:
            for ref in feature.qualifiers['db_xref']:
                if gene.id is not None:
                    splitrefs = ref.split(':')
                    ome_synonym = { 'type': 'gene' }
                    ome_synonym['ome_id'] = gene.id
                    ome_synonym['synonym'] = splitrefs[1]

                    try:
                        data_source_id = db_xref_data_source_id[splitrefs[0]]
                    except KeyError:
                        data_source = base.DataSource(name=splitrefs[0])
                        session.add(data_source)
                        session.flush()
                        data_source_id = data_source.id
                        db_xref_data_source_id[splitrefs[0]] = data_source_id

                    ome_synonym['synonym_data_source_id'] = data_source_id
                    found_synonym = (session
                                     .query(base.Synonym)
                                     .filter(base.Synonym.ome_id == gene.id)
                                     .filter(base.Synonym.synonym == splitrefs[1])
                                     .filter(base.Synonym.type == 'gene')
                                     .count() > 0)
                    if not found_synonym:
                        synonym = base.Synonym(**ome_synonym)
                        session.add(synonym)
                else:
                    print gene.id, " ->the gene id is none"

        if 'gene_synonym' in feature.qualifiers:
            for ref in feature.qualifiers['gene_synonym']:
                syn = ref.split(';')
                for value in syn:
                    ome_synonym = {'type': 'gene'}
                    ome_synonym['ome_id'] = gene.id
                    ome_synonym['synonym'] = value
                    ome_synonym['synonym_data_source_id'] = None
                    found_synonym = (session
                                     .query(base.Synonym)
                                     .filter(base.Synonym.ome_id == gene.id)
                                     .filter(base.Synonym.synonym == value)
                                     .filter(base.Synonym.type == 'gene')
                                     .count() > 0)
                    if not found_synonym:
                        synonym = base.Synonym(**ome_synonym)
                        session.add(synonym)
                        
        if 'note' in feature.qualifiers:
            for ref in feature.qualifiers['note']:
                syn = ref.split(';')
                for value in syn:
                    value  = value.split(':')
                    if value[0]== 'ORF_ID':
                        ome_synonym = {'type': 'gene'}
                        ome_synonym['ome_id'] = gene.id
                        ome_synonym['synonym'] = value[1]
                        ome_synonym['synonym_data_source_id'] = None
                        found_synonym = (session
                                         .query(base.Synonym)
                                         .filter(base.Synonym.ome_id == gene.id)
                                         .filter(base.Synonym.synonym == value[1])
                                         .filter(base.Synonym.type == 'gene')
                                         .count() > 0)
                        if not found_synonym:
                            synonym = base.Synonym(**ome_synonym)
                            session.add(synonym)

    session.commit()


@timing
def load_motifs(base, components):
    motifs = open(settings.data_directory + '/annotation/ec_annotation_from_metacyc_2010July19_wpseudo.gff','r')


@timing
def load_metacyc_proteins(base, components, genome, metacyc_proteins):
    session = base.Session()

    for unique_id,entry in metacyc_proteins.iteritems():

        vals = scrub_metacyc_entry(entry)
        if vals is None: continue

        if 'Protein-Complexes' in vals['TYPES'] or 'Protein-Small-Molecule-Complexes' in vals['TYPES']:
            get_or_create_metacyc_protein_complex(session, base, components, genome, entry)

        elif 'Polypeptides' in vals['TYPES']:
            update_protein_with_metacyc(session, base, components, genome, entry)

    session.close()


@timing
def load_metacyc_protein_cplxs(base, components, genome, metacyc_protein_cplxs):
    session = base.Session()

    for unique_id,entry in metacyc_protein_cplxs.iteritems():

        vals = scrub_metacyc_entry(entry)
        if vals is None: continue

        if 'Protein-Complexes' or 'Protein-Small-Molecule-Complexes' in vals['TYPES']:
            get_or_create_metacyc_protein_complex(session, base, components, entry)


@timing
def load_metacyc_transcription_units(base, components, genome, metacyc_tus,
                                     metacyc_genes, metacyc_promoters):
    session = base.Session()
    ##First load annotation file containing merger of metacyc and NCBI
    metacyc_ID = session.get_or_create(base.DataSource, name="metacyc").id

    for unique_id,entry in metacyc_tus.iteritems():

        get_or_create_metacyc_transcription_unit(session, base, components,
                                                 genome, entry, metacyc_genes,
                                                 metacyc_promoters)

    session.close()


@timing
def load_metacyc_bindsites(base, components, chromosome):
    session = base.Session()
    ##First load annotation file containing merger of metacyc and NCBI
    metacyc_ID = session.get_or_create(base.DataSource, name="metacyc").id
    metacyc_binding_sites = parse_metacyc_dat('dnabindsites.dat')
    metacyc_regulation = parse_metacyc_dat('regulation.dat')

    for unique_id,entry in metacyc_binding_sites.iteritems():

        vals = scrub_metacyc_entry(entry, args=['UNIQUE-ID','TYPES','ABS-CENTER-POS'])
        if vals is None: continue

        if 'DNA-Binding-Sites' in vals['TYPES']:

            try: centerpos = math.floor(float(vals['ABS-CENTER-POS'][0]))
            except: continue

            try:
                length = int(metacyc_binding_sites[unique_id]['SITE-LENGTH'][0])
                leftpos = math.floor(centerpos-(length/2))
                rightpos = math.floor(centerpos+(length/2))
            except:
                length = 0
                leftpos = centerpos
                rightpos = centerpos

            session.get_or_create(components.DnaBindingSite,
                                  bigg_id=vals['UNIQUE-ID'][0],
                                  name=vals['COMMON-NAME'][0], leftpos=leftpos,
                                  rightpos=rightpos, strand='+',
                                  chromosome_id=chromosome.id,
                                  centerpos=centerpos, width=length)


    for unique_id,entry in metacyc_regulation.iteritems():
        vals = scrub_metacyc_entry(entry, args=['UNIQUE-ID','TYPES','ASSOCIATED-BINDING-SITE','REGULATOR'])
        if vals is None: continue

        if 'Transcription-Factor-Binding' in vals['TYPES']:
            regulator = session.get_or_create(components.Component, name=vals['REGULATOR'][0])

            binding_site = session.query(components.DnaBindingSite).filter_by(name=vals['ASSOCIATED-BINDING-SITE'][0]).first()
            if binding_site is None: continue

            tf_binding_complex = session.get_or_create(components.Complex, name=vals['UNIQUE-ID'][0])

            session.get_or_create(components.ComplexComposition, complex_id=tf_binding_complex.id,\
                                                                     component_id=regulator.id,\
                                                                     stoichiometry=1.)

            session.get_or_create(components.ComplexComposition, complex_id=tf_binding_complex.id,\
                                                                     component_id=binding_site.id,\
                                                                     stoichiometry=1.)


@timing
def load_metacyc_chemicals(ome, metacyc_ID=None):
    if metacyc_ID is None:
        metacyc_ID = session.get_or_create(Dataset, name="metacyc").id
    chemicals = parseEco_dat('compounds.dat')
    markup = re.compile("<.+?>")
    for c in chemicals:
        try:
            ID = chemicals[c]['UNIQUE-ID'][0]
            name = chemicals[c]['COMMON-NAME'][0].replace('\'','[prime]')
            smiles = chemicals[c]['SMILES'][0]
        except:
            continue
        name = markup.sub("", name)
        result = ome.execute("INSERT INTO chemicals(name, smiles) VALUES " + \
                        "('%s', '%s') RETURNING id;" % (name, smiles))
        id = result.fetchone()[0]

        ome.execute("INSERT INTO id2otherID(id, otherID, type, dataset_id) VALUES " + \
                        "(%i, '%s', '%s', %i);"%(id, ID, 'chemical', metacyc_ID))
        make_citations(ome, chemicals[c], id)
        ome.commit()


@timing
def load_metacyc_promoters(ome, metacyc_ID=None):
    if metacyc_ID is None:
        metacyc_ID = session.get_or_create(Dataset, name="metacyc").id
    promoters = parseEco_dat('promoters.dat')

    for p in promoters:
        try:
            ID = promoters[p]['UNIQUE-ID'][0]
            position = int(promoters[p]['ABSOLUTE-PLUS-1-POS'][0])
        except: continue
        try: name = promoters[p]['COMMON-NAME'][0]
        except: name = ''

        result = ome.execute("INSERT INTO tss(name, position) VALUES ('%s', %i) RETURNING id;"%(name, position))
        id = result.fetchone()[0]

        ome.execute("INSERT INTO id2otherID(id, otherID, type, dataset_id) VALUES " + \
                        "(%i, '%s', '%s', %i);"%(id, ID, 'tss', metacyc_ID))
        make_citations(ome, promoters[p], id)
    ome.commit()


@timing
def load_kegg_pathways(base, components):

    session = base.Session()
    kegg_pathways = open(settings.data_directory+'/annotation/KEGG/ecoli_KEGG_pathways.txt','r')
    kegg_dataset = session.get_or_create(base.DataSource, name="kegg")
    for line in kegg_pathways.readlines():
        vals = line.rstrip('\r\n').split('\t')
        bnums = vals[2].split(',')
        pathway_ID = vals[0]
        pathway_name = vals[1].strip()
        kegg_pathway = session.get_or_create(components.GeneGroup, name = pathway_name)
        for bnum in bnums:
            gene = session.query(components.Gene).filter_by(bigg_id=bnum).first()
            if not gene: continue
            session.get_or_create(components.GeneGrouping, group_id = kegg_pathway.id, gene_id = gene.id)
    session.commit()
    session.close()


@timing
def load_regulatory_network(base, components, datasets, chromosome):

    datasets.RegulatoryNetwork.__table__.drop()
    datasets.RegulatoryNetwork.__table__.create()

    session = base.Session()
    regulatory_network = open(settings.data_directory+'/annotation/regulondb_gene_reg_network.txt','r')
    for line in regulatory_network.readlines():
        vals = line.split('\t')


        reg_gene = session.query(components.Gene).filter(and_(func.lower(components.Gene.name) == vals[0].lower(),
                                                              components.Gene.chromosome_id == chromosome.id)).first()

        regd_gene = session.query(components.Gene).filter(and_(func.lower(components.Gene.name) == vals[1].lower(),
                                                               components.Gene.chromosome_id == chromosome.id)).first()

        if reg_gene is None or regd_gene is None: continue

        session.get_or_create(datasets.RegulatoryNetwork, reg_gene_id=reg_gene.id, regd_gene_id=regd_gene.id,
                                                          direction=vals[2],
                                                          evidence=vals[3])

    session.flush()
    session.commit()
    session.close()


@timing
def write_chromosome_annotation_gff(base, components, chromosome):
    session = base.Session()

    genbank_fasta_string = 'gi|'+chromosome.genbank_id+'|ref|'+chromosome.ncbi_id+'|'

    with open(settings.data_directory+'/annotation/'+chromosome.ncbi_id+'.gff', 'wb') as gff_file:

        for gene in session.query(components.Gene).filter(components.Gene.chromosome_id == chromosome.id).all():

            info_string = 'gene_id "%s"; transcript_id "%s"; gene_name "%s";' % (gene.bigg_id, gene.bigg_id, gene.name)

            gff_string = '%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s\n' % (genbank_fasta_string, 'ome_db', 'exon', gene.leftpos,
                                                                                                         gene.rightpos,
                                                                                                         gene.strand,
                                                                                                         info_string)
            gff_file.write(gff_string)

    session.close()

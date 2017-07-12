# -*- coding: utf-8 -*-

from cobradb.models import *
from cobradb.component_loading import get_genbank_accessions

import pytest


def test_get_genbank_accessions(test_genbank_files):
    acc = get_genbank_accessions(test_genbank_files[0][1], fast=True)
    assert acc == {'ncbi_accession': 'NC_000913.2', 'ncbi_assembly': None, 'ncbi_bioproject': None}
    acc = get_genbank_accessions(test_genbank_files[0][1], fast=False)
    assert acc == {'ncbi_accession': 'NC_000913.2', 'ncbi_assembly': None, 'ncbi_bioproject': None}
    acc = get_genbank_accessions(test_genbank_files[1][1], fast=True)
    assert acc == {'ncbi_accession': 'NC_000913.2',
                   'ncbi_assembly': 'test_assembly.1',
                   'ncbi_bioproject': 'PRJNA57779-core-2'}
    acc = get_genbank_accessions(test_genbank_files[1][1], fast=False)
    assert acc == {'ncbi_accession': 'NC_000913.2',
                   'ncbi_assembly': 'test_assembly.1',
                   'ncbi_bioproject': 'PRJNA57779-core-2'}


@pytest.mark.usefixtures('load_genomes')
class TestWithGenomes():
    def test_genome_taxon(self, session):
        assert session.query(Genome).first().taxon_id == '511145'

    def test_genome_genes(self, session):
        assert (session.query(Gene)
                .filter(Gene.bigg_id == 'b0114')
                .count()) == 2

    def test_genome_synonyms_locus_tag(self, session):
        assert (session.query(Synonym)
                .join(DataSource)
                .join(Gene, Gene.id == Synonym.ome_id)
                .filter(Synonym.type == 'gene')
                .filter(DataSource.bigg_id == 'refseq_locus_tag')
                .filter(Synonym.synonym == 'b0114')
                .filter(Gene.bigg_id == 'b0114')
                .count()) == 2

    def test_genome_empty_name(self, session):
        gene_db = (session.query(Gene)
                   .filter(Gene.bigg_id == 'b0720')
                   .first())
        assert gene_db.name is None

    def test_gene_name_no_locus_id(self, session):
        """Keep the gene name when there is no locus ID, so that, if alternative
        transcripts are used later to match agains a synonym, the gene name is still
        around.

        """
        assert session.query(Gene).filter(Gene.bigg_id == 'appB').first().name == 'appB'

    def test_gene_name_no_locus_id_gene_id(self, session):
        """Use the Gene ID (ncbigene) if present."""
        assert session.query(Gene).filter(Gene.bigg_id == '945585').first().name == 'appC'

    def test_genome_synonyms_name(self, session):
        assert (session.query(Synonym)
                .join(DataSource)
                .filter(DataSource.bigg_id == 'refseq_name')
                .filter(Synonym.synonym == 'aceE')
                .count()) == 2

    def test_genome_synonyms_synonyms(self, session):
        assert (session.query(Synonym)
                .join(DataSource)
                .filter(DataSource.bigg_id == 'refseq_synonym')
                .filter(Synonym.synonym == 'ECK0113')
                .count()) == 2


    def test_genome_synonyms_db_xref(self, session):
        assert (session.query(Synonym)
                .join(DataSource)
                .filter(DataSource.bigg_id == 'ncbigi')
                .filter(Synonym.synonym == '16128107')
                .count()) == 2


    def test_genome_synonyms_db_xref_duplicate(self, session):
        # this causes an error when we are not dealing with duplicates correctly
        assert (session.query(Synonym)
                .join(DataSource)
                .filter(DataSource.bigg_id == 'tests_dup_syn')
                .filter(Synonym.synonym == 'b0114')
                .count()) == 1


    # only in core.gb:
    def test_genome_synonyms_old_locus_tag(self, session):
        assert (session.query(Synonym)
                .join(DataSource)
                .filter(DataSource.bigg_id == 'refseq_old_locus_tag')
                .filter(Synonym.synonym == 'test_b0114')
                .count()) == 1


    # only in core.gb:
    def test_genome_synonyms_orf_id(self, session):
        assert (session.query(Synonym)
                .join(DataSource)
                .filter(DataSource.bigg_id == 'refseq_orf_id')
                .filter(Synonym.synonym == 'test_orf')
                .count()) == 1

    def test_dna_sequence(self, session):
        assert (session.query(Gene.dna_sequence)
                .join(Chromosome)
                .join(Genome)
                .filter(Gene.bigg_id == 'b0008')
                .filter(Genome.accession_value == 'core')
                .one())[0] == 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACCTCGAATCTACCGTCGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGGCTGATCACATGGTGCTGATGGCAGGTTTCACCGCCGGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGT'

    def test_protein_sequence(self, session):
        assert (session.query(Gene.protein_sequence)
                .join(Chromosome)
                .join(Genome)
                .filter(Gene.bigg_id == 'b0008')
                .filter(Genome.accession_value == 'core')
                .one())[0] == 'MTDKLTSLRQYTTVVADTGDIAAMKLYQPQDATTNPSLILNAAQIPEYRKLIDDAVAWAKQQSNDRAQQIVDATDKLAVNIGLEILKLVPGRISTEVDARLSYDTEASIAKAKRLIKLYNDAGISNDRILIKLASTWQGIRAAEQLEKEGINCNLTLLFSFAQARACAEAGVFLISPFVGRILDWYKANTDKKEYAPAEDPGVVSVSEIYQYYKEHGYETVVMGASFRNIGEILELAGCDRLTIAPALLKELAESEGAIERKLSYTGEVKARPARITESEFLWQHNQDPMAVDKLAEGIRKFAIDQEKLEKMIGDLL'

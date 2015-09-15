# -*- coding: utf-8 -*-

from ome.util import *
from ome.base import *

import pytest


def test_increment_id():
    assert increment_id('ACALD_1') == 'ACALD_2'
    assert increment_id('ACALD_1a') == 'ACALD_1a_1'
    assert increment_id('ACALD') == 'ACALD_1'
    assert increment_id('ACALD_9') == 'ACALD_10'
    assert increment_id('ACALD_10') == 'ACALD_11'
    # name
    assert increment_id('ACALD_1', 'copy') == 'ACALD_1_copy1'
    assert increment_id('ACALD_copy1', 'copy') == 'ACALD_copy2'


def test_check_pseudoreaction():
    assert check_pseudoreaction('ATPM') is True
    assert check_pseudoreaction('ATPM_NGAM') is True
    assert check_pseudoreaction('ATPM1') is False
    assert check_pseudoreaction('EX_glc_e') is True
    assert check_pseudoreaction('aEX_glc_e') is False
    assert check_pseudoreaction('biomass_objective') is True
    assert check_pseudoreaction('BiomassEcoli') is True
    assert check_pseudoreaction('DM_8') is True


def test_find_data_source_url():
    url_prefs = [['KEGGID', 'http://identifiers.org/kegg.compound/']]
    assert find_data_source_url('KEGGID', url_prefs) == 'http://identifiers.org/kegg.compound/'


def test_read_data_source_preferences(test_prefs):
    settings.data_source_preferences = test_prefs['data_source_preferences']
    urls_dictionary = read_data_source_preferences()
    assert urls_dictionary[0][0].split()[0] == 'KEGGID'
    assert urls_dictionary[0][0].split()[1] == 'http://identifiers.org/kegg.compound/'


def test_get_or_create_data_source(test_db, session):
    get_or_create_data_source(session, 'my_data_source')
    assert (session
            .query(DataSource)
            .filter(DataSource.name == 'my_data_source')
            .count()) == 1


def test_scrub_gene_id():
    assert scrub_gene_id('1234.5') == '1234_AT5'
    assert scrub_gene_id('1234.56') == '1234_AT56'
    assert scrub_gene_id('1234.56a') == '1234_56a'
    assert scrub_gene_id('asdkf@#%*(@#$sadf') == 'asdkf________sadf'


def test_load_tsv(tmpdir):
    # test file
    a_file = tmpdir.join('temp.txt')
    a_file.write('# ignore\tignore\na\ttest  \n\n')
    # run the test
    rows = load_tsv(str(a_file))
    assert rows == [['a', 'test']]

    # with required_column_num
    rows = load_tsv(str(a_file), required_column_num=3)
    assert rows == []

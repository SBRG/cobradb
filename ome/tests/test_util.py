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


def test_make_reaction_copy_id():
    assert make_reaction_copy_id('ACALD', 3) == 'ACALD_copy3'


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


def test_get_or_create_data_source(test_db, session, test_prefs, tmpdir):
    prefsfile = str(tmpdir.join('data_source_preferences.txt'))
    with open(prefsfile, 'w') as f:
        f.write('my_data_source\tmy_url_prefix')

    settings.data_source_preferences = prefsfile

    get_or_create_data_source(session, 'my_data_source')
    assert (session
            .query(DataSource)
            .filter(DataSource.name == 'my_data_source')
            .filter(DataSource.url_prefix == 'my_url_prefix')
            .count()) == 1


def test_scrub_gene_id():
    assert scrub_gene_id('1234.5') == '1234_AT5'
    assert scrub_gene_id('1234.56') == '1234_AT56'
    assert scrub_gene_id('1234.56a') == '1234_56a'
    assert scrub_gene_id('asdkf@#%*(@#$sadf') == 'asdkf________sadf'


def test_scrub_name():
    assert scrub_name('R_ammonia_reversible_transport') == 'Ammonia reversible transport'
    assert scrub_name('_ammonia_reversible_transport') == 'Ammonia reversible transport'
    assert scrub_name(None) == None
    assert scrub_name('_ ') == None


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

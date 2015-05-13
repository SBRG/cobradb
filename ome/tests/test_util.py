from ome.util import increment_id, check_pseudoreaction

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
    assert check_pseudoreaction('ATPM(NGAM)') is True
    assert check_pseudoreaction('ATPM1') is False
    assert check_pseudoreaction('EX_glc_e') is True
    assert check_pseudoreaction('aEX_glc_e') is False
    assert check_pseudoreaction('biomass_objective') is True
    assert check_pseudoreaction('BiomassEcoli') is True
    assert check_pseudoreaction('DM_8') is True

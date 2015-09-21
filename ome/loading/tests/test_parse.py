from ome.loading.parse import *

from cobra.core import Reaction, Metabolite
from cobra.io import read_sbml_model
import pytest


@pytest.fixture(scope='session')
def example_model(test_model_files):
    return read_sbml_model(test_model_files[0]['path'])


@pytest.fixture(scope='session')
def convert_ids_model(example_model):
    example_model.id = 'A bad id'
    example_model.add_reaction(Reaction('DADA'))
    example_model.reactions.get_by_id('DADA').add_metabolites({
        Metabolite('dad_DASH_2_c'): -1
    })
    return convert_ids(example_model.copy())


def test_convert_ids_dad_2(convert_ids_model):
    returned, old_ids = convert_ids_model
    assert returned.id == 'A_bad_id'
    assert 'dad_2_c' in returned.metabolites
    assert 'dad_2_c' in [x.id for x in returned.reactions.get_by_id('DADA').metabolites]
    assert ('dad_2_c', ['dad_DASH_2_c']) in old_ids['metabolites'].items()


def test_convert_ids_repeats_reactions(convert_ids_model):
    returned, old_ids = convert_ids_model
    assert 'EX_glu__L_e' in returned.reactions
    assert 'EX_glu__L_e_1' in returned.reactions
    assert 'EX_glu__L_e_2' in returned.reactions
    assert 'EX_gln__L_e' in returned.reactions
    assert 'EX_gln__L_e_1' in returned.reactions

    old_ids_list = old_ids['reactions'].items()
    assert (
        ('EX_gln__L_e', ['EX_gln_L_e']) in old_ids_list and
        ('EX_gln__L_e_1', ['EX_gln__L_e']) in old_ids_list
    ) or (
        ('EX_gln__L_e', ['EX_gln__L_e']) in old_ids_list and
        ('EX_gln__L_e_1', ['EX_gln_L_e']) in old_ids_list
    )


def test_convert_ids_repeats_metabolites(convert_ids_model):
    # metabolites should get merged
    returned, old_ids = convert_ids_model
    assert set(old_ids['metabolites']['glc__D_e']) == {'glc_D_e', 'glc_DASH_D_e'}


def test_convert_ids_repeats_gene(convert_ids_model):
    # gene should get merged
    returned, old_ids = convert_ids_model
    assert '904_AT1' in [x.id for x in returned.reactions.get_by_id('FORt2').genes]
    assert '904_AT1' in [x.id for x in returned.reactions.get_by_id('FORti').genes]
    assert set(old_ids['genes']['904_AT1']) == {'904.1', '904_AT1'}


def test_convert_ids_genes(convert_ids_model, example_model):
    returned, old_ids = convert_ids_model

    # lost a gene to merging
    assert len(returned.genes) == len(example_model.genes) - 1

    assert 'gene_with_period_AT22' in [x.id for x in returned.genes]
    assert (
        returned.reactions.get_by_id('FRD7').gene_reaction_rule ==
        returned.reactions.get_by_id('FRD7').gene_reaction_rule
        .replace('.22', '_AT22').replace('.12', '_AT12')
    )
    assert old_ids['genes']['gene_with_period_AT22'] == ['gene_with_period.22']

    assert ['.22' not in x.id for x in returned.genes]
    assert ['.22' not in x.gene_reaction_rule for x in returned.reactions]


def test_id_for_new_id_style():
    """Test edge cases for the ID conversion."""
    cases = [x.strip() for x in """

M_sertrna_sec__c
M_lipidA_core_e_p
M_lipa_cold_e
M_lipa_cold_p
M_lipa_cold_c
M_sertrna_sec__c
M_lipidA_core_e_p
M_lipidA_core_e_p

    """.split('\n') if x.strip() != '']

    for case in cases:
        new = id_for_new_id_style(case, is_metabolite=True)
        met, compartment = split_compartment(new)

    # strip leading underscores
    assert id_for_new_id_style('_13dpg_c') == '13dpg_c'
    assert id_for_new_id_style('__13dpg_c') == '13dpg_c'

    # 2 character compartment
    assert id_for_new_id_style('abc(c1)') == 'abc_c1'

    # remove internal __
    assert id_for_new_id_style('26dap__Z_c') == '26dap_Z_c'
    assert id_for_new_id_style('26dap_Z_c') == '26dap_Z_c'
    # except with [LDSRM]
    assert id_for_new_id_style('26dap__M_c') == '26dap__M_c'
    assert id_for_new_id_style('26dap__M_c') == '26dap__M_c'

    # other characters
    assert id_for_new_id_style('ATPM(NGAM)') == 'ATPM_NGAM'
    assert id_for_new_id_style('a()[]c*&^%b') == 'a_c_b'


def test_hash_reaction(test_model_files):
    # there are no conflicts in model 2
    model, _ = load_and_normalize(test_model_files[1]['path'])

    # just the string
    string = hash_reaction(model.reactions.GAPD, string_only=True)
    assert string == '13dpg_c1.000g3p_c-1.000h_c1.000nad_c-1.000nadh_c1.000pi_c-1.000'

    # no conflicts
    num = 20
    hashes = {r.id: hash_reaction(r) for r in model.reactions[:20]}
    assert len(set(hashes.values())) == 20

    # repeatable
    k1, h1 = hashes.iteritems().next()
    assert h1 == hash_reaction(model.reactions.get_by_id(k1))

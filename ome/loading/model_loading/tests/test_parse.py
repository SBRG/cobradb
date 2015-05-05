from ome.loading.model_loading.parse import *

from cobra.core import Reaction, Metabolite
from cobra.io import read_sbml_model

def test_convert_ids(test_model):
    model_in = read_sbml_model(test_model[0]['path'])

    model_in.add_reaction(Reaction('DADA'))
    model_in.reactions.get_by_id('DADA').add_metabolites({
        Metabolite('dad_DASH_2_c'): -1
    })
    returned, old_ids = convert_ids(model_in, 'cobrapy')

    assert 'dad__2_c' in returned.metabolites
    assert 'dad__2_c' in [x.id for x in returned.reactions.get_by_id('DADA').metabolites]
    assert ('dad__2_c', 'dad_DASH_2_c') in old_ids['metabolites'].items()
    assert ('EX_gln__L_e', 'EX_gln_L_e') in old_ids['reactions'].items()


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
    assert id_for_new_id_style('_13dpg') == '13dpg'


def test_hash_reaction(test_model):
    # there are no conflicts in model 2
    model, _ = load_and_normalize(test_model[1]['path'])

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

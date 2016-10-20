from cobradb.loading.parse import *
from cobradb.loading.parse import (_has_gene_reaction_rule,
                                   _normalize_pseudoreaction)

from cobra.core import Reaction, Metabolite, Model
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


# --------------------------------------------------------------------
# pseudoreactions
# --------------------------------------------------------------------

def test__has_gene_reaction_rule():
    reaction = Reaction('rxn')
    assert _has_gene_reaction_rule(reaction) is False
    reaction.gene_reaction_rule = 'b1779'
    assert _has_gene_reaction_rule(reaction) is True
    reaction.gene_reaction_rule = ' '
    assert _has_gene_reaction_rule(reaction) is False


def test__normalize_pseudoreaction_exchange():
    reaction = Reaction('EX_gone')
    reaction.add_metabolites({Metabolite('glu__L_e'): -1})
    reaction.lower_bound = -1000
    reaction.upper_bound = 0
    _normalize_pseudoreaction(reaction)
    assert reaction.id == 'EX_glu__L_e'
    assert reaction.subsystem == 'Extracellular exchange'


def test__normalize_pseudoreaction_exchange_reversed():
    reaction = Reaction('EX_gone')
    reaction.add_metabolites({Metabolite('glu__L_e'): 1})
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    _normalize_pseudoreaction(reaction)
    assert reaction.id == 'EX_glu__L_e'
    assert reaction.lower_bound == -1000
    assert reaction.upper_bound == 0
    assert reaction.metabolites.values() == [-1]


def test__normalize_pseudoreaction_exchange_error_bad_coeff():
    reaction = Reaction('EX_gone')
    reaction.add_metabolites({Metabolite('glu__L_e'): -2})
    with pytest.raises(ConflictingPseudoreaction) as excinfo:
        _normalize_pseudoreaction(reaction)
    assert 'with coefficient' in str(excinfo.value)
    assert reaction.id == 'EX_gone'


def test__normalize_pseudoreaction_exchange_error_bad_name():
    reaction = Reaction('gone')
    reaction.add_metabolites({Metabolite('glu__L_e'): -1})
    _normalize_pseudoreaction(reaction)
    assert reaction.id == 'EX_glu__L_e'


def test__normalize_pseudoreaction_exchange_error_has_gpr():
    reaction = Reaction('EX_gone')
    reaction.add_metabolites({Metabolite('glu__L_e'): -1})
    reaction.gene_reaction_rule = 'b1779'
    with pytest.raises(ConflictingPseudoreaction) as excinfo:
        _normalize_pseudoreaction(reaction)
    assert 'has a gene_reaction_rule' in str(excinfo.value)
    assert reaction.id == 'EX_gone'


def test__normalize_pseudoreaction_demand():
    reaction = Reaction('DM_gone')
    reaction.add_metabolites({Metabolite('glu__L_c'): -1})
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    _normalize_pseudoreaction(reaction)
    assert reaction.id == 'DM_glu__L_c'
    assert reaction.subsystem == 'Intracellular demand'


def test__normalize_pseudoreaction_demand_reversed():
    reaction = Reaction('DM_gone')
    reaction.add_metabolites({Metabolite('glu__L_c'): 1})
    reaction.lower_bound = -1000
    reaction.upper_bound = 0
    _normalize_pseudoreaction(reaction)
    assert reaction.metabolites.values() == [-1]
    assert reaction.lower_bound == 0
    assert reaction.upper_bound == 1000
    assert reaction.id == 'DM_glu__L_c'


def test__normalize_pseudoreaction_demand_reversed_prefer_sink_name():
    reaction = Reaction('sink_gone')
    reaction.add_metabolites({Metabolite('glu__L_c'): 1})
    reaction.lower_bound = -1000
    reaction.upper_bound = 0
    _normalize_pseudoreaction(reaction)
    assert reaction.metabolites.values() == [-1]
    assert reaction.lower_bound == 0
    assert reaction.upper_bound == 1000
    assert reaction.id == 'SK_glu__L_c'


def test__normalize_pseudoreaction_demand_error_has_gpr():
    reaction = Reaction('DM_gone')
    reaction.add_metabolites({Metabolite('glu__L_c'): -1})
    reaction.gene_reaction_rule = 'b1779'
    with pytest.raises(ConflictingPseudoreaction) as excinfo:
        _normalize_pseudoreaction(reaction)
    assert 'has a gene_reaction_rule' in str(excinfo.value)
    assert reaction.id == 'DM_gone'


def test__normalize_pseudoreaction_sink():
    reaction = Reaction('SInk_gone')
    reaction.add_metabolites({Metabolite('glu__L_c'): -1})
    reaction.lower_bound = -1000
    reaction.upper_bound = 0
    _normalize_pseudoreaction(reaction)
    assert reaction.id == 'SK_glu__L_c'
    assert reaction.subsystem == 'Intracellular source/sink'


def test__normalize_pseudoreaction_sink_reversed():
    reaction = Reaction('Sink_gone')
    reaction.add_metabolites({Metabolite('glu__L_c'): 1})
    reaction.lower_bound = 0
    reaction.upper_bound = 50
    _normalize_pseudoreaction(reaction)
    assert reaction.metabolites.values() == [-1]
    assert reaction.lower_bound == -50
    assert reaction.upper_bound == 0
    assert reaction.id == 'SK_glu__L_c'


def test__normalize_pseudoreaction_biomass():
    reaction = Reaction('my_biomass_2')
    _normalize_pseudoreaction(reaction)
    assert reaction.id == 'BIOMASS_my_2'
    assert reaction.subsystem == 'Biomass and maintenance functions'


def test__normalize_pseudoreaction_biomass_has_gpr():
    reaction = Reaction('my_biomass_2')
    reaction.gene_reaction_rule = 'b1779'
    with pytest.raises(ConflictingPseudoreaction) as excinfo:
        _normalize_pseudoreaction(reaction)
    assert 'has a gene_reaction_rule' in str(excinfo.value)
    assert reaction.id == 'my_biomass_2'


def test__normalize_pseudoreaction_atpm():
    reaction = Reaction('notATPM')
    reaction.add_metabolites({Metabolite('atp_c'): -1,
                              Metabolite('h2o_c'): -1,
                              Metabolite('pi_c'): 1,
                              Metabolite('h_c'): 1,
                              Metabolite('adp_c'): 1})
    _normalize_pseudoreaction(reaction)
    assert reaction.id == 'ATPM'
    assert reaction.subsystem == 'Biomass and maintenance functions'


def test__normalize_pseudoreaction_atpm_reversed():
    reaction = Reaction('notATPM')
    reaction.add_metabolites({Metabolite('atp_c'): 1,
                              Metabolite('h2o_c'): 1,
                              Metabolite('pi_c'): -1,
                              Metabolite('h_c'): -1,
                              Metabolite('adp_c'): -1})
    reaction.lower_bound = -50
    reaction.upper_bound = 100
    _normalize_pseudoreaction(reaction)
    assert reaction.id == 'ATPM'
    assert reaction.lower_bound == -100
    assert reaction.upper_bound == 50


def test__normalize_pseudoreaction_atpm_has_gpr():
    reaction = Reaction('NPT1')
    reaction.add_metabolites({Metabolite('atp_c'): -1,
                              Metabolite('h2o_c'): -1,
                              Metabolite('pi_c'): 1,
                              Metabolite('h_c'): 1,
                              Metabolite('adp_c'): 1})
    reaction.gene_reaction_rule = 'b1779'
    _normalize_pseudoreaction(reaction)
    # should not change
    assert reaction.id == 'NPT1'


# --------------------------------------------------------------------
# ID fixes
# --------------------------------------------------------------------

def test_convert_ids_dad_2(convert_ids_model):
    returned, old_ids = convert_ids_model
    assert returned.id == 'A_bad_id'
    assert 'dad_2_c' in returned.metabolites
    assert 'dad_2_c' in [x.id for x in returned.reactions.get_by_id('DM_dad_2_c').metabolites]
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

def test_custom_hashes():
    # These hashes are generated from old IDs in models (in reaction strings),
    # and they match to these corrected BiGG reaction IDs
    cases = [
        ('39b5f90a1919aef07473e2f835ce63af', 'EX_frmd_e', 'foam_e <=>'),
        ('92f1047c72db0a36413d822863be514e', 'EX_phllqne_e', 'phyQ_e <=>'),
    ]
    model = Model()
    for reaction_hash, bigg_id, reaction_string in cases:
        reaction = Reaction(bigg_id)
        model.add_reaction(reaction)
        reaction.build_reaction_from_string(reaction_string)
        assert hash_reaction(reaction) == reaction_hash

def test_reverse_reaction():
    model = Model()
    reaction = Reaction('AB')
    model.add_reaction(reaction)
    reaction.build_reaction_from_string('a -> b')
    reversed_reaction = reverse_reaction(reaction)
    assert reversed_reaction.reaction == 'b --> a'

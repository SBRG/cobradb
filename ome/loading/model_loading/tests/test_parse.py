from ome.loading.model_loading.parse import *

def test_hash_reaction(test_model):
    model, _ = load_and_normalize(test_model[0]['path'])

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

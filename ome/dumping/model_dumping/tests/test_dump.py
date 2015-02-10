from ome.dumping.model_dumping import dump_model
from ome.loading import load_model, load_genome

import pytest
import logging

def test_dump_model(test_genbank, test_model, test_db, setup_logger):
    logging.info('Testing model dump')
    
    timestamp = '2014-9-16 14:26:22'
    pmid = '25575024'
    # load the test genome
    
    load_genome(test_genbank[0]['path'])

    # load the model
    load_model(test_model[0]['id'], test_model[0]['dir'], test_genbank[0]['genome_id'],
               timestamp, pmid)
    
    with pytest.raises(Exception):
        dump_model('C3PO')
        
    model = dump_model(test_model[0]['bigg_id'])
    
    assert len(model.reactions) == 95
    assert len(model.metabolites) == 72
    assert len(model.genes) == 137
    # TODO test gene names
    # assert model.genes.get_by_id('b0008').name != 'b0008'
    
    assert 'GAPD' in model.reactions
    assert model.reactions.get_by_id('GAPD').name == 'glyceraldehyde-3-phosphate dehydrogenase'
    assert model.genes.get_by_id('b3528').name == 'dctA'
    # test solve
    model.reactions.get_by_id('EX_glc_e').lower_bound = -10
    assert model.optimize().f > 0

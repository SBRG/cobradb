from ome.dumping.model_dumping import dump_model
from ome.loading import load_model, load_genome
from ome import base

import pytest
import logging
import cobra.io
from os.path import join, exists
import shutil

def test_dump_model(test_genbank, test_model, test_db, setup_logger, tmpdir):
    logging.info('Testing model dump')

    session = base.Session()
    
    timestamp = '2014-9-16 14:26:22'
    pmid = '25575024'
    # load the test genome
    
    load_genome(test_genbank[0]['path'], session)

    # load the model
    bigg_id = load_model(test_model[0]['path'], test_genbank[0]['genome_id'],
               timestamp, pmid, session, dump_directory=str(tmpdir))
    assert exists(join(str(tmpdir), bigg_id+'.xml'))
    
    with pytest.raises(Exception):
        dump_model('C3PO', session)
        
    model = dump_model(bigg_id)
    
    assert len(model.reactions) == 95
    assert len(model.metabolites) == 72
    assert len(model.genes) == 137
    assert model.genes.get_by_id('b0114').name == 'aceE'
    assert model.genes.get_by_id('b3528').name == 'dctA'
    
    assert 'GAPD' in model.reactions
    assert model.reactions.get_by_id('GAPD').name == 'glyceraldehyde-3-phosphate dehydrogenase'

    # test solve
    model.reactions.get_by_id('EX_glc_e').lower_bound = -10
    assert model.optimize().f > 0

    # temp sbml dump
    m_file = join(str(tmpdir), 'test_model.sbml')
    cobra.io.write_sbml_model(model, m_file)
    loaded_model = cobra.io.read_sbml_model(m_file)

    session.close()
    
    # delete test directory
    shutil.rmtree(str(tmpdir))

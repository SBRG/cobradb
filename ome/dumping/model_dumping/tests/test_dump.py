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
    pub_ref = 'pmid:25575024'
    # load the test genome
    
    load_genome(test_genbank[0]['path'], session)

    # load the model
    bigg_id = load_model(test_model[0]['path'], test_genbank[0]['genome_id'],
                         timestamp, pub_ref, session, dump_directory=str(tmpdir),
                         published_directory=str(tmpdir),
                         polished_directory=None)
    assert exists(join(str(tmpdir), bigg_id+'.xml'))
    
    with pytest.raises(Exception):
        dump_model('C3PO', session)
        
    model = dump_model(bigg_id)
    # COBRApy uses the description as the ID sometimes. See https://github.com/opencobra/cobrapy/pull/152
    assert model.id == bigg_id
    assert model.description == bigg_id
    
    assert len(model.reactions) == 96
    assert len(model.metabolites) == 73
    assert len(model.genes) == 140
    assert model.genes.get_by_id('b0114').name == 'aceE'
    assert model.genes.get_by_id('b3528').name == 'dctA'
    
    assert 'GAPD' in model.reactions
    assert model.reactions.get_by_id('GAPD').name == 'glyceraldehyde-3-phosphate dehydrogenase'

    # make sure ATPM_NGAM and NTP1 are represented
    r1 = model.reactions.get_by_id('ATPM_NGAM')
    r2 = model.reactions.get_by_id('NTP1')
    # either order is OK
    r1 = model.reactions.get_by_id('NTP1')
    r2 = model.reactions.get_by_id('ATPM_NGAM')
    assert r1.lower_bound == 3.15
    assert r2.lower_bound == 8.39
    assert r1.notes['original_bigg_id'] == 'NTP1'
    assert r2.notes['original_bigg_id'] == 'ATPM(NGAM)'

    # test solve
    model.reactions.get_by_id('EX_glc_e').lower_bound = -10
    sol = model.optimize()
    assert sol.f > 0.5 and sol.f < 1.0

    # temp sbml dump
    m_file = join(str(tmpdir), 'test_model.sbml')
    cobra.io.write_sbml_model(model, m_file)
    loaded_model = cobra.io.read_sbml_model(m_file)
    # check model id
    assert loaded_model.id == model.id
    assert loaded_model.description == model.description

    # Check for buggy metabolite. This was being saved without an ID in the SBML
    # model for some reason.
    assert 'formmfr_b_c' in [x.id for x in model.metabolites]
    assert 'formmfr_b_c' in [x.id for x in loaded_model.metabolites]
    assert '' not in [x.id for x in loaded_model.metabolites]

    session.close()
    
    # delete test directory
    shutil.rmtree(str(tmpdir))

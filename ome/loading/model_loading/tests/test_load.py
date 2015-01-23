# -*- coding: utf-8 -*-

from ome.loading.model_loading import load_model
from ome.loading.component_loading import load_genome
from ome import base
from ome.models import *
from ome.components import *

import pytest
import os

@pytest.mark.usefixtures('test_db', 'setup_logger')
def test_load_model(test_genbank, test_model):
    timestamp = '2014-9-16 14:26:22'
    pmid = '25575024'

    # can't load the model without the genome
    with pytest.raises(Exception):
        load_model(test_model, test_genbank['genome_id'], timestamp, pmid)

    # load the test genome
    load_genome(test_genbank['path'])

    # load the model
    load_model(test_model['id'], test_model['dir'], test_genbank['genome_id'],
               timestamp, pmid)
    
    # test the model
    session = base.Session()
    assert session.query(Model).count() == 1
    assert session.query(Reaction).count() == 95
    assert session.query(ModelReaction).count() == 95
    assert session.query(CompartmentalizedComponent).count() == 72
    assert session.query(ModelCompartmentalizedComponent).count() == 72
    assert session.query(Metabolite).count() == 54
    assert session.query(Gene).count() == 137
    assert session.query(ModelGene).count() == 137
    session.close()

    # can't load the same model twice
    with pytest.raises(Exception):
        load_model(test_model, genome_id, timestamp, pmid)

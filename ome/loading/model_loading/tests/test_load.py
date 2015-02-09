# -*- coding: utf-8 -*-

from ome.loading.model_loading import load_model
from ome.loading.component_loading import load_genome
from ome import base
from ome.models import *
from ome.components import *

import pytest
import os

def test_load_model(test_genbank, test_model, test_db, setup_logger):
    session = base.Session()

    timestamp = '2014-9-16 14:26:22'
    pmid = '25575024'
    models_count = 2
    # can't load the model without the genome
    with pytest.raises(Exception):
        for x in range(0, models_count):
            load_model(test_model[x], test_genbank[x]['genome_id'], timestamp,
                       pmid, session)
    
    for x in range(0, models_count):
    # load the test genome
        load_genome(test_genbank[x]['path'], session)

    # load the model
        load_model(test_model[x]['id'], test_model[x]['dir'], test_genbank[x]['genome_id'],
                    timestamp, pmid, session)
    
    # test the model
    assert session.query(Model).count() == 2
    assert session.query(Reaction).count() == 95
    assert session.query(ModelReaction).count() == 95 * 2 
    assert session.query(CompartmentalizedComponent).count() == 72
    assert session.query(ModelCompartmentalizedComponent).count() == 72 * 2
    assert session.query(Metabolite).count() == 54
    assert session.query(LinkOut).count() == 16
    assert session.query(Gene).count() == 137 * 2
    assert session.query(ModelGene).count() == 137 * 2
    
    # test linkouts
    result = (session
              .query(LinkOut.external_source, LinkOut.external_id, LinkOut.ome_id)
              .join(Metabolite, Metabolite.id == LinkOut.ome_id)
              .filter(Metabolite.bigg_id == '13dpg')
              .filter(LinkOut.external_source == 'KEGGID')
              .all())
    assert len(result) == 2

    r_db =  (session.query(ModelReaction)
             .join(Reaction)
             .filter(Reaction.bigg_id == 'GAPD')
             .first())
    assert r_db.objective_coefficient == 0
    assert r_db.upper_bound == 1000

    # can't load the same model twice
    with pytest.raises(Exception):
        load_model(test_model, genome_id, timestamp, pmid, session)

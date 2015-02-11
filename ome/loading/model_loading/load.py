# -*- coding: utf-8 -*-

from ome import base, settings, components, timing
from ome.models import Model
from ome.loading.model_loading import independent, dependent, parse

import os
import logging

def get_model_list():
    """Get the models that are available, as SBML, in ome_data/models"""
    return [x.replace('.xml', '').replace('.mat', '') for x in
            os.listdir(join(settings.data_directory, 'models'))
            if '.xml' in x or '.mat' in x]

def check_for_model(name):
    """Check for model, case insensitive, and ignore periods and underscores"""
    def min_name(n):
        return n.lower().replace('.','').replace(' ','').replace('_','')
    for x in get_model_list():
        if min_name(name)==min_name(x):
            return x
    return None

@timing
def load_model(model_id, model_dir, genome_id, model_timestamp, pmid, session):
    """Load a model into the database.

    Arguments
    ---------

    model_id: bigg_id for the model
    
    model_dir: the directory where models are stored.

    genome_id: id for the loaded genome annotation.

    model_timestamp: a timestamp for the model.

    pmid: a publication PMID for the model.

    """
    
    # check for a genome annotation for this model
    genome_db = session.query(base.Genome).filter_by(bioproject_id=genome_id).first()
    if genome_db is None:
        raise Exception('Genbank file %s for model %s not found in the database' %
                        (genome_id, model_id))

    # check that the model doesn't already exist
    if session.query(Model).filter_by(bigg_id=model_id).count() > 0:
        raise Exception('Model %s already loaded' % model_id)

    # apply id normalization
    logging.debug('Parsing SBML')
    model, old_ids = parse.load_and_normalize(model_id, model_dir)
    # Load the model components. Remember: ORDER MATTERS! So don't mess around.
    logging.debug('Loading independent objects')
    independent.loadModel(session, model, genome_db.id, model_timestamp, pmid)
    independent.loadComponents(session, [model])
    independent.loadReactions(session, [model])

    logging.debug('Loading dependent objects')
    dependent.loadModelGenes(session, [model])
    dependent.loadModelCompartmentalizedComponent(session, [model])
    dependent.loadModelReaction(session, [model])
    dependent.loadGeneReactionMatrix(session, [model])
    dependent.loadReactionMatrix(session, [model])
    dependent.loadModelCount(session, [model])
    dependent.loadOldIdtoSynonyms(session, old_ids)

    session.commit()

# -*- coding: utf-8 -*-

from ome import base, settings, components, timing
from ome.models import Model
import cobra.io
from ome.loading import AlreadyLoadedError
from ome.loading.model_loading import loading_methods, parse
from ome.dumping.model_dumping import dump_model
import os
from os.path import join
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
def load_model(model_filepath, genome_id, model_timestamp, pmid, session, 
               dump_directory=settings.model_dump_directory):
    """Load a model into the database. Returns the bigg_id for the new model.

    Arguments
    ---------

    model_filepath: the path to the file where model is stored.

    genome_id: id for the loaded genome annotation.

    model_timestamp: a timestamp for the model.

    pmid: a publication PMID for the model.

    """

    # apply id normalization
    logging.debug('Parsing SBML')
    model, old_parsed_ids = parse.load_and_normalize(model_filepath)
    model_bigg_id = model.id
    
    # check that the model doesn't already exist
    if session.query(Model).filter_by(bigg_id=model_bigg_id).count() > 0:
        raise AlreadyLoadedError('Model %s already loaded' % model_bigg_id)

    # check for a genome annotation for this model
    genome_db = session.query(base.Genome).filter_by(bioproject_id=genome_id).first()
    if genome_db is None:
        raise Exception('Genbank file %s for model %s not found in the database' %
                        (genome_id, model_bigg_id))

    # Load the model objects. Remember: ORDER MATTERS! So don't mess around.
    logging.debug('Loading objects for model {}'.format(model.id))
    model_database_id = loading_methods.load_model(session, model, genome_db.id,
                                                   model_timestamp, pmid)

    # metabolites/components and linkouts
    # get compartment names
    if os.path.exists(settings.compartment_names_file):
        with open(settings.compartment_names_file, 'r') as f:
            compartment_names = {}
            for line in f.readlines():
                sp = [x.strip() for x in line.split('\t')]
                try:
                    compartment_names[sp[0]] = sp[1]
                except IndexError:
                    continue
    else:
        logging.warn('No compartment names file')
        compartment_names = {}
    loading_methods.load_metabolites(session, model_database_id, model,
                                     compartment_names,
                                     old_parsed_ids['metabolites'])

    # reactions
    model_db_rxn_ids = loading_methods.load_reactions(session, model_database_id,
                                                      model, old_parsed_ids['reactions'])

    # genes
    loading_methods.load_genes(session, model_database_id, model,
                               model_db_rxn_ids)

    # count model objects for the model summary web page
    loading_methods.load_model_count(session, model_database_id)

    session.commit()
    
    if dump_directory:
        cobra_model = dump_model(model_bigg_id)
        # make folder if it doesn't exist
        try:
            os.makedirs(dump_directory)
        except OSError:
            pass
        cobra.io.write_sbml_model(cobra_model, join(dump_directory, model_bigg_id + '.xml'))
    
    return model_bigg_id

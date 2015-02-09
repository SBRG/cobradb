from ome import base
from ome.models import *
from ome.loading.model_loading import parse

from tornado.escape import url_escape
import json
import escher
import logging
import sys

def load_maps_from_server(session, drop_maps=False):
    if drop_maps:
        logging.info('Dropping Escher maps')
        connection = base.engine.connect()
        trans = connection.begin()
        try:
            connection.execute('DROP TABLE escher_map, escher_map_matrix CASCADE;')
            trans.commit()
        except:
            logging.warn('Could not drop Escher tables')
            trans.rollback()
        logging.info('Creating tables')
        base.Base.metadata.create_all()

    logging.info('Getting index')
    index = escher.plots.server_index()

    loaded_models = (session
                     .query(Model.bigg_id, Model.id)
                     .all())
    matching_models = [x for x in loaded_models
                       if x[0] in [m['model_name'] for m in index['models']]]
                       
    for model_bigg_id, model_id in matching_models:
        maps = [(m['map_name'], m['organism']) for m in index['maps'] if
                m['map_name'].split('.')[0] == model_bigg_id]
        for map_name, org in maps:
            map_json = escher.plots.map_json_for_name(map_name)
            # map_url = (escher.urls.get_url('map_download', source='web', protocol='https') +
            #            '/'.join([url_escape(x, plus=False) for x in [org, map_name + '.json']]))
            load_the_map(session, model_id, map_name, map_json)
            
def load_the_map(session, model_id, map_name, map_json):
    high_priority = ['central', 'glycolysis']
    priority = 5 if any([s in map_name for s in high_priority]) else 1

    escher_map_db = (session
                     .query(EscherMap)
                     .filter(EscherMap.map_name == map_name)
                     .first())
    if escher_map_db is None:
        logging.info('Creating map %s' % map_name)
        escher_map_db = EscherMap(map_name=map_name, model_id=model_id,
                                  priority=priority, map_data=bytes(map_json))
        session.add(escher_map_db)
        session.commit()
    else:
        logging.info('Map %s already in the database' % map_name)
        if escher_map_db.model_id != model_id:
            model_bigg_id = (session
                             .query(Model.bigg_id)
                             .filter(Model.id == model_id)
                             .first())[0]
            logging.warn('Map %s does not match model %s' % (map_name,
                                                             model_bigg_id))

    map_object = json.loads(map_json)
            
    logging.info('Adding reactions')
    for element_id, reaction in map_object[1]['reactions'].iteritems():
        # check for an existing mat row
        mat_db = (session
                  .query(EscherMapMatrix)
                  .join(ModelReaction, ModelReaction.id == EscherMapMatrix.ome_id)
                  .join(Reaction)
                  .filter(EscherMapMatrix.escher_map_id == escher_map_db.id)
                  .filter(Reaction.bigg_id == reaction['bigg_id'])
                  .first())
        if mat_db is None:
            # find the model reaction
            model_reaction_id = (session
                                 .query(ModelReaction.id)
                                 .join(Reaction)
                                 .filter(Reaction.bigg_id == reaction['bigg_id'])
                                 .filter(ModelReaction.model_id == model_id)
                                 .first())[0]
            mat_db = EscherMapMatrix(escher_map_id=escher_map_db.id,
                                     ome_id=model_reaction_id, 
                                     escher_map_element_id=element_id,
                                     type='model_reaction')
            session.add(mat_db)
            
    logging.info('Adding metabolites')
    for element_id, node in map_object[1]['nodes'].iteritems():
        if node['node_type'] != 'metabolite':
            continue
        metabolite = node
            
        # split the bigg_id
        try:
            met_id, comp_id = parse.split_compartment(metabolite['bigg_id'])
        except Exception:
            logging.warn('Could not split compartment for metabolite %s' % metabolite['bigg_id'])
        # check for an existing mat row
        mat_db = (session
                  .query(EscherMapMatrix)
                  .join(ModelCompartmentalizedComponent,
                        ModelCompartmentalizedComponent.id == EscherMapMatrix.ome_id)
                  .join(CompartmentalizedComponent,
                        CompartmentalizedComponent.id == ModelCompartmentalizedComponent.compartmentalized_component_id)
                  .join(Metabolite,
                        Metabolite.id == CompartmentalizedComponent.component_id)
                  .join(Compartment,
                        Compartment.id == CompartmentalizedComponent.compartment_id)
                  .filter(EscherMapMatrix.escher_map_id == escher_map_db.id)
                  .filter(Metabolite.bigg_id == met_id)
                  .filter(Compartment.bigg_id == comp_id)
                  .first())
        if mat_db is None:
            # find the compartmentalized compartment
            model_comp_comp_id = (session
                                  .query(ModelCompartmentalizedComponent.id)
                                  .join(CompartmentalizedComponent,
                                        CompartmentalizedComponent.id == ModelCompartmentalizedComponent.compartmentalized_component_id)
                                  .join(Metabolite,
                                        Metabolite.id == CompartmentalizedComponent.component_id)
                                  .join(Compartment,
                                        Compartment.id == CompartmentalizedComponent.compartment_id)
                                  .join(Model)
                                  .filter(Compartment.bigg_id == comp_id)
                                  .filter(Metabolite.bigg_id == met_id)
                                  .filter(Model.id == model_id)
                                  .first())[0]
            mat_db = EscherMapMatrix(escher_map_id=escher_map_db.id,
                                     ome_id=model_comp_comp_id,
                                     escher_map_element_id=element_id,
                                     type='model_compartmentalized_component')
            session.add(mat_db)
    session.commit()

if __name__=="__main__":
    logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)

    session = base.Session()
    load_maps_from_server(session, drop_maps=True)
    session.close()

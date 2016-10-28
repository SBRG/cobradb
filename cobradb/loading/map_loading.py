# -*- coding: utf-8 -*-

from cobradb import base
from cobradb.models import *
from cobradb.loading import parse

from tornado.escape import url_escape
import json
import escher
import logging
import sys
import re
import six

def load_maps_from_server(session, drop_maps=False):
    if drop_maps:
        logging.info('Dropping Escher maps')
        connection = base.engine.connect()
        trans = connection.begin()
        try:
            connection.execute('TRUNCATE escher_map, escher_map_matrix CASCADE;')
            trans.commit()
        except:
            logging.warn('Could not drop Escher tables')
            trans.rollback()

    logging.info('Getting index')
    index = escher.plots.server_index()

    loaded_models = (session
                     .query(Model.bigg_id, Model.id)
                     .all())
    matching_models = [x for x in loaded_models
                       if x[0] in [m['model_name'] for m in index['models']]
                          # TODO remove: trick for matching E coli core to e_coli_core
                          or x[0] == 'e_coli_core' and 'E coli core' in [m['model_name'] for m in index['models']]]

    for model_bigg_id, model_id in matching_models:
        maps = [(m['map_name'], m['organism']) for m in index['maps'] if
                m['map_name'].split('.')[0] == model_bigg_id or
                # TODO remove: trick for matching E coli core to e_coli_core
                m['map_name'].split('.')[0] == 'E coli core' and model_bigg_id == 'e_coli_core']
        for map_name, org in maps:
            map_json = escher.plots.map_json_for_name(map_name)
            load_the_map(session, model_id, map_name, map_json)

def load_the_map(session, model_id, map_name, map_json):
    size = sys.getsizeof(map_json)
    if size > 1e6:
        logging.info('Skipping Escher map {} because it is too large ({:.2e} bytes)'
                     .format(map_name, size))
        return 1

    warning_num = 5

    high_priority = ['central', 'glycolysis']
    priority = (5 if any([s in map_name.lower() for s in high_priority]) else 1)

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
    reaction_warnings = 0
    for element_id, reaction in six.iteritems(map_object[1]['reactions']):
        # deal with reaction copies
        map_reaction_bigg_id = re.sub(r'_copy[0-9]+$', '', reaction['bigg_id'])
        # check for an existing mat row
        mat_db = (session
                  .query(EscherMapMatrix)
                  .join(ModelReaction, ModelReaction.id == EscherMapMatrix.ome_id)
                  .join(Reaction)
                  .filter(EscherMapMatrix.escher_map_id == escher_map_db.id)
                  .filter(Reaction.bigg_id == map_reaction_bigg_id)
                  .first())
        if mat_db is None:
            # find the model reaction
            model_reaction_db = (session
                                 .query(ModelReaction.id)
                                 .join(Reaction)
                                 .filter(Reaction.bigg_id == map_reaction_bigg_id)
                                 .filter(ModelReaction.model_id == model_id)
                                 .first())
            if model_reaction_db is None:
                if reaction_warnings <= warning_num:
                    msg = ('Could not find reaction %s in model for map %s' % (map_reaction_bigg_id,
                                                                               map_name))
                    if reaction_warnings == warning_num:
                        msg += ' (Warnings limited to %d)' % warning_num
                    logging.warn(msg)
                    reaction_warnings += 1
                continue
            model_reaction_id = model_reaction_db[0]
            mat_db = EscherMapMatrix(escher_map_id=escher_map_db.id,
                                     ome_id=model_reaction_id,
                                     escher_map_element_id=element_id,
                                     type='model_reaction')
            session.add(mat_db)

    logging.info('Adding metabolites')
    comp_comp_warnings = 0
    for element_id, node in six.iteritems(map_object[1]['nodes']):
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
            model_comp_comp_db = (session
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
                                  .first())
            if model_comp_comp_db is None:
                if comp_comp_warnings <= warning_num:
                    msg = ('Could not find compartmentalized component %s in model for map %s' %
                           ('%s_%s' % (met_id, comp_id), map_name))
                    if comp_comp_warnings == warning_num:
                        msg += ' (Warnings limited to %d)' % warning_num
                    logging.warn(msg)
                    comp_comp_warnings += 1
                continue
            model_comp_comp_id = model_comp_comp_db[0]
            mat_db = EscherMapMatrix(escher_map_id=escher_map_db.id,
                                     ome_id=model_comp_comp_id,
                                     escher_map_element_id=element_id,
                                     type='model_compartmentalized_component')
            session.add(mat_db)
    session.commit()

    return 0

if __name__=="__main__":
    logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)

    session = base.Session()
    load_maps_from_server(session, drop_maps=True)
    session.close()

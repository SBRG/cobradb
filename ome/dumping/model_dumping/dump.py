# -*- coding: utf-8 -*-

from ome import base, settings, components, timing
from ome.base import Session
from ome.models import *

import cobra.core
import logging
from itertools import izip, repeat

try:
    from ipy_progressbar import ProgressBar
except ImportError:
    pass

@timing
def dump_model(bigg_id, session):
    
    # find the model
    model_db = (session
                .query(Model)
                .filter(Model.bigg_id == bigg_id)
                .first())

    if model_db is None:
        session.commit()
        raise Exception('Could not find model %s' % bigg_id)

    model = cobra.core.Model(bigg_id)

    # reactions
    logging.info('Dumping reactions')
    reactions_db = (session
                    .query(Reaction, ModelReaction)
                    .join(ModelReaction, ModelReaction.reaction_id == Reaction.id)
                    .filter(ModelReaction.model_id == model_db.id)
                    .all())
    reactions = []
    for r_db, mr_db in reactions_db:
        r = cobra.core.Reaction(r_db.bigg_id)
        r.name = r_db.name
        r.gene_reaction_rule = mr_db.gene_reaction_rule
        r.lower_bound = mr_db.lower_bound
        r.upper_bound = mr_db.upper_bound
        r.objective_coefficient = mr_db.objective_coefficient
        reactions.append(r)
    model.add_reactions(reactions) 

    # metabolites
    logging.info('Dumping metabolites')
    metabolites_db = (session
                      .query(Component.bigg_id, Compartment.bigg_id)
                      .join(CompartmentalizedComponent,
                            Component.id == CompartmentalizedComponent.component_id)
                      .join(Compartment,
                            Compartment.id == CompartmentalizedComponent.compartment_id)
                      .join(ModelCompartmentalizedComponent,
                            CompartmentalizedComponent.id == ModelCompartmentalizedComponent.compartmentalized_component_id)
                      .filter(ModelCompartmentalizedComponent.model_id == model_db.id)
                      .filter(Component.type == 'metabolite')
                      .all())
    metabolites = []
    for component_id, compartment_id in metabolites_db:
        m = cobra.core.Metabolite(id=component_id + '_' + compartment_id)
        metabolites.append(m)
    model.add_metabolites(metabolites) 

    # reaction matrix
    logging.info('Dumping reaction matrix')
    matrix_db = (session
                 .query(ReactionMatrix.stoichiometry, Reaction.bigg_id,
                        Component.bigg_id, Compartment.bigg_id)
                 # component, compartment
                 .join(CompartmentalizedComponent,
                       ReactionMatrix.compartmentalized_component_id == CompartmentalizedComponent.id)
                 .join(Component,
                       CompartmentalizedComponent.component_id == Component.id)
                 .join(Compartment,
                       CompartmentalizedComponent.compartment_id == Compartment.id)
                 # reaction
                 .join(Reaction,
                       ReactionMatrix.reaction_id == Reaction.id)
                 .join(ModelReaction,
                       Reaction.id == ModelReaction.reaction_id)
                 .filter(ModelReaction.model_id == model_db.id)
                 .filter(Component.type == 'metabolite')
                 .all())

    for stoich, reaction_id, component_id, compartment_id in matrix_db:
        r = model.reactions.get_by_id(reaction_id)
        m = model.metabolites.get_by_id(component_id + '_' + compartment_id)
        r.add_metabolites({ m: stoich }) 

    session.commit()

    return model

def dump_universal_model(session):
    
    model = cobra.core.Model('Universal model')

    # reaction matrix
    logging.info('Dumping reaction matrix')
    matrix_db = (session
                 .query(ReactionMatrix.stoichiometry, Reaction, ModelReaction,
                        Model.bigg_id, Component, Compartment.bigg_id)
                 # component, compartment
                 .join(CompartmentalizedComponent,
                       ReactionMatrix.compartmentalized_component_id == CompartmentalizedComponent.id)
                 .join(Component,
                       CompartmentalizedComponent.component_id == Component.id)
                 .join(Compartment,
                       CompartmentalizedComponent.compartment_id == Compartment.id)
                 # reaction
                 .join(Reaction,
                       ReactionMatrix.reaction_id == Reaction.id)
                 .join(ModelReaction,
                       Reaction.id == ModelReaction.reaction_id)
                 .join(Model,
                       ModelReaction.model_id == Model.id)
                 .filter(Component.type == 'metabolite')
                 # .limit(100)
                 .all())

    def assign_reaction(reaction, db_reaction, db_model_reaction, model_bigg_id):
        reaction.bigg_id = '%s (%s)' % (db_reaction.bigg_id, model_bigg_id)
        reaction.name = db_reaction.name
        reaction.gene_reaction_rule = db_model_reaction.gpr
        reaction.lower_bound = float(db_model_reaction.lower_bound)
        reaction.upper_bound = float(db_model_reaction.upper_bound)
        reaction.objective_coefficient = float(db_model_reaction.objective_coefficient)

    try:
        it = izip(matrix_db, ProgressBar(len(matrix_db)))
    except NameError:
        it = izip(matrix_db, repeat(None))
    
    for (stoich, r_db, mr_db, model_bigg_id, component, compartment_id), _ in it:
        # assign the reaction
        try:
            r = model.reactions.get_by_id(r_db.bigg_id)
            # choose the most lenient reversibility
            if mr_db.lower_bound < 0 and r.lower_bound >= 0:
                assign_reaction(r, r_db, mr_db, model_bigg_id)
        except KeyError:
            r = cobra.core.Reaction(r_db.bigg_id)
            assign_reaction(r, r_db, mr_db, model_bigg_id)
            model.add_reaction(r)

        # add mets
        met_bigg_id = component.bigg_id + '_' + compartment_id
        try:
            m = model.metabolites.get_by_id(met_bigg_id)
        except KeyError:
            m = cobra.core.Metabolite(met_bigg_id)
            m.bigg_id = '%s (%s)' % (component.bigg_id, component.kegg_id)
        r.add_metabolites({ m: float(stoich) }) 

    # TODO genes

    session.commit()

    return model

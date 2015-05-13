# -*- coding: utf-8 -*-

from ome import base, settings, components, timing
from ome.base import Session
from ome.models import *
from ome.util import increment_id

import cobra.core
import logging
from itertools import izip, repeat
from collections import defaultdict
from copy import copy

try:
    from ipy_progressbar import ProgressBar
except ImportError:
    pass

@timing
def dump_model(bigg_id):
    session = Session()
    
    # find the model
    model_db = (session
                .query(Model)
                .filter(Model.bigg_id == bigg_id)
                .first())

    if model_db is None:
        session.commit()
        session.close()
        raise Exception('Could not find model %s' % bigg_id)

    model = cobra.core.Model(str(bigg_id))

    # reactions
    logging.debug('Dumping reactions')
    reactions_db = (session
                    .query(Reaction, ModelReaction, Synonym)
                    .join(ModelReaction, ModelReaction.reaction_id == Reaction.id)
                    .join(OldIDSynonym, OldIDSynonym.ome_id == ModelReaction.id)
                    .join(Synonym, Synonym.id == OldIDSynonym.synonym_id)
                    .filter(ModelReaction.model_id == model_db.id)
                    .all())

    # make dictionaries and cast results
    result_dicts = []
    for r_db, mr_db, synonym_db in reactions_db:
        d = {}
        d['bigg_id'] = str(r_db.bigg_id)
        d['name'] = str(r_db.name)
        d['gene_reaction_rule'] = str(mr_db.gene_reaction_rule)
        d['lower_bound'] = float(mr_db.lower_bound)
        d['upper_bound'] = float(mr_db.upper_bound)
        d['objective_coefficient'] = float(mr_db.objective_coefficient)
        d['original_bigg_id'] = str(synonym_db.synonym)
        result_dicts.append(d)

    def filter_duplicates(result_dicts):
        """Find the reactions with multiple ModelReactions and increment names."""
        tups_by_bigg_id = defaultdict(list)
        # for each ModelReaction
        for d in result_dicts:
            # add to duplicates
            tups_by_bigg_id[d['bigg_id']].append(d)
        # duplicates have multiple ModelReactions
        duplicates = {k: v for k, v in tups_by_bigg_id.iteritems() if len(v) > 1}
        for bigg_id, dup_dicts in duplicates.iteritems():
            # add -copy1, copy2, etc. to the bigg ids for the duplicates
            last = bigg_id
            for d in dup_dicts:
                last = increment_id(last, 'copy')
                d['bigg_id'] = last
        return result_dicts

    # fix duplicates
    result_filtered = filter_duplicates(result_dicts)

    reactions = []
    for result_dict in result_filtered:
        r = cobra.core.Reaction(result_dict['bigg_id'])
        r.name = result_dict['name']
        r.gene_reaction_rule = result_dict['gene_reaction_rule']
        r.lower_bound = result_dict['lower_bound']
        r.upper_bound = result_dict['upper_bound']
        r.objective_coefficient = result_dict['objective_coefficient']
        r.notes = {'original_bigg_id': result_dict['original_bigg_id']}
        reactions.append(r)
    model.add_reactions(reactions) 

    # metabolites
    logging.debug('Dumping metabolites')
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
        if component_id is not None and compartment_id is not None:
            m = cobra.core.Metabolite(id=str(component_id + '_' + compartment_id), compartment=str(compartment_id))
            metabolites.append(m)
    model.add_metabolites(metabolites) 

    # reaction matrix
    logging.debug('Dumping reaction matrix')
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
                 .distinct() # make sure we don't duplicate
                 .all())

    for stoich, reaction_id, component_id, compartment_id in matrix_db:
        try:
            m = model.metabolites.get_by_id(str(component_id + '_' + compartment_id))
        except KeyError:
            logging.warn('Metabolite not found %s in compartment %s for reaction %s' % \
                         (component_id, compartment_id, reaction_id))
            continue
        if reaction_id in model.reactions:
            # check again that we don't duplicate
            r = model.reactions.get_by_id(reaction_id)
            if m not in r.metabolites:
                r.add_metabolites({ m: float(stoich) }) 
        else:
            # try incremented ids
            while True:
                reaction_id = increment_id(reaction_id, 'copy')
                try:
                    # check again that we don't duplicate
                    r = model.reactions.get_by_id(reaction_id)
                    if m not in r.metabolites:
                        r.add_metabolites({ m: float(stoich) }) 
                except KeyError:
                    break

    # dump the gene names
    gene_names = (session
                .query(Gene.bigg_id, Gene.name)
                .join(ModelGene)
                .filter(ModelGene.model_id == model_db.id)
                .all())
    
    for gene_id, gene_name in gene_names:
        try:
            model.genes.get_by_id(gene_id).name = str(gene_name)
        except:
            try:
                model.genes.get_by_id(gene_name).name = str(gene_id)
            except:
                logging.warn('Gene not found %s with name %s' %
                            (gene_id, gene_name))

    session.commit()
    session.close()

    return model


# def dump_universal_model():
#     session = Session()
    
#     model = cobra.core.Model('Universal model')

#     # reaction matrix
#     logging.debug('Dumping reaction matrix')
#     matrix_db = (session
#                  .query(ReactionMatrix.stoichiometry, Reaction, ModelReaction,
#                         Model.bigg_id, Component, Compartment.bigg_id)
#                  # component, compartment
#                  .join(CompartmentalizedComponent,
#                        ReactionMatrix.compartmentalized_component_id == CompartmentalizedComponent.id)
#                  .join(Component,
#                        CompartmentalizedComponent.component_id == Component.id)
#                  .join(Compartment,
#                        CompartmentalizedComponent.compartment_id == Compartment.id)
#                  # reaction
#                  .join(Reaction,
#                        ReactionMatrix.reaction_id == Reaction.id)
#                  .join(ModelReaction,
#                        Reaction.id == ModelReaction.reaction_id)
#                  .join(Model,
#                        ModelReaction.model_id == Model.id)
#                  .filter(Component.type == 'metabolite')
#                  # .limit(100)
#                  .all())

#     def assign_reaction(reaction, db_reaction, db_model_reaction, model_bigg_id):
#         reaction.bigg_id = '%s (%s)' % (db_reaction.bigg_id, model_bigg_id)
#         reaction.name = str(db_reaction.name)
#         reaction.gene_reaction_rule = str(db_model_reaction.gpr)
#         reaction.lower_bound = float(db_model_reaction.lower_bound)
#         reaction.upper_bound = float(db_model_reaction.upper_bound)
#         reaction.objective_coefficient = float(db_model_reaction.objective_coefficient)

#     try:
#         it = izip(matrix_db, ProgressBar(len(matrix_db)))
#     except NameError:
#         it = izip(matrix_db, repeat(None))
    
#     for (stoich, r_db, mr_db, model_bigg_id, component, compartment_id), _ in it:
#         # assign the reaction
#         try:
#             r = model.reactions.get_by_id(r_db.bigg_id)
#             # choose the most lenient reversibility
#             if mr_db.lower_bound < 0 and r.lower_bound >= 0:
#                 assign_reaction(r, r_db, mr_db, model_bigg_id)
#         except KeyError:
#             r = cobra.core.Reaction(r_db.bigg_id)
#             assign_reaction(r, r_db, mr_db, model_bigg_id)
#             model.add_reaction(r)

#         # add mets
#         met_bigg_id = component.bigg_id + '_' + compartment_id
#         try:
#             m = model.metabolites.get_by_id(met_bigg_id)
#         except KeyError:
#             m = cobra.core.Metabolite(met_bigg_id)
#             m.bigg_id = '%s (%s)' % (component.bigg_id, component.kegg_id)
#         r.add_metabolites({ m: float(stoich) }) 

#     session.commit()
#     session.close()

#     return model

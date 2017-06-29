# -*- coding: utf-8 -*-

from ome import base, settings, components, timing
from ome.base import Session
from ome.models import *
from ome.util import increment_id, make_reaction_copy_id

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

    model = cobra.core.Model(bigg_id)

    # genes
    logging.debug('Dumping genes')
    # get genes and original bigg ids (might be multiple)
    genes_db = (session
                .query(Gene.bigg_id, Gene.name, Synonym.synonym)
                .join(ModelGene, ModelGene.gene_id == Gene.id)
                .join(OldIDSynonym, OldIDSynonym.ome_id == ModelGene.id)
                .filter(OldIDSynonym.type == 'model_gene')
                .join(Synonym, Synonym.id == OldIDSynonym.synonym_id)
                .filter(ModelGene.model_id == model_db.id))
    gene_names = []
    old_gene_ids_dict = defaultdict(list)
    for gene_id, gene_name, old_id in genes_db:
        if gene_id not in old_gene_ids_dict:
            gene_names.append((gene_id, gene_name))
        old_gene_ids_dict[gene_id].append(old_id)

    for gene_id, gene_name in gene_names:
        gene = cobra.core.Gene(gene_id)
        gene.name = gene_name
        gene.notes = {'original_bigg_ids': old_gene_ids_dict[gene_id]}
        model.genes.append(gene)

    # reactions
    logging.debug('Dumping reactions')
    # get original bigg ids (might be multiple)
    reactions_db = (session
                    .query(ModelReaction, Reaction, Synonym.synonym)
                    .join(Reaction)
                    .join(OldIDSynonym, OldIDSynonym.ome_id == ModelReaction.id)
                    .filter(OldIDSynonym.type == 'model_reaction')
                    .filter(Synonym.type == 'reaction')
                    .join(Synonym, Synonym.id == OldIDSynonym.synonym_id)
                    .filter(ModelReaction.model_id == model_db.id))
    reactions_model_reactions = []
    found_model_reactions = set()
    old_reaction_ids_dict = defaultdict(list)
    for model_reaction, reaction, old_id in reactions_db:
        # there may be multiple model reactions for a given bigg_id
        if model_reaction.id not in found_model_reactions:
            reactions_model_reactions.append((model_reaction, reaction))
            found_model_reactions.add(model_reaction.id)
        old_reaction_ids_dict[reaction.bigg_id].append(old_id)

    # make dictionaries and cast results
    result_dicts = []
    for mr_db, r_db in reactions_model_reactions:
        d = {}
        d['bigg_id'] = r_db.bigg_id
        d['name'] = r_db.name
        d['gene_reaction_rule'] = mr_db.gene_reaction_rule
        d['lower_bound'] = float(mr_db.lower_bound)
        d['upper_bound'] = float(mr_db.upper_bound)
        d['objective_coefficient'] = float(mr_db.objective_coefficient)
        d['original_bigg_ids'] = old_reaction_ids_dict[r_db.bigg_id]
        d['subsystem'] = mr_db.subsystem
        d['copy_number'] = int(mr_db.copy_number)
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
            # add _copy1, copy2, etc. to the bigg ids for the duplicates
            for d in dup_dicts:
                d['bigg_id'] = make_reaction_copy_id(bigg_id, d['copy_number'])

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
        r.notes = {'original_bigg_ids': result_dict['original_bigg_ids']}
        r.subsystem = result_dict['subsystem']
        reactions.append(r)
    model.add_reactions(reactions)

    # metabolites
    logging.debug('Dumping metabolites')
    # get original bigg ids (might be multiple)
    metabolites_db = (session
                      .query(Metabolite.bigg_id,
                             Metabolite.name,
                             ModelCompartmentalizedComponent.formula,
                             ModelCompartmentalizedComponent.charge,
                             Compartment.bigg_id,
                             Synonym.synonym)
                      .join(CompartmentalizedComponent)
                      .join(Compartment)
                      .join(ModelCompartmentalizedComponent)
                      .join(OldIDSynonym, OldIDSynonym.ome_id == ModelCompartmentalizedComponent.id)
                      .filter(OldIDSynonym.type == 'model_compartmentalized_component')
                      .filter(Synonym.type == 'compartmentalized_component')
                      .join(Synonym)
                      .filter(ModelCompartmentalizedComponent.model_id == model_db.id))
    metabolite_names = []
    old_metabolite_ids_dict = defaultdict(list)
    for metabolite_id, metabolite_name, formula, charge, compartment_id, old_id in metabolites_db:
        if metabolite_id + '_' + compartment_id not in old_metabolite_ids_dict:
            metabolite_names.append((metabolite_id, metabolite_name, formula, charge, compartment_id))
        old_metabolite_ids_dict[metabolite_id + '_' + compartment_id].append(old_id)

    metabolites = []
    compartments = set()
    for component_id, component_name, formula, charge, compartment_id in metabolite_names:
        if component_id is not None and compartment_id is not None:
            m = cobra.core.Metabolite(id=component_id + '_' + compartment_id,
                                      compartment=compartment_id,
                                      formula=formula)
            m.charge = charge
            m.name = component_name
            m.notes = {'original_bigg_ids': old_metabolite_ids_dict[component_id + '_' + compartment_id]}
            compartments.add(compartment_id)
            metabolites.append(m)
    model.add_metabolites(metabolites)

    # compartments
    compartment_db = (session.query(Compartment)
                      .filter(Compartment.bigg_id.in_(compartments)))
    model.compartments = {i.bigg_id: i.name for i in compartment_db}

    # reaction matrix
    logging.debug('Dumping reaction matrix')
    matrix_db = (session
                 .query(ReactionMatrix.stoichiometry, Reaction.bigg_id,
                        Component.bigg_id, Compartment.bigg_id)
                 # component, compartment
                 .join(CompartmentalizedComponent)
                 .join(Component)
                 .join(Compartment)
                 # reaction
                 .join(Reaction)
                 .join(ModelReaction)
                 .filter(ModelReaction.model_id == model_db.id)
                 .distinct())  # make sure we don't duplicate

    # load metabolites
    for stoich, reaction_id, component_id, compartment_id in matrix_db:
        try:
            m = model.metabolites.get_by_id(component_id + '_' + compartment_id)
        except KeyError:
            logging.warn('Metabolite not found %s in compartment %s for reaction %s' % \
                         (component_id, compartment_id, reaction_id))
            continue
        # add to reactions
        if reaction_id in model.reactions:
            # check again that we don't duplicate
            r = model.reactions.get_by_id(reaction_id)
            if m not in r.metabolites:
                r.add_metabolites({m: float(stoich)})
        else:
            # try incremented ids
            while True:
                reaction_id = increment_id(reaction_id, 'copy')
                try:
                    # check again that we don't duplicate
                    r = model.reactions.get_by_id(reaction_id)
                    if m not in r.metabolites:
                        r.add_metabolites({m: float(stoich)})
                except KeyError:
                    break

    session.commit()
    session.close()

    cobra.manipulation.annotate.add_SBO(model)

    return model

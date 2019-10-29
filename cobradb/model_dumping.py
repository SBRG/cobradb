# -*- coding: utf-8 -*-

from cobradb.models import *
from cobradb.util import increment_id, make_reaction_copy_id, timing

from sqlalchemy import and_
import cobra.core
import logging
from itertools import repeat
from collections import defaultdict
from copy import copy
import six

try:
    from ipy_progressbar import ProgressBar
except ImportError:
    pass


def _none_to_str(val):
    return '' if val is None else val


def _make_annotation_lookup(db_links):
    """Make a lookup dictionary from a list of flat external DB links"""
    lookup = defaultdict(lambda: defaultdict(set))
    for res in db_links:
        # skip old_bigg_id because it will be in notes
        if res[0] not in ['old_bigg_id', 'deprecated']:
            lookup[res[2]][res[0]].add(res[1])
    # return lists instead of sets
    return {
        bigg_id: {source: list(vals) for source, vals in links.items()}
        for bigg_id, links in lookup.items()
    }


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

    # get gene annotations
    gene_db_links = _make_annotation_lookup(
        session
        .query(DataSource.bigg_id, Synonym.synonym, Gene.bigg_id)
        .join(Synonym)
        .join(Gene, Gene.id == Synonym.ome_id)
        .filter(Synonym.type == 'gene')
    )

    for gene_id, gene_name in gene_names:
        gene = cobra.core.Gene(gene_id)
        gene.name = _none_to_str(gene_name)
        gene.notes = {'original_bigg_ids': old_gene_ids_dict[gene_id]}
        gene.annotation = gene_db_links.get(gene_id, {})
        # add SBO terms
        gene.annotation['sbo'] = 'SBO:0000243'
        model.genes.append(gene)

    # reactions
    logging.debug('Dumping reactions')
    # get original bigg ids (might be multiple)
    reactions_db = (session
                    .query(ModelReaction, Reaction, Synonym.synonym)
                    .join(Reaction)
                    .outerjoin(OldIDSynonym,
                               and_(OldIDSynonym.ome_id == ModelReaction.id,
                                    OldIDSynonym.type == 'model_reaction'))
                    .outerjoin(Synonym,
                               and_(Synonym.id == OldIDSynonym.synonym_id,
                                    Synonym.type == 'reaction'))
                    .filter(ModelReaction.model_id == model_db.id))
    reactions_model_reactions = []
    found_model_reactions = set()
    old_reaction_ids_dict = defaultdict(list)
    for model_reaction, reaction, old_id in reactions_db:
        # there may be multiple model reactions for a given bigg_id
        if model_reaction.id not in found_model_reactions:
            reactions_model_reactions.append((model_reaction, reaction))
            found_model_reactions.add(model_reaction.id)
        if old_id is not None:
            old_reaction_ids_dict[reaction.bigg_id].append(old_id)

    # get reaction annotations
    reaction_db_links = _make_annotation_lookup(
        session
        .query(DataSource.bigg_id, Synonym.synonym, Reaction.bigg_id)
        .join(Synonym)
        .join(Reaction, Reaction.id == Synonym.ome_id)
        .filter(Synonym.type == 'reaction')
    )

    # make dictionaries and cast results
    result_dicts = []
    for mr_db, r_db in reactions_model_reactions:
        d = {}
        d['bigg_id'] = r_db.bigg_id
        d['name'] = r_db.name
        d['gene_reaction_rule'] = mr_db.gene_reaction_rule
        d['lower_bound'] = mr_db.lower_bound
        d['upper_bound'] = mr_db.upper_bound
        d['objective_coefficient'] = mr_db.objective_coefficient
        d['original_bigg_ids'] = old_reaction_ids_dict[r_db.bigg_id]
        d['subsystem'] = mr_db.subsystem
        d['annotation'] = reaction_db_links.get(r_db.bigg_id, {})
        # add SBO terms
        if r_db.bigg_id.startswith('BIOMASS_'):
            d['annotation']['sbo'] = 'SBO:0000629'
        elif r_db.bigg_id.startswith('EX_'):
            d['annotation']['sbo'] = 'SBO:0000627'
        elif r_db.bigg_id.startswith('DM_'):
            d['annotation']['sbo'] = 'SBO:0000628'
        elif r_db.bigg_id.startswith('SK_'):
            d['annotation']['sbo'] = 'SBO:0000632'
        else:
            # assume non-transport. will update for transporters later
            d['annotation']['sbo'] = 'SBO:0000176'
        # specify bigg id
        d['annotation']['bigg.reaction'] = [r_db.bigg_id]
        d['copy_number'] = mr_db.copy_number
        result_dicts.append(d)

    def filter_duplicates(result_dicts):
        """Find the reactions with multiple ModelReactions and increment names."""
        tups_by_bigg_id = defaultdict(list)
        # for each ModelReaction
        for d in result_dicts:
            # add to duplicates
            tups_by_bigg_id[d['bigg_id']].append(d)
        # duplicates have multiple ModelReactions
        duplicates = {k: v for k, v in six.iteritems(tups_by_bigg_id) if len(v) > 1}
        for bigg_id, dup_dicts in six.iteritems(duplicates):
            # add _copy1, copy2, etc. to the bigg ids for the duplicates
            for d in dup_dicts:
                d['bigg_id'] = make_reaction_copy_id(bigg_id, d['copy_number'])

        return result_dicts

    # fix duplicates
    result_filtered = filter_duplicates(result_dicts)

    reactions = []
    objectives = {}
    for result_dict in result_filtered:
        r = cobra.core.Reaction(result_dict['bigg_id'])
        r.name = _none_to_str(result_dict['name'])
        r.gene_reaction_rule = result_dict['gene_reaction_rule']
        r.lower_bound = result_dict['lower_bound']
        r.upper_bound = result_dict['upper_bound']
        r.notes = {'original_bigg_ids': result_dict['original_bigg_ids']}
        r.subsystem = result_dict['subsystem']
        r.annotation = result_dict['annotation']
        reactions.append(r)

        objectives[r.id] = result_dict['objective_coefficient']
    model.add_reactions(reactions)

    for k, v in six.iteritems(objectives):
        model.reactions.get_by_id(k).objective_coefficient = v

    # metabolites
    logging.debug('Dumping metabolites')
    # get original bigg ids (might be multiple)
    metabolites_db = (session
                      .query(Component.bigg_id,
                             Component.name,
                             ModelCompartmentalizedComponent.formula,
                             ModelCompartmentalizedComponent.charge,
                             Compartment.bigg_id,
                             Synonym.synonym)
                      .join(CompartmentalizedComponent,
                            CompartmentalizedComponent.component_id == Component.id)  # noqa
                      .join(Compartment,
                            Compartment.id == CompartmentalizedComponent.compartment_id)  # noqa
                      .join(ModelCompartmentalizedComponent,
                            ModelCompartmentalizedComponent.compartmentalized_component_id == CompartmentalizedComponent.id)  # noqa
                      .join(OldIDSynonym, OldIDSynonym.ome_id == ModelCompartmentalizedComponent.id)  # noqa
                      .filter(OldIDSynonym.type == 'model_compartmentalized_component')  # noqa
                      .filter(Synonym.type == 'compartmentalized_component')
                      .join(Synonym)
                      .filter(ModelCompartmentalizedComponent.model_id == model_db.id))  # noqa
    metabolite_names = []
    old_metabolite_ids_dict = defaultdict(list)
    for metabolite_id, metabolite_name, formula, charge, compartment_id, old_id in metabolites_db:
        if metabolite_id + '_' + compartment_id not in old_metabolite_ids_dict:
            metabolite_names.append((metabolite_id, metabolite_name, formula, charge, compartment_id))
        old_metabolite_ids_dict[metabolite_id + '_' + compartment_id].append(old_id)

    # get metabolite annotations
    metabolite_db_links = _make_annotation_lookup(
        session
        .query(DataSource.bigg_id, Synonym.synonym, Component.bigg_id)
        .join(Synonym)
        .join(Component, Component.id == Synonym.ome_id)
        .filter(Synonym.type == 'component')
    )

    metabolites = []
    compartments = set()
    for component_id, component_name, formula, charge, compartment_id in metabolite_names:
        if component_id is not None and compartment_id is not None:
            m = cobra.core.Metabolite(id=component_id + '_' + compartment_id,
                                      compartment=compartment_id,
                                      formula=formula)
            m.charge = charge
            m.name = _none_to_str(component_name)
            m.notes = {'original_bigg_ids': old_metabolite_ids_dict[component_id + '_' + compartment_id]}
            m.annotation = metabolite_db_links.get(component_id, {})
            # specify bigg id
            m.annotation['bigg.metabolite'] = [component_id]
            m.annotation['sbo'] = 'SBO:0000247'
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
                 .join(CompartmentalizedComponent,
                       CompartmentalizedComponent.id == ReactionMatrix.compartmentalized_component_id)  # noqa
                 .join(Component,
                       Component.id == CompartmentalizedComponent.component_id)
                 .join(Compartment,
                       Compartment.id == CompartmentalizedComponent.compartment_id)  # noqa
                 # reaction
                 .join(Reaction, Reaction.id == ReactionMatrix.reaction_id)
                 .join(ModelReaction, ModelReaction.reaction_id == Reaction.id)
                 .filter(ModelReaction.model_id == model_db.id)
                 .distinct())  # make sure we don't duplicate

    # load metabolites
    compartments_for_reaction = defaultdict(set)
    for stoich, reaction_id, component_id, compartment_id in matrix_db:
        try:
            m = model.metabolites.get_by_id(component_id + '_' + compartment_id)
        except KeyError:
            logging.warning('Metabolite not found %s in compartment %s for reaction %s' %
                            (component_id, compartment_id, reaction_id))
            continue
        # add to reactions
        if reaction_id in model.reactions:
            # check again that we don't duplicate
            r = model.reactions.get_by_id(reaction_id)
            if m not in r.metabolites:
                r.add_metabolites({m: float(stoich)})
                # check for transporters and update sbo term
                compartments_for_reaction[reaction_id].add(m.compartment)
                if len(compartments_for_reaction[reaction_id]) > 1 and r.annotation['sbo'] == 'SBO:0000176':  # noqa
                    r.annotation['sbo'] = 'SBO:0000185'
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

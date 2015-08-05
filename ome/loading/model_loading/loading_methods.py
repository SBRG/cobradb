# -*- coding: utf-8 -*-

from ome import base
from ome.base import NotFoundError
from ome.models import *
from ome.components import *
from ome.loading.model_loading import parse
from ome.util import increment_id, check_pseudoreaction, load_tsv

from sqlalchemy.orm.exc import MultipleResultsFound, NoResultFound
from sqlalchemy import func
import re
import logging
from collections import defaultdict
import os
from itertools import ifilter


def _fix_name(name):
    """Clean up descriptive names"""
    return name.strip('_ ')


def _get_data_source(session, name):
    """Get or create a DataSource."""
    data_source_db = (session
                      .query(base.DataSource)
                      .filter(base.DataSource.name == name)
                      .first())
    if data_source_db is None:
        data_source_db = base.DataSource(name=name)
        session.add(data_source_db)
        session.commit()
    return data_source_db.id


def load_new_model(session, model, genome_db_id, pub_ref, published_filename):
    """Load the model.

    Arguments:
    ---------

    session: A SQLAlchemy session.

    model: A COBRApy model.

    genome_db_id: The database ID of the genome. Can be None.

    pub_ref: a publication PMID or doi for the model, as a string like this:

        doi:10.1128/ecosalplus.10.2.1

        pmid:21988831

        Can be None

    Returns:
    -------

    The database ID of the new model row.

    """
    model_db = Model(bigg_id=model.id, genome_id=genome_db_id,
                     description=model.description,
                     published_filename=published_filename)
    session.add(model_db)
    if pub_ref is not None:
        try:
            ref_type, ref_id = pub_ref.split(':')
        except ValueError:
            logging.warn('Bad publication reference {}'.format(pub_ref))
        else:
            publication_db = (session
                              .query(base.Publication)
                              .filter(base.Publication.reference_type==ref_type)
                              .filter(base.Publication.reference_id==ref_id)
                              .first())
            if publication_db is None:
                publication_db = base.Publication(reference_type=ref_type,
                                                  reference_id=ref_id)
                session.add(publication_db)
                session.commit()
            publication_model_db = (session
                                    .query(base.PublicationModel)
                                    .filter(base.PublicationModel.publication_id == publication_db.id)
                                    .filter(base.PublicationModel.model_id == model_db.id)
                                    .first())
            if publication_model_db is None:
                publication_model_db = base.PublicationModel(model_id=model_db.id,
                                                             publication_id=publication_db.id)
                session.add(publication_model_db)
    session.commit()
    return model_db.id


def _load_metabolite_linkouts(session, cobra_metabolite, metabolite_database_id):
    """Load new linkouts even ones that are pointing to previously created universal
    metabolites.

    The only scenario where we don't load a linkout is if the external id and
    metabolite is exactly the same as a previous linkout.

    """

    # parse the notes
    def parse_linkout_str(id):
        if id is None:
            return None
        id_string = str(id)
        for s in ['{', '}', '[', ']', '&apos;', "'",]:
            id_string = id_string.replace(s, '')
        return id_string.strip()

    linkouts = ['KEGGID',
                'CASNUMBER',
                'SEED',
                'METACYC',
                'CHEBI',
                'BRENDA',
                'UPA',
                'HMDB',
                'BIOPATH',
                'REACTOME',
                'LIPIDMAPS',
                'CASID',
                'PUBCHEM ID']

    for external_source, v in cobra_metabolite.notes.iteritems():
        # ignore formulas
        if external_source.lower() in ['formula', 'formula1', 'none']:
            continue
        # check if linkout matches the list
        external_source = external_source.upper()
        v = v[0]
        if not external_source in linkouts:
            logging.warning('The linkout type {} is not in our list'
                            .format(external_source))
            continue
        if '&apos' in v:
            ids = [parse_linkout_str(x) for x in v.split(',')]
        else:
            ids = [parse_linkout_str(v)]

        for external_id in ids:
            if external_id.lower() in ['na', 'none']:
                continue
            exists = (session
                      .query(LinkOut)
                      .filter(LinkOut.external_id == external_id)
                      .filter(LinkOut.external_source == external_source)
                      .filter(LinkOut.type == 'metabolite')
                      .filter(LinkOut.ome_id == metabolite_database_id)
                      .count() > 0)
            if not exists:
                linkout = LinkOut(external_id=external_id,
                                  external_source=external_source,
                                  type='metabolite',
                                  ome_id=metabolite_database_id)
                session.add(linkout)


def load_metabolites(session, model_id, model, compartment_names,
                     old_metabolite_ids):
    """Load the metabolites as components and model components.

    Arguments:
    ---------

    session: An SQLAlchemy session.

    model_id: The database ID for the model.

    model: The COBRApy model.

    old_metabolite_ids: A dictionary where keys are new IDs and values are old
    IDs for compartmentalized metabolites.

    """

    # only grab this once
    data_source_id = _get_data_source(session, 'old_id')

    # for each metabolite in the model
    for metabolite in model.metabolites:
        try:
            component_bigg_id, compartment_bigg_id = parse.split_compartment(metabolite.id)
        except Exception:
            logging.error(('Could not find compartment for metabolite %s in'
                            'model %s' % (metabolite.id, model.id)))
            continue

        # If there is no metabolite, add a new one.
        # TODO we could also double check these ID matches with linkouts and formula
        metabolite_db = (session
                         .query(Metabolite)
                         .filter(Metabolite.bigg_id == component_bigg_id)
                         .first())

        # Look for the formula in these places
        formula_fns = [lambda m: getattr(m, 'formula', None), # support cobra v0.3 and 0.4
                       lambda m: m.notes.get('FORMULA', None),
                       lambda m: m.notes.get('FORMULA1', None)]
        # Cast to string, but not for None
        strip_str_or_none = lambda v: str(v).strip() if v is not None else None
        # Ignore the empty string
        ignore_empty_str = lambda s: s if s != '' else None
        # Use a generator for lazy evaluation
        values = (ignore_empty_str(strip_str_or_none(formula_fn(metabolite)))
                  for formula_fn in formula_fns)
        # Get the first non-null result. Otherwise _formula = None.
        _formula = next(ifilter(None, values), None)

        # if necessary, add the new metabolite, and keep track of the ID
        if metabolite_db is None:
            # check for missing info
            if metabolite.name.strip() == '':
                metabolite.name = ''
                logging.warn('No name for metabolite {} in model {}. TODO Solution: Add name from other models.'
                             .format(metabolite.id, model.id))

            # make the new metabolite
            metabolite_db = Metabolite(bigg_id=component_bigg_id,
                                       name=_fix_name(metabolite.name),
                                       formula=_formula)
            session.add(metabolite_db)
            session.commit()

        elif metabolite_db.formula is None:
            metabolite_db.formula = _formula

        # load the linkouts for the universal metabolite
        _load_metabolite_linkouts(session, metabolite, metabolite_db.id)

        # if there is no compartment, add a new one
        compartment_db = (session
                          .query(Compartment)
                          .filter(Compartment.bigg_id == compartment_bigg_id)
                          .first())
        if compartment_db is None:
            try:
                name = compartment_names[compartment_bigg_id]
            except KeyError:
                logging.warn('No name found for compartment %s' % compartment_bigg_id)
                name = ''
            compartment_db = Compartment(bigg_id=compartment_bigg_id, name=name)
            session.add(compartment_db)
            session.commit()

        # if there is no compartmentalized compartment, add a new one
        comp_component_db = (session
                             .query(CompartmentalizedComponent)
                             .filter(CompartmentalizedComponent.component_id == metabolite_db.id)
                             .filter(CompartmentalizedComponent.compartment_id == compartment_db.id)
                             .first())
        if comp_component_db is None:
            comp_component_db = CompartmentalizedComponent(component_id=metabolite_db.id,
                                                           compartment_id=compartment_db.id)
            session.add(comp_component_db)
            session.commit()

        # if there is no model compartmentalized compartment, add a new one
        model_comp_comp_db = (session
                              .query(ModelCompartmentalizedComponent)
                              .filter(ModelCompartmentalizedComponent.compartmentalized_component_id == comp_component_db.id)
                              .filter(ModelCompartmentalizedComponent.model_id == model_id)
                              .first())
        if model_comp_comp_db is None:
            model_comp_comp_db = ModelCompartmentalizedComponent(model_id=model_id,
                                                                 compartmentalized_component_id=comp_component_db.id)
            session.add(model_comp_comp_db)
            session.commit()

        # add synonyms
        old_bigg_id_c = old_metabolite_ids[metabolite.id]
        synonym_db = (session
                      .query(base.Synonym)
                      .filter(base.Synonym.ome_id == metabolite_db.id)
                      .filter(base.Synonym.synonym == old_bigg_id_c)
                      .filter(base.Synonym.type == 'metabolite')
                      .first())
        if synonym_db is None:
            synonym_db = base.Synonym(type='metabolite',
                                      ome_id=metabolite_db.id,
                                      synonym=old_bigg_id_c,
                                      synonym_data_source_id=data_source_id)
            session.add(synonym_db)
            session.commit()

        # add OldIDSynonym
        old_id_db = (session
                     .query(OldIDSynonym)
                     .filter(OldIDSynonym.ome_id == model_comp_comp_db.id)
                     .filter(OldIDSynonym.synonym_id == synonym_db.id)
                     .first())
        if old_id_db is None:
            old_id_db = base.OldIDSynonym(type='model_compartmentalized_component',
                                          ome_id=model_comp_comp_db.id,
                                          synonym_id=synonym_db.id)
            session.add(old_id_db)
            session.commit()


def _new_reaction(session, reaction, bigg_id, reaction_hash, model_db_id, model,
                  is_pseudoreaction):
    """Add a new universal reaction with reaction matrix rows."""

    reaction_db = Reaction(bigg_id=bigg_id, name=_fix_name(reaction.name),
                           reaction_hash=reaction_hash,
                           pseudoreaction=is_pseudoreaction)
    session.add(reaction_db)
    session.commit

    # for each reactant, add to the reaction matrix
    for metabolite, stoich in reaction.metabolites.iteritems():
        try:
            component_bigg_id, compartment_bigg_id = parse.split_compartment(metabolite.id)
        except NotFoundError:
            logging.error('Could not split metabolite %s in model %s' % (metabolite.id, model.id))
            continue

        # get the component in the model
        comp_comp_db = (session
                        .query(CompartmentalizedComponent)
                        .join(Component,
                              Component.id == CompartmentalizedComponent.component_id)
                        .join(Compartment,
                              Compartment.id == CompartmentalizedComponent.compartment_id)
                        .join(ModelCompartmentalizedComponent,
                              ModelCompartmentalizedComponent.compartmentalized_component_id == CompartmentalizedComponent.id)
                        .filter(Component.bigg_id == component_bigg_id)
                        .filter(Compartment.bigg_id == compartment_bigg_id)
                        .filter(ModelCompartmentalizedComponent.model_id == model_db_id)
                        .first())
        if comp_comp_db is None:
            logging.error('Could not find metabolite {!s} for model {!s} in the database'
                        .format(metabolite.id, model.id))
            continue

        # check if the reaction matrix row already exists
        found_reaction_matrix = (session
                                .query(ReactionMatrix)
                                .filter(ReactionMatrix.reaction_id == reaction_db.id)
                                .filter(ReactionMatrix.compartmentalized_component_id == comp_comp_db.id)
                                .count() > 0)
        if not found_reaction_matrix:
            new_object = ReactionMatrix(reaction_id=reaction_db.id,
                                        compartmentalized_component_id=comp_comp_db.id,
                                        stoichiometry=stoich)
            session.add(new_object)
        else:
            logging.debug('ReactionMatrix row already present for model {!s} metabolite {!s} reaction {!s}'
                        .format(model.id, metabolite.id, reaction.id))

    return reaction_db

def load_reactions(session, model_db_id, model, old_reaction_ids):
    """Load the reactions and stoichiometries into the model.

    TODO if the reaction is already loaded, we need to check the stoichometry
    has. If that doesn't match, then add a new reaction with an incremented ID
    (e.g. ACALD_1)

    Arguments
    ---------

    session: An SQLAlchemy session.

    model_db_id: The database ID for the model.

    model: The COBRApy model.

    old_reaction_ids: A dictionary where keys are new IDs and values are old IDs
    for reactions.

    Returns
    -------

    A dictionary with keys for reactions in the model and values for the
    associated bigg_id in the database.

    """

    # only grab this once
    data_source_id = _get_data_source(session, 'old_id')

    # get reaction id_prefs
    id_prefs = load_tsv(settings.reaction_id_prefs)
    def _check_id_prefs(an_id, versus_id):
        """Return True if an_id is preferred over versus_id."""
        for row in id_prefs:
            try:
                idx1 = row.index(an_id)
                idx2 = row.index(versus_id)
            except ValueError:
                continue

            return idx1 < idx2
        return False

    # get reaction hash_prefs
    hash_prefs = load_tsv(settings.reaction_hash_prefs)
    def _check_hash_prefs(a_hash):
        """Return the preferred BiGG ID for a_hash, or None."""
        for row in hash_prefs:
            if row[0] == a_hash:
                return row[1]
        return None

    model_db_rxn_ids = {}
    for reaction in model.reactions:
        # get the reaction
        reaction_db = (session
                       .query(Reaction)
                       .filter(Reaction.bigg_id == reaction.id)
                       .first())

        # check for pseudoreaction
        is_pseudoreaction = check_pseudoreaction(reaction.id)

        # calculate the hash
        reaction_hash = parse.hash_reaction(reaction)
        hash_db = (session
                   .query(Reaction)
                   .filter(Reaction.reaction_hash == reaction_hash)
                   .filter(Reaction.pseudoreaction == is_pseudoreaction)
                   .first())

        # bigg_id match  hash match b==h  pseudoreaction  example                   function
        #  n               n               n            first GAPD                _new_reaction (1)
        #  n               n               y            first EX_glc_e            _new_reaction (1)
        #  y               n               n            incorrect GAPD            _new_reaction & increment (2)
        #  y               n               y            incorrect EX_glc_e        _new_reaction & increment (2)
        #  n               y               n            GAPDH after GAPD          reaction = hash_reaction (3a)
        #  n               y               y            EX_glc__e after EX_glc_e  reaction = hash_reaction (3a)
        #  y               y         n     n            ?                         reaction = hash_reaction (3a)
        #  y               y         n     y            ?                         reaction = hash_reaction (3a)
        #  y               y         y     n            second GAPD               reaction = bigg_reaction (3b)
        #  y               y         y     y            second EX_glc_e           reaction = bigg_reaction (3b)
        # NOTE: only check pseudoreaction hash against other pseudoreactions

        def _find_new_incremented_id(session, original_id):
            """Look for a reaction bigg_id that is not already taken."""
            new_id = increment_id(original_id)
            while True:
                if session.query(Reaction).filter(Reaction.bigg_id == new_id).first() is None:
                    return new_id
                new_id = increment_id(new_id)

        preferred_id = _check_hash_prefs(reaction_hash)
        # (0) If there is a preferred ID, make that the new ID, and increment any old IDs
        if preferred_id is not None:
            # if the reaction already matches, just continue
            if hash_db is not None and hash_db.bigg_id == preferred_id:
                reaction_db = hash_db
            # otherwise, make the new reaction
            else:
                # if existing reactions match the preferred reaction find a new,
                # incremented id for the existing match
                preferred_id_db = session.query(Reaction).filter(Reaction.bigg_id == preferred_id).first()
                if preferred_id_db is not None:
                    new_id = _find_new_incremented_id(session, preferred_id)
                    logging.warn('Incrementing database reaction {} to {} and prefering {} (from model {}) based on hash preferences'
                                .format(preferred_id, new_id, preferred_id, model.id))
                    preferred_id_db.bigg_id = new_id
                    session.commit()

                # make a new reaction for the preferred_id
                reaction_db = _new_reaction(session, reaction, preferred_id,
                                            reaction_hash, model_db_id, model,
                                            is_pseudoreaction)

        # (1) no bigg_id matches, no stoichiometry match or pseudoreaction, then
        # make a new reaction
        elif reaction_db is None and hash_db is None:
            reaction_db = _new_reaction(session, reaction, reaction.id,
                                        reaction_hash, model_db_id, model,
                                        is_pseudoreaction)

        # (2) bigg_id matches, but not the hash, then increment the BIGG_ID
        elif reaction_db is not None and hash_db is None:
            # loop until we find a non-matching find non-matching ID
            new_id = _find_new_incremented_id(session, reaction.id)
            logging.warn('Incrementing bigg_id {} to {} (from model {}) based on conflicting reaction hash'
                        .format(reaction.id, new_id, model.id))
            reaction_db = _new_reaction(session, reaction, new_id,
                                        reaction_hash, model_db_id, model,
                                        is_pseudoreaction)

        # (3) but found a stoichiometry match, then use the hash reaction match.
        elif hash_db is not None:
            # WARNING TODO this requires that loaded metabolites always match on
            # bigg_id, which should be the case.

            # (3a)
            if reaction_db is None or reaction_db.id != hash_db.id:
                is_preferred = _check_id_prefs(reaction.id, hash_db.bigg_id)
                if is_preferred:
                    logging.warn('Switching database reaction {} to bigg_id {} based on reaction hash and id_prefs file'
                                .format(hash_db.bigg_id, reaction.id, model.id))
                    hash_db.bigg_id = reaction.id
                    session.commit()
                reaction_db = hash_db
            # (3b) BIGG ID matches a reaction with the same hash, then just continue
            else:
                pass

        else:
            raise Exception('Should not get here')

        # remember the changed ids
        model_db_rxn_ids[reaction.id] = reaction_db.bigg_id

        # get the model reaction
        model_reaction_db = (session
                             .query(ModelReaction)
                             .filter(ModelReaction.reaction_id == reaction_db.id)
                             .filter(ModelReaction.model_id == model_db_id)
                             .filter(ModelReaction.lower_bound == reaction.lower_bound)
                             .filter(ModelReaction.upper_bound == reaction.upper_bound)
                             .filter(ModelReaction.gene_reaction_rule == reaction.gene_reaction_rule)
                             .filter(ModelReaction.objective_coefficient == reaction.objective_coefficient)
                             .first())
        if model_reaction_db is None:
            # make a new reaction
            model_reaction_db = ModelReaction(model_id=model_db_id,
                                              reaction_id=reaction_db.id,
                                              gene_reaction_rule=reaction.gene_reaction_rule,
                                              original_gene_reaction_rule=reaction.gene_reaction_rule,
                                              upper_bound=reaction.upper_bound,
                                              lower_bound=reaction.lower_bound,
                                              objective_coefficient=reaction.objective_coefficient)
            session.add(model_reaction_db)
            session.commit()

        # add synonyms
        #
        # get the id from the published model
        old_bigg_id = old_reaction_ids[reaction.id]
        # add a synonym
        synonym_db = (session
                      .query(base.Synonym)
                      .filter(base.Synonym.ome_id == reaction_db.id)
                      .filter(base.Synonym.synonym == old_bigg_id)
                      .first())
        if synonym_db is None:
            synonym_db = base.Synonym(type='reaction',
                                      ome_id=reaction_db.id,
                                      synonym=old_bigg_id,
                                      synonym_data_source_id=data_source_id)
            session.add(synonym_db)
            session.commit()

        # add OldIDSynonym
        old_id_db = (session
                     .query(OldIDSynonym)
                     .filter(OldIDSynonym.ome_id == model_reaction_db.id)
                     .filter(OldIDSynonym.synonym_id == synonym_db.id)
                     .first())
        if old_id_db is None:
            old_id_db = base.OldIDSynonym(type='model_reaction',
                                          ome_id=model_reaction_db.id,
                                          synonym_id=synonym_db.id)
            session.add(old_id_db)
            session.commit()

    return model_db_rxn_ids


# find gene functions
def _match_gene_by_fns(fn_list, session, gene_id, chromosome_ids):
    """Go through each funciton and look for a match.

    """
    for fn in fn_list:
        match, is_alternative_transcript = fn(session, gene_id, chromosome_ids)
        if len(match) > 0:
            if len(match) > 1:
                logging.warn('Multiple matches for gene {} with function {}. Using the first match.'
                            .format(gene_id, fn.__name__))
            return match[0], is_alternative_transcript
    return None, False


def _by_bigg_id(session, gene_id, chromosome_ids):
    # look for a matching model gene
    gene_db = (session
               .query(Gene)
               .filter(func.lower(Gene.bigg_id) == func.lower(gene_id))
               .filter(Gene.chromosome_id.in_(chromosome_ids))
               .all())
    return gene_db, False


def _by_name(session, gene_id, chromosome_ids):
    gene_db = (session
               .query(Gene)
               .filter(func.lower(Gene.name) == func.lower(gene_id))
               .filter(Gene.chromosome_id.in_(chromosome_ids))
               .all())
    return gene_db, False


def _by_synonym(session, gene_id, chromosome_ids):
    gene_db = (session
               .query(Gene)
               .join(Synonym, Synonym.ome_id == Gene.id)
               .filter(Gene.chromosome_id.in_(chromosome_ids))
               .filter(func.lower(Synonym.synonym) == func.lower(gene_id))
               .all())
    return gene_db, False


def _by_alternative_transcript(session, gene_id, chromosome_ids):
    """Function to check for the alternative transcript match."""
    check = re.match(r'(.*)_AT[0-9]{1,2}$', gene_id)
    if not check:
        gene_db = []
    else:
        # find the old gene
        gene_db = (session
                   .query(Gene)
                   .filter(Gene.chromosome_id.in_(chromosome_ids))
                   .filter(func.lower(Gene.bigg_id) == func.lower(check.group(1)))
                   .filter(Gene.alternative_transcript_of.is_(None))
                   .all())
    return gene_db, True


def _by_alternative_transcript_name(session, gene_id, chromosome_ids):
    """Function to check for the alternative transcript match."""
    check = re.match(r'(.*)_AT[0-9]{1,2}$', gene_id)
    if not check:
        gene_db = []
    else:
        # find the old gene
        gene_db = (session
                   .query(Gene)
                   .filter(Gene.chromosome_id.in_(chromosome_ids))
                   .filter(func.lower(Gene.name) == func.lower(check.group(1)))
                   .filter(Gene.alternative_transcript_of.is_(None))
                   .all())
    return gene_db, True


def _by_alternative_transcript_synonym(session, gene_id, chromosome_ids):
    """Function to check for the alternative transcript match."""
    check = re.match(r'(.*)_AT[0-9]{1,2}$', gene_id)
    if not check:
        gene_db = []
    else:
        # find the old gene
        gene_db = (session
                   .query(Gene)
                   .join(Synonym, Synonym.ome_id == Gene.id)
                   .filter(Gene.chromosome_id.in_(chromosome_ids))
                   .filter(func.lower(Synonym.synonym) == func.lower(check.group(1)))
                   .filter(Gene.alternative_transcript_of.is_(None))
                   .all())
    return gene_db, True


def _by_bigg_id_no_underscore(session, gene_id, chromosome_ids):
    """Matches for T maritima genes"""
    # look for a matching model gene
    gene_db = (session
               .query(Gene)
               .filter(func.lower(Gene.bigg_id) == func.lower(gene_id.replace('_', '')))
               .filter(Gene.chromosome_id.in_(chromosome_ids))
               .all())
    return gene_db, False


def _replace_gene_str(rule, old_gene, new_gene):
    return re.sub(r'\b'+old_gene+r'\b', new_gene, rule)


def load_genes(session, model_db_id, model, model_db_rxn_ids):
    """Load the genes for this model.

    Arguments:
    ---------

    session: An SQLAlchemy session.

    model_db_id: The database ID for the model.

    model: The COBRApy model.

    model_db_rxn_ids: A dictionary with keys for reactions in the model and
    values for the associated bigg_id in the database.

    """
    # find the model in the db
    model_db = session.query(Model).get(model_db_id)

    # find the chromosomes in the db
    chromosome_ids = (session
                       .query(Chromosome.id)
                       .filter(Chromosome.genome_id == model_db.genome_id)
                       .all())
    if len(chromosome_ids) == 0:
        logging.warn('No chromosomes for model %s' % model_db.bigg_id)

    # keep track of the gene-reaction associations
    gene_bigg_id_to_model_reaction_db_ids = defaultdict(set)
    for reaction in model.reactions:
        # find the ModelReaction that corresponds to this particular reaction in
        # the model
        model_reaction_db = (session
                             .query(ModelReaction)
                             .join(Reaction)
                             .filter(ModelReaction.model_id == model_db_id)
                             .filter(Reaction.bigg_id == model_db_rxn_ids[reaction.id])
                             .all())
        if model_reaction_db is None:
            logging.error('Could not find ModelReaction of {} (originally {}) in model {}. Cannot load GeneReactionMatrix entries'
                          .format(model_db_rxn_ids[reaction.id], reaction.id, model.id))
            continue
        for gene in reaction.genes:
            for mr in model_reaction_db:
                gene_bigg_id_to_model_reaction_db_ids[gene.id].add(mr.id)

    # load the genes
    for gene in model.genes:
        if len(chromosome_ids) == 0:
            gene_db = None; is_alternative_transcript = False
        else:
            # find a matching gene
            fns = [_by_bigg_id, _by_name, _by_synonym, _by_alternative_transcript,
                   _by_alternative_transcript_name, _by_alternative_transcript_synonym,
                   _by_bigg_id_no_underscore]
            gene_db, is_alternative_transcript = _match_gene_by_fns(fns, session,
                                                                    gene.id,
                                                                    chromosome_ids)

        if not gene_db:
            # add
            if len(chromosome_ids) > 0:
                logging.warn('Gene not in genbank file: {} from model {}'
                            .format(gene.id, model.id))
            ome_gene = {}
            ome_gene['bigg_id'] = gene.id
            ome_gene['name'] = (gene.name if gene.name.strip() != gene.id.strip() else None)
            ome_gene['leftpos'] = None
            ome_gene['rightpos'] = None
            ome_gene['chromosome_id'] = None
            ome_gene['strand'] = None
            ome_gene['info'] = str(gene.annotation)
            ome_gene['mapped_to_genbank'] = False
            gene_db = Gene(**ome_gene)
            session.add(gene_db)
            session.commit()

        elif is_alternative_transcript:
            # duplicate gene for the alternative transcript
            old_gene_db = gene_db
            ome_gene = {}
            ome_gene['bigg_id'] = gene.id
            ome_gene['name'] = old_gene_db.name
            ome_gene['leftpos'] = old_gene_db.leftpos
            ome_gene['rightpos'] = old_gene_db.rightpos
            ome_gene['chromosome_id'] = old_gene_db.chromosome_id
            ome_gene['strand'] = old_gene_db.strand
            ome_gene['info'] = old_gene_db.info
            ome_gene['mapped_to_genbank'] = True
            ome_gene['alternative_transcript_of'] = old_gene_db.id
            gene_db = Gene(**ome_gene)
            session.add(gene_db)
            session.commit()

            # duplicate all the synonyms
            synonyms_db = (session
                           .query(Synonym)
                           .filter(Synonym.ome_id == old_gene_db.id)
                           .all())
            for syn_db in synonyms_db:
                # add a new synonym
                ome_synonym = {}
                ome_synonym['type'] = syn_db.type
                ome_synonym['ome_id'] = gene_db.id
                ome_synonym['synonym'] = syn_db.synonym
                ome_synonym['synonym_data_source_id'] = syn_db.synonym_data_source_id
                synonym_object = Synonym(**ome_synonym)
                session.add(synonym_object)

        # add model gene
        model_gene_db = (session
                         .query(ModelGene)
                         .filter(ModelGene.gene_id == gene_db.id)
                         .filter(ModelGene.model_id == model_db_id)
                         .first())
        if model_gene_db is None:
            model_gene_db = ModelGene(gene_id=gene_db.id,
                                      model_id=model_db_id)
            session.add(model_gene_db)
            session.commit()

        # find model reaction
        try:
            model_reaction_db_ids = gene_bigg_id_to_model_reaction_db_ids[gene.id]
        except KeyError:
            # error message above
            continue

        for mr_db_id in model_reaction_db_ids:
            # add to the GeneReactionMatrix, if not already present
            found_gene_reaction_row = (session
                                       .query(GeneReactionMatrix)
                                       .filter(GeneReactionMatrix.model_gene_id == model_gene_db.id)
                                       .filter(GeneReactionMatrix.model_reaction_id == mr_db_id)
                                       .count() > 0)
            if not found_gene_reaction_row:
                new_object = GeneReactionMatrix(model_gene_id=model_gene_db.id,
                                                model_reaction_id=mr_db_id)
                session.add(new_object)

            # update the gene_reaction_rule if the gene id has changed
            if gene.id != gene_db.bigg_id:
                mr = session.query(ModelReaction).get(mr_db_id)
                new_rule = _replace_gene_str(mr.gene_reaction_rule, gene.id,
                                             gene_db.bigg_id)
                (session
                .query(ModelReaction)
                .filter(ModelReaction.id == mr_db_id)
                .update({ModelReaction.gene_reaction_rule: new_rule}))


def load_model_count(session, model_db_id):
    metabolite_count = (session
                        .query(ModelCompartmentalizedComponent.id)
                        .filter(ModelCompartmentalizedComponent.model_id == model_db_id)
                        .count())
    reaction_count = (session
                      .query(ModelReaction.id)
                      .filter(ModelReaction.model_id == model_db_id)
                      .count())
    gene_count = (session
                  .query(ModelGene.id)
                  .filter(ModelGene.model_id == model_db_id)
                  .count())
    mc = ModelCount(model_id=model_db_id,
                    gene_count=gene_count,
                    metabolite_count=metabolite_count,
                    reaction_count=reaction_count)
    session.add(mc)

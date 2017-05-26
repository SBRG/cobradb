# -*- coding: utf-8 -*-

from ome import base, settings, components, timing
from ome.loading import AlreadyLoadedError
from ome.dumping.model_dumping import dump_model
from ome.base import *
from ome.models import *
from ome.components import *
from ome.loading import parse
from ome.util import (increment_id, check_pseudoreaction, load_tsv,
                      get_or_create_data_source, format_formula, scrub_name,
                      check_none, get_or_create)

from sqlalchemy.orm.exc import MultipleResultsFound, NoResultFound
from sqlalchemy import func
import re
import logging
from collections import defaultdict
import os
from os.path import join, basename, abspath, dirname
from itertools import ifilter
import cobra.io


class GenbankNotFound(Exception):
    pass


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
def load_model(model_filepath, pub_ref, genome_ref, session):
    """Load a model into the database. Returns the bigg_id for the new model.

    Arguments
    ---------

    model_filepath: the path to the file where model is stored.

    pub_ref: a publication PMID or doi for the model, as a tuple like this:

        ('doi', '10.1128/ecosalplus.10.2.1')

        ('pmid', '21988831')

        Can be None

    genome_ref: A tuple specifying the genome accession type and value. The
    first element can be ncbi_accession, ncbi_assembly, or organism.

    session: An instance of base.Session.

    """
    # apply id normalization
    logging.debug('Parsing SBML')
    model, old_parsed_ids = parse.load_and_normalize(model_filepath)
    model_bigg_id = model.id

    # check that the model doesn't already exist
    if session.query(Model).filter_by(bigg_id=model_bigg_id).count() > 0:
        raise AlreadyLoadedError('Model %s already loaded' % model_bigg_id)

    # check for a genome annotation for this model
    if genome_ref is not None and genome_ref[0] == 'organism':
        genome_id = None
        organism = genome_ref[1]
    elif genome_ref is not None and genome_ref[0] in ['ncbi_accession', 'ncbi_assembly']:
        genome_db = (session
                     .query(Genome)
                     .filter(Genome.accession_type==genome_ref[0])
                     .filter(Genome.accession_value==genome_ref[1])
                     .first())
        if genome_db is None:
            raise GenbankNotFound('Genome for model {} not found with genome_ref {}'
                                  .format(model_bigg_id, genome_ref))
        genome_id = genome_db.id
        organism = genome_db.organism
    else:
        logging.info('No Genome reference or organism provided for model {}'.format(model_bigg_id))
        genome_id = None
        organism = None

    # Load the model objects. Remember: ORDER MATTERS! So don't mess around.
    logging.debug('Loading objects for model {}'.format(model.id))
    published_filename = os.path.basename(model_filepath)
    model_database_id = load_new_model(session, model, genome_id, pub_ref,
                                       published_filename, organism)

    # metabolites/components and linkouts
    # get compartment names
    if os.path.exists(settings.compartment_names):
        with open(settings.compartment_names, 'r') as f:
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
    comp_comp_db_ids = load_metabolites(session, model_database_id, model,
                                        compartment_names,
                                        old_parsed_ids['metabolites'])

    # reactions
    model_db_rxn_ids = load_reactions(session, model_database_id, model,
                                      old_parsed_ids['reactions'],
                                      comp_comp_db_ids)

    # genes
    load_genes(session, model_database_id, model, model_db_rxn_ids,
               old_parsed_ids['genes'])

    # count model objects for the model summary web page
    load_model_count(session, model_database_id)

    session.commit()

    return model_bigg_id


def load_new_model(session, model, genome_db_id, pub_ref, published_filename,
                   organism):
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

    organism: The organism. Can be None.

    Returns:
    -------

    The database ID of the new model row.

    """
    model_db = Model(bigg_id=model.id, genome_id=genome_db_id,
                     published_filename=published_filename, organism=organism)
    session.add(model_db)
    if pub_ref is not None:
        ref_type, ref_id = pub_ref
        publication_db = (session
                            .query(Publication)
                            .filter(Publication.reference_type==ref_type)
                            .filter(Publication.reference_id==ref_id)
                            .first())
        if publication_db is None:
            publication_db = Publication(reference_type=ref_type,
                                                reference_id=ref_id)
            session.add(publication_db)
            session.commit()
        publication_model_db = (session
                                .query(PublicationModel)
                                .filter(PublicationModel.publication_id == publication_db.id)
                                .filter(PublicationModel.model_id == model_db.id)
                                .first())
        if publication_model_db is None:
            publication_model_db = PublicationModel(model_id=model_db.id,
                                                            publication_id=publication_db.id)
            session.add(publication_model_db)
    session.commit()
    return model_db.id


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

    Returns
    -------

    comp_comp_db_ids: A dictionary where keys are the original compartmentalized
    metabolite ids and the values are the database IDs for the compartmentalized
    components.

    """
    comp_comp_db_ids = {}

    # only grab this once
    data_source_id = get_or_create_data_source(session, 'old_bigg_id')

    # get metabolite id duplicates
    met_dups = load_tsv(settings.metabolite_duplicates)
    def _check_metabolite_duplicates(bigg_id):
        """Return a new ID if there is a preferred ID, otherwise None."""
        for row in met_dups:
            if bigg_id in row[1:]:
                return row[0]
        return None

    # for each metabolite in the model
    for metabolite in model.metabolites:
        try:
            component_bigg_id, compartment_bigg_id = parse.split_compartment(metabolite.id)
        except Exception:
            logging.error(('Could not find compartment for metabolite %s in'
                            'model %s' % (metabolite.id, model.id)))
            continue

        original_id = None
        preferred = _check_metabolite_duplicates(component_bigg_id)
        if preferred:
            component_bigg_id = preferred
            original_id = component_bigg_id

        # look for the formula in these places
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
        _formula = format_formula(next(ifilter(None, values), None))

        # get charge
        try:
            charge = int(metabolite.charge)
        except Exception:
            if hasattr(metabolite, 'charge') and metabolite.charge is not None:
                logging.debug('Could not convert charge to integer for metabolite {} in model {}: {}'
                              .format(metabolite.id, model.id, metabolite.charge))
            charge = None

        # If there is no metabolite, add a new one.
        metabolite_db = (session
                         .query(Metabolite)
                         .filter(Metabolite.bigg_id == component_bigg_id)
                         .first())

        # if necessary, add the new metabolite, and keep track of the ID
        if metabolite_db is None:
            # make the new metabolite
            metabolite_db = Metabolite(bigg_id=component_bigg_id,
                                       name=scrub_name(getattr(metabolite, 'name', None)))
            session.add(metabolite_db)
            session.commit()

        # add the deprecated id if necessary
        if original_id:
            get_or_create(session, DeprecatedID, deprecated_id=original_id,
                          type='component', ome_id=metabolite_db.id)

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

        # remember for adding the reaction
        comp_comp_db_ids[metabolite.id] = comp_component_db.id

        # if there is no model compartmentalized compartment, add a new one
        model_comp_comp_db = (session
                              .query(ModelCompartmentalizedComponent)
                              .filter(ModelCompartmentalizedComponent.compartmentalized_component_id == comp_component_db.id)
                              .filter(ModelCompartmentalizedComponent.model_id == model_id)
                              .first())
        if model_comp_comp_db is None:
            model_comp_comp_db = ModelCompartmentalizedComponent(model_id=model_id,
                                                                 compartmentalized_component_id=comp_component_db.id,
                                                                 formula=_formula,
                                                                 charge=charge)
            session.add(model_comp_comp_db)
            session.commit()
        else:
            if model_comp_comp_db.formula is None:
                model_comp_comp_db.formula = _formula
            if model_comp_comp_db.charge is None:
                model_comp_comp_db.charge = charge
            session.commit()

        # add synonyms
        for old_bigg_id_c in old_metabolite_ids[metabolite.id]:
            synonym_db = (session
                          .query(Synonym)
                          .filter(Synonym.type == 'compartmentalized_component')
                          .filter(Synonym.ome_id == comp_component_db.id)
                          .filter(Synonym.synonym == old_bigg_id_c)
                          .filter(Synonym.data_source_id == data_source_id)
                          .first())
            if synonym_db is None:
                synonym_db = Synonym(type='compartmentalized_component',
                                     ome_id=comp_component_db.id,
                                     synonym=old_bigg_id_c,
                                     data_source_id=data_source_id)
                session.add(synonym_db)
                session.commit()

            # add OldIDSynonym
            old_id_db = (session
                         .query(OldIDSynonym)
                         .filter(OldIDSynonym.type == 'model_compartmentalized_component')
                         .filter(OldIDSynonym.ome_id == model_comp_comp_db.id)
                         .filter(OldIDSynonym.synonym_id == synonym_db.id)
                         .first())
            if old_id_db is None:
                old_id_db = OldIDSynonym(type='model_compartmentalized_component',
                                         ome_id=model_comp_comp_db.id,
                                         synonym_id=synonym_db.id)
                session.add(old_id_db)
                session.commit()

    return comp_comp_db_ids

def _new_reaction(session, reaction, bigg_id, reaction_hash, model_db_id, model,
                  is_pseudoreaction, comp_comp_db_ids):
    """Add a new universal reaction with reaction matrix rows."""

    # name is optional in cobra 0.4b2. This will probably change back.
    name = check_none(getattr(reaction, 'name', None))
    reaction_db = Reaction(bigg_id=bigg_id, name=scrub_name(name),
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
        try:
            comp_comp_db_id = comp_comp_db_ids[metabolite.id]
        except KeyError:
            logging.error('Could not find metabolite {!s} for model {!s} in the database'
                          .format(metabolite.id, model.id))
            continue

        # check if the reaction matrix row already exists
        found_reaction_matrix = (session
                                .query(ReactionMatrix)
                                .filter(ReactionMatrix.reaction_id == reaction_db.id)
                                .filter(ReactionMatrix.compartmentalized_component_id == comp_comp_db_id)
                                .count() > 0)
        if not found_reaction_matrix:
            new_object = ReactionMatrix(reaction_id=reaction_db.id,
                                        compartmentalized_component_id=comp_comp_db_id,
                                        stoichiometry=stoich)
            session.add(new_object)
        else:
            logging.debug('ReactionMatrix row already present for model {!s} metabolite {!s} reaction {!s}'
                        .format(model.id, metabolite.id, reaction.id))

    return reaction_db

def load_reactions(session, model_db_id, model, old_reaction_ids, comp_comp_db_ids):
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

    comp_comp_db_ids: A dictionary where keys are the original compartmentalized
    metabolite ids and the values are the database IDs for the compartmentalized
    components.

    Returns
    -------

    A dictionary with keys for reaction BiGG IDs in the model and values for the
    associated ModelReaction.id in the database.

    """

    # only grab this once
    data_source_id = get_or_create_data_source(session, 'old_bigg_id')

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
        # if there wasn't a match for the forward hash, also check the reverse hash
        if hash_db:
            reverse_hash_db = None
        else:
            reverse_hash = parse.hash_reaction(parse.reverse_reaction(reaction))
            reverse_hash_db = (session
                               .query(Reaction)
                               .filter(Reaction.reaction_hash == reverse_hash)
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
        # 4a and 4b are 3a and 3b with a reversed reaction

        def _find_new_incremented_id(session, original_id):
            """Look for a reaction bigg_id that is not already taken."""
            new_id = increment_id(original_id)
            while True:
                if session.query(Reaction).filter(Reaction.bigg_id == new_id).first() is None:
                    return new_id
                new_id = increment_id(new_id)

        # Check for a preferred ID in the preferences, based on the forward
        # hash. Don't check the reverse hash in preferences.
        preferred_id = _check_hash_prefs(reaction_hash)

        # no reversed by default
        is_reversed = False

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
                                            is_pseudoreaction, comp_comp_db_ids)

        # (1) no bigg_id matches, no stoichiometry match or pseudoreaction, then
        # make a new reaction
        elif reaction_db is None and hash_db is None and reverse_hash_db is None:
            reaction_db = _new_reaction(session, reaction, reaction.id,
                                        reaction_hash, model_db_id, model,
                                        is_pseudoreaction, comp_comp_db_ids)

        # (2) bigg_id matches, but not the hash, then increment the BIGG_ID
        elif reaction_db is not None and hash_db is None and reverse_hash_db is None:
            # loop until we find a non-matching find non-matching ID
            new_id = _find_new_incremented_id(session, reaction.id)
            logging.warn('Incrementing bigg_id {} to {} (from model {}) based on conflicting reaction hash'
                        .format(reaction.id, new_id, model.id))
            reaction_db = _new_reaction(session, reaction, new_id,
                                        reaction_hash, model_db_id, model,
                                        is_pseudoreaction, comp_comp_db_ids)

        # (3) but found a stoichiometry match, then use the hash reaction match.
        elif hash_db is not None:
            # WARNING TODO this requires that loaded metabolites always match on
            # bigg_id, which should be the case.

            # (3a)
            if reaction_db is None or reaction_db.id != hash_db.id:
                is_preferred = _check_id_prefs(reaction.id, hash_db.bigg_id)
                if is_preferred:
                    logging.warn('Switching database reaction {} to bigg_id {} based on reaction hash and id_prefs file'
                                .format(hash_db.bigg_id, reaction.id))
                    hash_db.bigg_id = reaction.id
                    session.commit()
                reaction_db = hash_db
            # (3b) BIGG ID matches a reaction with the same hash, then just continue
            else:
                pass

        # (4) but found a stoichiometry match, then use the hash reaction match.
        elif reverse_hash_db is not None:
            # WARNING TODO this requires that loaded metabolites always match on
            # bigg_id, which should be the case.

            # Remember to switch upper and lower bounds
            is_reversed = True
            logging.info('Matched {} to {} based on reverse hash'
                         .format(reaction.id, reverse_hash_db.bigg_id))

            # (4a)
            if reaction_db is None or reaction_db.id != reverse_hash_db.id:
                is_preferred = _check_id_prefs(reaction.id, reverse_hash_db.bigg_id)
                if is_preferred:
                    logging.warn('Switching database reaction {} to bigg_id {} based on reversed reaction hash and id_prefs file'
                                .format(reverse_hash_db.bigg_id, reaction.id))
                    reverse_hash_db.bigg_id = reaction.id
                    session.commit()
                reaction_db = reverse_hash_db
            # (4b) BIGG ID matches a reaction with the same hash, then just continue
            else:
                pass

        else:
            raise Exception('Should not get here')

        # If the reaction is reversed, then switch upper and lower bound
        lower_bound = -reaction.upper_bound if is_reversed else reaction.lower_bound
        upper_bound = -reaction.lower_bound if is_reversed else reaction.upper_bound

        # subsystem
        subsystem = check_none(reaction.subsystem.strip())

        # get the model reaction
        model_reaction_db = (session
                             .query(ModelReaction)
                             .filter(ModelReaction.reaction_id == reaction_db.id)
                             .filter(ModelReaction.model_id == model_db_id)
                             .filter(ModelReaction.lower_bound == lower_bound)
                             .filter(ModelReaction.upper_bound == upper_bound)
                             .filter(ModelReaction.gene_reaction_rule == reaction.gene_reaction_rule)
                             .filter(ModelReaction.objective_coefficient == reaction.objective_coefficient)
                             .filter(ModelReaction.subsystem == subsystem)
                             .first())
        if model_reaction_db is None:
            # get the number of existing copies of this reaction in the model
            copy_number = (session
                           .query(ModelReaction)
                           .filter(ModelReaction.reaction_id == reaction_db.id)
                           .filter(ModelReaction.model_id == model_db_id)
                           .count()) + 1
            # make a new reaction
            model_reaction_db = ModelReaction(model_id=model_db_id,
                                              reaction_id=reaction_db.id,
                                              gene_reaction_rule=reaction.gene_reaction_rule,
                                              original_gene_reaction_rule=reaction.gene_reaction_rule,
                                              upper_bound=upper_bound,
                                              lower_bound=lower_bound,
                                              objective_coefficient=reaction.objective_coefficient,
                                              copy_number=copy_number,
                                              subsystem=subsystem)
            session.add(model_reaction_db)
            session.commit()

        # remember the changed ids
        model_db_rxn_ids[reaction.id] = model_reaction_db.id

        # add synonyms
        #
        # get the id from the published model
        for old_bigg_id in old_reaction_ids[reaction.id]:
            # add a synonym
            synonym_db = (session
                          .query(Synonym)
                          .filter(Synonym.type == 'reaction')
                          .filter(Synonym.ome_id == reaction_db.id)
                          .filter(Synonym.synonym == old_bigg_id)
                          .filter(Synonym.data_source_id == data_source_id)
                          .first())
            if synonym_db is None:
                synonym_db = Synonym(type='reaction',
                                     ome_id=reaction_db.id,
                                     synonym=old_bigg_id,
                                     data_source_id=data_source_id)
                session.add(synonym_db)
                session.commit()

            # add OldIDSynonym
            old_id_db = (session
                         .query(OldIDSynonym)
                         .filter(OldIDSynonym.type == 'model_reaction')
                         .filter(OldIDSynonym.ome_id == model_reaction_db.id)
                         .filter(OldIDSynonym.synonym_id == synonym_db.id)
                         .first())
            if old_id_db is None:
                old_id_db = OldIDSynonym(type='model_reaction',
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


def load_genes(session, model_db_id, model, model_db_rxn_ids, old_gene_ids):
    """Load the genes for this model.

    Arguments:
    ---------

    session: An SQLAlchemy session.

    model_db_id: The database ID for the model.

    model: The COBRApy model.

    model_db_rxn_ids: A dictionary with keys for reactions in the model and
    values for the associated ModelReaction.id in the database.

    old_gene_ids: A dictionary where keys are new IDs and values are old IDs for
    genes.

    """
    # only grab this once
    data_source_id = get_or_create_data_source(session, 'old_bigg_id')

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
                             .get(model_db_rxn_ids[reaction.id]))
        if model_reaction_db is None:
            logging.error('Could not find ModelReaction {} for {} in model {}. Cannot load GeneReactionMatrix entries'
                          .format(model_db_rxn_ids[reaction.id], reaction.id, model.id))
            continue
        for gene in reaction.genes:
            gene_bigg_id_to_model_reaction_db_ids[gene.id].add(model_reaction_db.id)

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
            gene_db = Gene(bigg_id=gene.id,
                           # name is optional in cobra 0.4b2. This will probably change back.
                           name=scrub_name(getattr(gene, 'name', None)),
                           mapped_to_genbank=False)
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
                ome_synonym['data_source_id'] = syn_db.data_source_id
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

        # add old gene synonym
        for old_bigg_id in old_gene_ids[gene.id]:
            synonym_db = (session
                          .query(Synonym)
                          .filter(Synonym.type == 'gene')
                          .filter(Synonym.ome_id == gene_db.id)
                          .filter(Synonym.synonym == old_bigg_id)
                          .filter(Synonym.data_source_id == data_source_id)
                          .first())
            if synonym_db is None:
                synonym_db = Synonym(type='gene',
                                     ome_id=gene_db.id,
                                     synonym=old_bigg_id,
                                     data_source_id=data_source_id)
                session.add(synonym_db)
                session.commit()
            # add OldIDSynonym
            old_id_db = (session
                         .query(OldIDSynonym)
                         .filter(OldIDSynonym.type == 'model_gene')
                         .filter(OldIDSynonym.ome_id == model_gene_db.id)
                         .filter(OldIDSynonym.synonym_id == synonym_db.id)
                         .first())
            if old_id_db is None:
                old_id_db = OldIDSynonym(type='model_gene',
                                        ome_id=model_gene_db.id,
                                        synonym_id=synonym_db.id)
                session.add(old_id_db)
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

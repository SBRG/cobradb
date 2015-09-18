# -*- coding: utf-8 -*-

from ome.base import NotFoundError
from ome.util import scrub_gene_id, load_tsv
from ome import settings

import re
import cobra
import cobra.io
from cobra.core import Formula
from os.path import join
import hashlib
import logging
from collections import defaultdict


def hash_reaction(reaction, string_only=False):
    """Generate a unique hash for the metabolites and coefficients of the
    reaction.

    """
    def sorted_mets(reaction):
        return sorted([(m.id, v) for m, v in reaction.metabolites.iteritems()],
                      key=lambda x: x[0])

    if string_only:
        hash_fn = lambda s: s
    else:
        hash_fn = lambda s: hashlib.md5(s).hexdigest()

    return hash_fn(''.join(['%s%.3f' % t for t in sorted_mets(reaction)]))


def load_and_normalize(model_filepath):
    """Load a model, and give it a particular id style"""

    # load the model
    if model_filepath.endswith('.xml'):
        model = cobra.io.read_sbml_model(model_filepath)
    elif model_filepath.endswith('.mat'):
        model = cobra.io.load_matlab_model(model_filepath)
    else:
       raise Exception('The %s file is not a valid filetype', model_filepath)
    # convert the ids
    model, old_ids = convert_ids(model)

    # extract metabolite formulas from names (e.g. for iAF1260)
    model = get_formulas_from_names(model)

    return model, old_ids


def _get_rule_prefs():
    """Get gene_reaction_rule prefs."""
    return load_tsv(settings.gene_reaction_rule_prefs, required_column_num=2)


def _check_rule_prefs(rule_prefs, rule):
    """Check the gene_reaction_rule against the prefs file, and return an existing
    rule or the fixed one."""
    for row in rule_prefs:
        old_rule, new_rule = row
        if old_rule == rule:
            return new_rule
    return rule


def convert_ids(model):
    """Converts metabolite and reaction ids to the new style.

    Returns a tuple with the new model and a dictionary of old ids set up like this:

    {'reactions': {'new_id': 'old_id'},
     'metabolites': {'new_id': 'old_id'},
     'genes': {'new_id': 'old_id'}}

    """
    # loop through the ids:
    metabolite_id_dict = defaultdict(list)
    reaction_id_dict = defaultdict(list)
    gene_id_dict = defaultdict(list)

    # fix metabolites
    for metabolite in model.metabolites:
        new_id = id_for_new_id_style(fix_legacy_id(metabolite.id, use_hyphens=False),
                                     is_metabolite=True)
        metabolite_id_dict[new_id].append(metabolite.id)
        metabolite.id = new_id
    model.metabolites._generate_index()

    # remove boundary metabolites (end in _b and only present in exchanges). Be
    # sure to loop through a static list of ids so the list does not get shorter
    # as the metabolites are deleted
    for metabolite_id in [str(x) for x in model.metabolites]:
        metabolite = model.metabolites.get_by_id(metabolite_id)
        if not metabolite.id.endswith("_b"):
            continue
        for reaction in list(metabolite._reaction):
            if reaction.id.startswith("EX_"):
                metabolite.remove_from_model()
                break
    model.metabolites._generate_index()

    # load fixes for gene_reaction_rule's
    rule_prefs = _get_rule_prefs()

    # separate ids and compartments, and convert to the new_id_style
    for reaction in model.reactions:
        new_id = id_for_new_id_style(fix_legacy_id(reaction.id, use_hyphens=False))
        reaction_id_dict[new_id].append(reaction.id)
        reaction.id = new_id
        # fix the gene reaction rules
        reaction.gene_reaction_rule = _check_rule_prefs(rule_prefs, reaction.gene_reaction_rule)
    model.reactions._generate_index()

    # update the genes
    for gene in list(model.genes):
        new_id = scrub_gene_id(gene.id)
        gene_id_dict[new_id].append(gene.id)
        for reaction in gene.reactions:
            reaction.gene_reaction_rule = re.sub(r'\b'+gene.id+r'\b', new_id,
                                                 reaction.gene_reaction_rule)

    # remove old genes
    for gene in list(model.genes):
        if len(gene.reactions) == 0:
            gene.remove_from_model()

    # fix the model id
    bigg_id = re.sub(r'[^a-zA-Z0-9_]', '_', model.id)
    model.id = bigg_id
    model.description = bigg_id

    old_ids = {'metabolites': metabolite_id_dict,
               'reactions': reaction_id_dict,
               'genes': gene_id_dict}
    return model, old_ids

# the regex to separate the base id, the chirality ('_L') and the compartment ('_c')
reg_compartment = re.compile(r'(.*?)[_\(\[]([a-z][a-z0-9]?)[_\)\]]?$')
reg_chirality = re.compile(r'(.*?)_?_([LDSRM])$')
def id_for_new_id_style(old_id, is_metabolite=False):
    """ Get the new style id"""
    new_id = old_id

    def _join_parts(the_id, the_compartment):
        if the_compartment:
            the_id = the_id + '_' + the_compartment
        return the_id

    def _remove_d_underscore(s):
        """Removed repeated, leading, and trailing underscores."""
        s = re.sub(r'_+', '_', s)
        s = re.sub(r'^_+', '', s)
        s = re.sub(r'_+$', '', s)
        return s

    # remove parentheses and brackets, for SBML & BiGG spec compatibility
    new_id = re.sub(r'[^a-zA-Z0-9_]', '_', new_id)

    # strip leading and trailing underscores
    # new_id = re.sub(r'^_+', '', new_id)
    # new_id = re.sub(r'_+$', '', new_id)

    compartment_match = reg_compartment.match(new_id)
    if compartment_match is None:
        # still remove double underscores
        new_id = _remove_d_underscore(new_id)
    else:
        base, compartment = compartment_match.groups()
        chirality_match = reg_chirality.match(base)
        if chirality_match is None:
            new_id = _join_parts(_remove_d_underscore(base), compartment)
        else:
            new_base = '%s__%s' % (_remove_d_underscore(chirality_match.group(1)),
                                   chirality_match.group(2))
            new_id = _join_parts(new_base, compartment)

    return new_id


def get_formulas_from_names(model):
    reg = re.compile(r'.*_([A-Za-z0-9]+)$')
    # support cobra 0.3 and 0.4
    for metabolite in model.metabolites:
        if (metabolite.formula is not None and str(metabolite.formula) != '' and getattr(metabolite, 'formula', None) is not None):
            continue
        m = reg.match(metabolite.name)
        if m:
            try:
                metabolite.formula = Formula(m.group(1))
            except TypeError:
                metabolite.formula = str(m.group(1))
    return model


def setup_model(model, substrate_reactions, aerobic=True, sur=10, max_our=10,
                id_style='cobrapy', fix_iJO1366=False):
    """Set up the model with environmntal parameters.

    model: a cobra model
    substrate_reactions: A single reaction id, list of reaction ids, or dictionary with reaction
    ids as keys and max substrate uptakes as keys. If a list or single id is
    given, then each substrate will be limited to /sur/
    aerobic: True or False
    sur: substrate uptake rate. Ignored if substrate_reactions is a dictionary.
    max_our: Max oxygen uptake rate.
    id_style: 'cobrapy' or 'simpheny'.

    """
    if id_style=='cobrapy': o2 = 'EX_o2_e'
    elif id_style=='simpheny': o2 = 'EX_o2(e)'
    else: raise Exception('Invalid id_style')

    if isinstance(substrate_reactions, dict):
        for r, v in substrate_reactions.iteritems():
            model.reactions.get_by_id(r).lower_bound = -abs(v)
    elif isinstance(substrate_reactions, list):
        for r in substrate_reactions:
            model.reactions.get_by_id(r).lower_bound = -abs(sur)
    elif isinstance(substrate_reactions, str):
        model.reactions.get_by_id(substrate_reactions).lower_bound = -abs(sur)
    else: raise Exception('bad substrate_reactions argument')

    if aerobic:
        model.reactions.get_by_id(o2).lower_bound = -abs(max_our)
    else:
        model.reactions.get_by_id(o2).lower_bound = 0

    # model specific setup
    if str(model)=='iJO1366' and aerobic==False:
        for r in ['CAT', 'SPODM', 'SPODMpp']:
            model.reactions.get_by_id(r).lower_bound = 0
            model.reactions.get_by_id(r).upper_bound = 0
    if fix_iJO1366 and str(model)=='iJO1366':
        for r in ['ACACT2r']:
            model.reactions.get_by_id(r).upper_bound = 0
        print 'made ACACT2r irreversible'

    # TODO hydrogen reaction for ijo

    if str(model)=='iMM904' and aerobic==False:
        necessary_ex = ['EX_ergst(e)', 'EX_zymst(e)', 'EX_hdcea(e)',
                        'EX_ocdca(e)', 'EX_ocdcea(e)', 'EX_ocdcya(e)']
        for r in necessary_ex:
            rxn = model.reactions.get_by_id(r)
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000

    return model

def turn_on_subsystem(model, subsytem):
    raise NotImplementedError()
    for reaction in model.reactions:
        if reaction.subsystem.strip('_') == subsytem.strip('_'):
            reaction.lower_bound = -1000 if reaction.reversibility else 0
            reaction.upper_bound = 1000
    return model

def carbons_for_exchange_reaction(reaction):
    if len(reaction._metabolites) > 1:
        raise Exception('%s not an exchange reaction' % str(reaction))

    metabolite = reaction._metabolites.iterkeys().next()
    try:
        return metabolite.formula.elements['C']
    except KeyError:
        return 0
    # match = re.match(r'C([0-9]+)', str(metabolite.formula))
    # try:
    #     return int(match.group(1))
    # except AttributeError:
    #     return 0

def fix_legacy_id(id, use_hyphens=False):
    id = id.replace('_DASH_', '__')
    id = id.replace('_FSLASH_', '/')
    id = id.replace('_BSLASH_', "\\")
    id = id.replace('_LPAREN_', '(')
    id = id.replace('_LSQBKT_', '[')
    id = id.replace('_RSQBKT_', ']')
    id = id.replace('_RPAREN_', ')')
    id = id.replace('_COMMA_', ',')
    id = id.replace('_PERIOD_', '.')
    id = id.replace('_APOS_', "'")
    id = id.replace('&amp;', '&')
    id = id.replace('&lt;', '<')
    id = id.replace('&gt;', '>')
    id = id.replace('&quot;', '"')
    if use_hyphens:
        id = id.replace('__', '-')
    else:
        id = id.replace("-", "__")
    return id

def split_compartment(component_id):
    """Split the metabolite bigg_id into a metabolite and a compartment id.

    Arguments
    ---------

    component_id: the bigg_id of the metabolite.

    """
    match = re.search(r'_[a-z][a-z0-9]?$', component_id)
    if match is None:
        raise NotFoundError("No compartment found for %s" % component_id)
    met = component_id[0:match.start()]
    compartment = component_id[match.start()+1:]
    return met, compartment

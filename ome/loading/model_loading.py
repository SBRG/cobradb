# -*- coding: utf-8 -*-
# This code is primarily a merger of theseus.models by @zakandrewking and LoadTheseus by @jslu9
from ome import base, settings, components, timing
from ome.models import *
from ome.loading import component_loading


import cobra
import cobra.io
from cobra.core.Formula import Formula
import os
from os.path import join, abspath, dirname
import re
import cPickle as pickle

from sqlalchemy import create_engine, Table, MetaData, update


data_path = settings.data_directory


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

# the regex to separate the base id, the chirality ('_L') and the compartment ('_c')
reg = re.compile(r'(.*?)(?:(.*[^_])_([LDSR]))?[_\(\[]([a-z])[_\)\]]?$')
def id_for_new_id_style(old_id, is_metabolite=False, new_id_style='cobrapy'):
    """ Get the new style id"""

    def join_parts(the_id, the_compartment):
        if (new_id_style.lower()=='cobrapy'):
            if the_compartment:
                the_id = the_id+'_'+the_compartment
            the_id = the_id.replace('-', '__')
        elif (new_id_style.lower()=='simpheny'):
            if the_compartment and is_metabolite:
                the_id = the_id+'['+the_compartment+']'
            elif the_compartment:
                the_id = the_id+'('+the_compartment+')'
            the_id = the_id.replace('__', '-')
        else:
            raise Exception('Invalid id style')
        return the_id

    # separate the base id, the chirality ('_L') and the compartment ('_c')
    m = reg.match(old_id)
    if m is None:
        # still change the underscore/dash
        new_id = join_parts(old_id, None)
    elif m.group(2) is None:
        new_id = join_parts(m.group(1), m.group(4))
    else:
        # if the chirality is not joined by two underscores, then fix that
        a = "__".join(m.groups()[1:3])
        new_id = join_parts(a, m.group(4))

    # deal with inconsistent notation of (sec) vs. [sec] in iJO1366 versions
    new_id = new_id.replace('[sec]', '_sec_').replace('(sec)', '_sec_')

    return new_id

def convert_ids(model, new_id_style):
    """Converts metabolite and reaction ids to the new style. Style options:

    cobrapy: EX_lac__L_e
    simpheny: EX_lac-L(e)

    """
    # loop through the ids:

    # this code comes from cobra.io.sbml
    # legacy_ids add special characters to the names again
    for metabolite in model.metabolites:
        metabolite.id = fix_legacy_id(metabolite.id, use_hyphens=False)
    model.metabolites._generate_index()
    for reaction in model.reactions:
        reaction.id = fix_legacy_id(reaction.id, use_hyphens=False)
    model.reactions._generate_index()
    # remove boundary metabolites (end in _b and only present in exchanges) . Be
    # sure to loop through a static list of ids so the list does not get
    # shorter as the metabolites are deleted
    for metabolite_id in [str(x) for x in model.metabolites]:
        metabolite = model.metabolites.get_by_id(metabolite_id)
        if not metabolite.id.endswith("_b"):
            continue
        for reaction in list(metabolite._reaction):
            if reaction.id.startswith("EX_"):
                metabolite.remove_from_model()
                break
    model.metabolites._generate_index()

    # separate ids and compartments, and convert to the new_id_style
    for reaction in model.reactions:
        reaction.id = id_for_new_id_style(reaction.id, new_id_style=new_id_style)
    model.reactions._generate_index()
    for metabolite in model.metabolites:
        metabolite.id = id_for_new_id_style(metabolite.id, is_metabolite=True, new_id_style=new_id_style)
    model.metabolites._generate_index()

    return model

def parse_model(name, id_style='cobrapy'):
    """Load a model, and give it a particular id style"""

    # check for model
    name = check_for_model(name)
    if not name:
        raise Exception('Could not find model')

    # load the model pickle, or, if not, the sbml
    try:
        with open(join(settings.data_directory, 'models', 'model_pickles', name+'.pickle'), 'r') as f:
            model = pickle.load(f)
    except:
        try:
            model = cobra.io.load_matlab_model(join(settings.data_directory, 'models', name+'.mat'))
        except:
            model = cobra.io.read_sbml_model(join(settings.data_directory, 'models', name+'.xml'))
        with open(join(settings.data_directory, 'models', 'model_pickles', name+'.pickle'), 'w') as f:
            pickle.dump(model, f)

    # convert the ids
    model = convert_ids(model, id_style)

    # extract metabolite formulas from names (e.g. for iAF1260)
    model = get_formulas_from_names(model)

    # turn off carbon sources
    model = turn_off_carbon_sources(model)

    return model

def get_formulas_from_names(model):
    reg = re.compile(r'.*_([A-Za-z0-9]+)$')
    for metabolite in model.metabolites:
        if (metabolite.formula is not None
            and metabolite.formula.formula!=''
            and metabolite.formula.formula is not None): continue
        m = reg.match(metabolite.name)
        if m:
            metabolite.formula = Formula(m.group(1))
    return model

def turn_off_carbon_sources(model):
    for reaction in model.reactions:
        if 'EX_' not in str(reaction): continue
        if carbons_for_exchange_reaction(reaction) > 0:
            reaction.lower_bound = 0
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

def add_pathway(model, new_metabolites, new_reactions, subsystems, bounds,
                check_mass_balance=False, ignore_repeats=False):
    """Add a pathway to the model. Reversibility defaults to reversible (1).

    new_metabolites: e.g. { 'ggpp_c': {'formula': 'C20H33O7P2', 'name': 'name'},
                            'phyto_c': {'formula': 'C40H64'}},
                            'lyco_c': {'formula': 'C40H56'},
                            'lyco_e': {'formula': 'C40H56'} }
    new_reactions: e.g. { 'FPS': { 'ipdp_c': -2,
                                   'ppi_c': 1,
                                   'grdp_c': 1 },
                          'CRTE': { 'ipdp_c': -1,
                                    'frdp_c': -1,
                                    'ggpp_c': 1,
                                    'ppi_c': 1 } }
    subsystems: e.g. { 'FPS': 'Lycopene production',
                       'CRTE': 'Lycopene production' }
    bound: e.g. { 'FPS': (0, 0),
                  'CRTE': (0, 1000) }

    """

    for k, v in new_metabolites.iteritems():
        formula = Formula(v['formula']) if 'formula' in v else None
        name = v['name'] if 'name' in v else None
        m = cobra.Metabolite(id=k, formula=formula, name=name)
        try:
            model.add_metabolites([m])
        except Exception as err:
            if (not ignore_repeats or
                "already in the model" not in str(err)):
                raise(err)

    for name, mets in new_reactions.iteritems():
        r = cobra.Reaction(name=name)
        m_obj = {}
        for k, v in mets.iteritems():
            m_obj[model.metabolites.get_by_id(k)] = v
        r.add_metabolites(m_obj)
        if bounds and (name in bounds):
            r.lower_bound, r.upper_bound = bounds[name]
        else:
            r.upper_bound = 1000
            r.lower_bound = -1000
        if subsystems and (name in subsystems):
            r.subsystem = subsystems[name]
        try:
            model.add_reaction(r)
        except Exception as err:
            if (not ignore_repeats or
                "already in the model" not in str(err)):
                raise(err)
        if check_mass_balance and 'EX_' not in name:
            balance = model.reactions.get_by_id(name).check_mass_balance()
            if balance != []:
                raise Exception('Bad balance: %s' % str(balance))
    return model

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


class IndependentObjects:

    def loadGenes(self, modellist, session):
        for model in modellist:
            for gene in model.genes:
                if not session.query(Gene).filter(Gene.name == gene.id).count():
                    geneObject = Gene(locus_id = gene.id)
                    session.add(geneObject)

    def loadModel(self, model, session, genome_id, first_created):
        modelObject = Model(bigg_id = model.id, first_created = first_created, genome_id = genome_id, notes = '')
        session.add(modelObject)


        """
            if(model.id == 'iSF1195'):
                modelObject = Model(bigg_id = model.id, firstcreated = '2014-9-16 14:26:22', genome_id = 7)
            if(model.id == 'iSB619'):
                modelObject = Model(bigg_id = model.id, firstcreated = '2013-10-21 14:26:22', genome_id = 3)
            if(model.id == 'iJN746'):
                modelObject = Model(bigg_id = model.id, firstcreated = '2013-10-21 14:26:22', genome_id = 6)
            if(model.id == 'iIT341'):
                modelObject = Model(bigg_id = model.id, firstcreated = '2013-10-21 14:26:22', genome_id = 4)
            if(model.id == 'iNJ661'):
                modelObject = Model(bigg_id = model.id, firstcreated = '2013-10-21 14:26:22', genome_id = 3)
            if(model.id == 'iJO1366'):
                modelObject = Model(bigg_id = model.id, firstcreated = '2013-10-21 14:26:22', genome_id = 15)
            if(model.id == 'iAF692'):
                modelObject = Model(bigg_id = model.id, firstcreated = '2013-10-21 14:26:22', genome_id = 5)
            if(model.id == 'model'):
                modelObject = Model(bigg_id = model.id, firstcreated = '2013-10-21 14:26:22', genome_id = 1)
            if(model.id == 'iAPECO1_1312 '):

            session.add(modelObject)
        """

    def loadComponents(self, modellist, session):
        for model in modellist:
            for component in model.metabolites:
                metabolite = session.query(Metabolite).filter(Metabolite.name == component.id.split("_")[0])
                #metabolite = session.query(Metabolite).filter(Metabolite.kegg_id == component.notes.get("KEGGID")[0])
                if not metabolite.count():
                    try: kegg_id = component.notes.get("KEGGID")[0]
                    except: kegg_id = None
                    try: cas_number = component.notes.get("CASNUMBER")[0]
                    except: cas_number = None
                    try: formula = component.notes.get("FORMULA1")[0]
                    except: formula = None

                    metaboliteObject = Metabolite(name = component.id.split("_")[0],
                                                  long_name = component.name,
                                                  kegg_id = kegg_id,
                                                  cas_number = cas_number,
                                                  formula = formula,
                                                  flag = bool(kegg_id))

                    session.add(metaboliteObject)
                else:
                    metaboliteObject = metabolite.first()
                    if metaboliteObject.kegg_id == None or metaboliteObject.kegg_id == '':
                        if 'KEGGID' in component.notes.keys():
                            metaboliteObject.kegg_id = component.notes.get("KEGGID")[0]
                        #metabolite.update({Metabolite.kegg_id: str(component.notes.get("KEGGID"))})
                    if metaboliteObject.cas_number == None or metaboliteObject.cas_number == '':
                        if 'CASNUMBER' in component.notes.keys():
                            metaboliteObject.cas_number = component.notes.get("CASNUMBER")[0]
                        #metabolite.update({Metabolite.cas_number: str(component.notes.get("CASNUMBER"))})
                    if metaboliteObject.formula == None or metaboliteObject.formula == '':
                        if 'FORMULA1' in component.notes.keys():
                            metaboliteObject.formula = component.notes.get("FORMULA1")[0]
                        #metabolite.update({Metabolite.formula: str(component.notes.get("FORMULA1"))})

    def loadReactions(self , modellist, session):
        for model in modellist:
            for reaction in model.reactions:
                if not session.query(Reaction).filter(Reaction.name == reaction.id).count():
                    reactionObject = Reaction(name = reaction.id, long_name = reaction.name, notes = '')
                    session.add(reactionObject)

    def loadCompartments(self, modellist, session):
        for model in modellist:
            for component in model.metabolites:
                if component.id is not None:
                    if not session.query(Compartment).filter(Compartment.name == component.id[-1:len(component.id)]).count():
                        compartmentObject = Compartment(name = component.id[-1:len(component.id)])
                        session.add(compartmentObject)



class DependentObjects:
    def loadModelGenes(self, modellist, session):
        for model in modellist:
            for gene in model.genes:
                if gene.id != 's0001':
                    if session.query(Gene).filter(Gene.locus_id == gene.id).first() != None:
                        genequery = session.query(Gene).filter(Gene.locus_id == gene.id).first()
                        modelquery = session.query(Model).filter(Model.bigg_id == model.id).first()
                        #genequery = session.query(Gene).filter(Gene.locus_id == gene.id).filter(Gene.genome_id == modelquery.genome_id).first()
                        object = Model_Gene(model_id = modelquery.id, gene_id = genequery.id)
                        session.add(object)
                    elif session.query(Gene).filter(Gene.name == gene.id).first() != None:
                        genequery = session.query(Gene).filter(Gene.name == gene.id).first()
                        modelquery = session.query(Model).filter(Model.bigg_id == model.id).first()
                        object = Model_Gene(model_id = modelquery.id, gene_id = genequery.id)
                        session.add(object)
                    else:
                        #geneObject = Gene(locus_id = gene.id, leftpos=None, rightpos=None, strand=None, name=gene.id)
                        #session.add(geneObject)
                        synonymquery = session.query(Synonyms).filter(Synonyms.synonym == gene.id.split(".")[0]).filter(Synonyms.type == 'gene').first()
                        if synonymquery != None:
                            modelquery = session.query(Model).filter(Model.bigg_id == model.id).first()

                            genecheck = session.query(Gene).filter(Gene.id == synonymquery.ome_id).first()
                            if genecheck:
                                object = Model_Gene(model_id = modelquery.id, gene_id = synonymquery.ome_id)
                                session.add(object)

                                if modelquery.bigg_id == "RECON1":
                                    genequery = session.query(Gene).filter(Gene.id == synonymquery.ome_id).first()
                                    genequery.locus_id = gene.id
                            else:
                                print synonymquery.ome_id
                        else:
                            print gene.id, model.id


    def loadCompartmentalizedComponent(self, modellist, session):
        for model in modellist:
            for metabolite in model.metabolites:
                identifier = session.query(Compartment).filter(Compartment.name == metabolite.id[-1:len(metabolite.id)]).first()
                m = session.query(Metabolite).filter(Metabolite.name == metabolite.id.split("_")[0]).first()
                #m = session.query(Metabolite).filter(Metabolite.kegg_id == metabolite.notes.get("KEGGID")[0]).first()
                object = Compartmentalized_Component(component_id = m.id, compartment_id = identifier.id)
                session.add(object)

    def loadModelCompartmentalizedComponent(self, modellist, session):
        for model in modellist:
            for metabolite in model.metabolites:
                componentquery = session.query(Metabolite).filter(Metabolite.name == metabolite.id.split("_")[0]).first()
                #componentquery = session.query(Metabolite).filter(Metabolite.kegg_id == metabolite.notes.get("KEGGID")[0]).first()
                compartmentquery = session.query(Compartment).filter(Compartment.name == metabolite.id[-1:len(metabolite.id)]).first()
                compartmentalized_component_query = session.query(Compartmentalized_Component).filter(Compartmentalized_Component.component_id == componentquery.id).filter(Compartmentalized_Component.compartment_id == compartmentquery.id).first()
                modelquery = session.query(Model).filter(Model.bigg_id == model.id).first()
                if modelquery is None:
                    print "model query is none", model.id
                    from IPython import embed; embed()

                if compartmentalized_component_query is None:
                    print "compartmentalized_component_query is none", metabolite.id
                object = Model_Compartmentalized_Component(model_id = modelquery.id, compartmentalized_component_id = compartmentalized_component_query.id, compartment_id = compartmentquery.id)
                session.add(object)


    def loadModelReaction(self, modellist, session):
        for model in modellist:
            for reaction in model.reactions:
                reactionquery = session.query(Reaction).filter(Reaction.name == reaction.id).first()
                modelquery = session.query(Model).filter(Model.bigg_id == model.id).first()
                object = Model_Reaction(reaction_id = reactionquery.id, model_id = modelquery.id, name = reaction.id, upperbound = reaction.upper_bound, lowerbound = reaction.lower_bound, gpr = reaction.gene_reaction_rule)
                session.add(object)


    def loadGPRMatrix(self, modellist, session):
        for model in modellist:
            for reaction in model.reactions:
                for gene in reaction._genes:
                    if gene.id != 's0001':

                        model_query = session.query(Model).filter(Model.bigg_id == model.id).first()
                        model_gene_query = session.query(Model_Gene).join(Gene).filter(Gene.locus_id == gene.id).filter(Model_Gene.model_id == model_query.id).first()

                        if model_gene_query != None:
                            model_reaction_query = session.query(Model_Reaction).filter(Model_Reaction.name == reaction.id).filter(Model_Reaction.model_id == model_query.id).first()
                            object = GPR_Matrix(model_gene_id = model_gene_query.id, model_reaction_id = model_reaction_query.id)
                            session.add(object)
                        else:
                            model_gene_query = session.query(Model_Gene).join(Gene).filter(Gene.name == gene.id).filter(Model_Gene.model_id == model_query.id).first()
                            if model_gene_query != None:
                                model_reaction_query = session.query(Model_Reaction).filter(Model_Reaction.name == reaction.id).filter(Model_Reaction.model_id == model_query.id).first()
                                object = GPR_Matrix(model_gene_id = model_gene_query.id, model_reaction_id = model_reaction_query.id)
                                session.add(object)
                            else:
                                synonymquery = session.query(Synonyms).filter(Synonyms.synonym == gene.id.split(".")[0]).first()
                                if synonymquery != None:
                                    if synonymquery.ome_id != None:
                                        model_gene_query = session.query(Model_Gene).join(Gene).filter(Gene.id == synonymquery.ome_id).filter(Model_Gene.model_id == model_query.id).first()
                                        model_reaction_query = session.query(Model_Reaction).filter(Model_Reaction.name == reaction.id).filter(Model_Reaction.model_id == model_query.id).first()
                                        object = GPR_Matrix(model_gene_id = model_gene_query.id, model_reaction_id = model_reaction_query.id)
                                        session.add(object)
                                    else:
                                        print "ome id is null " + synonymquery.ome_id
                                else:
                                    print "mistake", gene.id, reaction.id

    def loadReactionMatrix(self, modellist, session):
        for model in modellist:
            for reaction in model.reactions:
                reactionquery = session.query(Reaction).filter(Reaction.name == reaction.id).first()
                for metabolite in reaction._metabolites:

                    componentquery = session.query(Metabolite).filter(Metabolite.name == metabolite.id.split("_")[0]).first()
                    #componentquery = session.query(Metabolite).filter(Metabolite.kegg_id == metabolite.notes.get("KEGGID")[0]).first()
                    compartmentquery = session.query(Compartment).filter(Compartment.name == metabolite.id[-1:len(metabolite.id)]).first()
                    compartmentalized_component_query = session.query(Compartmentalized_Component).filter(Compartmentalized_Component.component_id == componentquery.id).filter(Compartmentalized_Component.compartment_id == compartmentquery.id).first()
                    if not session.query(Reaction_Matrix).filter(Reaction_Matrix.reaction_id == reactionquery.id).filter(Reaction_Matrix.compartmentalized_component_id == compartmentalized_component_query.id).count():
                        for stoichKey in reaction._metabolites.keys():
                            if str(stoichKey) == metabolite.id:
                                stoichiometryobject = reaction._metabolites[stoichKey]
                        object = Reaction_Matrix(reaction_id = reactionquery.id, compartmentalized_component_id = compartmentalized_component_query.id, stoichiometry = stoichiometryobject)
                        session.add(object)

    def loadEscher(self, session):
        m = models.parse_model('iJO1366')
        for reaction in m.reactions:
            escher = Escher_Map(bigg_id = reaction.id, category = "reaction", model_name = m.id)
            session.add(escher)




@timing
def load_model(model_id, genome_id, model_creation_timestamp):
    with create_Session() as session:

        try: genome = session.query(base.Genome).filter_by(bioproject_id=genome_id).one()
        except:
            print 'Genbank file %s for model %s was not uploaded' % (genome_id, model_id)
            return

        model = parse_model(model_id)

        IndependentObjects().loadModel(model, session, genome.id, model_creation_timestamp)
        IndependentObjects().loadComponents([model], session)
        IndependentObjects().loadCompartments([model], session)
        DependentObjects().loadCompartmentalizedComponent([model], session)
        IndependentObjects().loadReactions([model], session)
        DependentObjects().loadModelGenes([model], session)
        DependentObjects().loadModelCompartmentalizedComponent([model], session)
        DependentObjects().loadModelReaction([model], session)
        DependentObjects().loadGPRMatrix([model], session)
        DependentObjects().loadReactionMatrix([model], session)
        #DependentObjects().loadEscher(session)




# -*- coding: utf-8 -*-

"""
@author: P0N3Y (Pierre SALVY)
Submitted by Amyris (Amyris, Inc). for SBRG@UCSD
"""

import sqlalchemy
import cobra
import cobra.io
from cobra.core.Formula import Formula
import os
from os.path import join, abspath, dirname
import re
import cPickle

from sqlalchemy import create_engine, Table, MetaData, update
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Float
from sqlalchemy.orm import sessionmaker
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref, joinedload
from sqlalchemy.orm import aliased
from sqlalchemy.sql import func

from ome import base, settings, components, timing
#from ome.base import *
from ome.browsing_tools import make_dict, make_light_dict, get_model_or_latest
from ome.models import *
from ome.components import *
from ome.loading import component_loading
from ome import gpr_rpn_tree as GPRTools

'''
This is the module that enables you to pull stuff from the DB in an efficient manner. Hopefully.
There are were originally modes :
[#] gossip_girls == querying everything to the db when needed
[#] crazy_hermit == fetches everything and reconstructs afterwards from
local copies.
The latter is much more efficient, with a loading time around 14s (vs 135s
for the former) for a 1000 mets/1200 rxns  model. crazy hermit is the only implementation you'll find here.
I hope the function names are obvious enough to let you guess what they do.
(Like for instance feed_the_cat() is likely not going to make you some coffee)
'''

class PullModel:

    def __init__(self, session, modelObject):

        self.elements = make_light_dict(engine.execute(select([Element.__table__.c.id, Element.__table__.c.symbol])))
        self.component_chemical_formulas = make_dict(session.query(ComponentChemicalFormula))
        self.compartments = make_dict(session.query(Compartment))
        self.metabolites = make_dict(session.query(Metabolite))
        self.reactions = make_dict(session.query(Reaction).\
            join(ModelReaction, ModelReaction.reaction_id == Reaction.id).\
            filter(ModelReaction.model_id == modelObject.id))
        self.reaction_matrixes = make_dict(session.query(ReactionMatrix))
        self.gpr_trees = make_dict(session.query(GPRTree))
        self.model_reactions = make_dict(session.query(ModelReaction).\
            filter(ModelReaction.model_id == modelObject.id))
        self.compartmentalized_components = make_dict(session.query(CompartmentalizedComponent).\
            #join(ReactionMatrix, CompartmentalizedComponent.id == \
            #    ReactionMatrix.compartmentalized_component_id).\
            #join(Reaction, Reaction.id == ReactionMatrix.reaction_id).\
            #join(ModelReaction, ModelReaction.reaction_id == Reaction.id).\
            #filter(ModelReaction.model_id == modelObject.id))
            join(ModelCompartmentalizedComponent, CompartmentalizedComponent.id == \
                ModelCompartmentalizedComponent.compartmentalized_component_id).\
            filter(ModelCompartmentalizedComponent.model_id == modelObject.id))
        self.genes = make_dict(session.query(Gene).\
            join(ModelGene, Gene.id == ModelGene.gene_id).\
            filter(ModelGene.model_id == modelObject.id))


        self.chem_dict = {'C':0, 'H':1, 'O':2, 'N':3, 'P':4, 'S':5}
        self.chem_order = lambda x:self.chem_dict[x] if self.chem_dict.has_key(x) else 6

        def make_full_formula(r):
            '''
            compiles the full formula in a nice CHONPS ordering.
            '''
            try:
                db_style_f = formula_matrixes_by_met.loc[r]
                f_dict = {}
                for bit in db_style_f.iterrows():
                    count = int(bit[1]['count']) if float.is_integer(bit[1]['count']) else bit[1]['count']
                    symbol = elements_table.loc[int(bit[1].element_id)].symbol
                    bit = symbol + ('<sub>' + str(count) + '</sub>') if count > 1 else ''
                    f_dict[bit] = chem_order(symbol)
                imfedupwiththis = sorted(f_dict.items(), key = lambda x:x[1])
                formula = ''.join([x[0] for x in imfedupwiththis])
            except Exception as e:
                #print e
                formula = ''


    def make_formula(self, fbits):
        f = ''
        f_dict = {}
        for x in fbits:
            count = int(x.count) if float.is_integer(x.count) else x.count
            symbol = self.elements[x.element_id][1]
            bit = symbol + (str(count) if count > 1 else '')
            f_dict[bit] = self.chem_order(symbol)
        sorted_f = sorted(f_dict.items(), key = lambda x:x[1])
        f = ''.join([x[0] for x in sorted_f])
        '''
            #elementObject = elements[x.element_id]
            c = x.count
            if c.is_integer():
                c = int(c)
            #f = f + str(elementObject.symbol) + str(c)
            f = f + str(self.elements[x.element_id][1]) + str(c)
        '''
        return str(f)


    def make_metabolite(self, compartmentalized_component_id):
        ccObject = self.compartmentalized_components[compartmentalized_component_id]
        compartmentObject = self.compartments[ccObject.compartment_id]
        metaboliteObject = self.metabolites[ccObject.component_id]
        formula_bits = [x for x in self.component_chemical_formulas.itervalues() if x.component_id == metaboliteObject.id]
        formula = self.make_formula(formula_bits)
        mid = metaboliteObject.name + '_' + compartmentObject.symbol
        name = metaboliteObject.long_name + '|' + metaboliteObject.short_name \
            if metaboliteObject.short_name \
            else metaboliteObject.long_name
        cobra_metabolite = cobra.core.Metabolite(
            id = mid,
            name = name,
            formula = formula,
            compartment = compartmentObject.name)
        cobra_metabolite.notes={}
        cobra_metabolite.notes['CHEBI'] = metaboliteObject.chebi
        return cobra_metabolite

    def make_all_metabolites(self):
        self.cobra_metabolites = {}
        for cc in self.compartmentalized_components:
            self.cobra_metabolites[cc] = self.make_metabolite(cc)


    def pull_rxn_equation(self, rxn_id):
        db_style_equation = [x for x in self.reaction_matrixes.itervalues() if x.reaction_id == rxn_id]
        d = {}
        for x in db_style_equation:
            m = self.cobra_metabolites[x.compartmentalized_component_id]
            d[m] = x.stoichiometry
        return d

    def pull_gpr(self, gpr_id, standard_names = True):
        gpr_trees = self.gpr_trees
        genes = self.genes
        try:
            tree_query = gpr_trees[gpr_id]
        except KeyError:
            tree_query = False
        node_value = ""
        if tree_query:
            treeObject = tree_query
            if treeObject.operator == 'leaf':
                geneObject = genes[treeObject.left]
                if standard_names and geneObject.name:
                    ret = geneObject.name #standard name
                else:
                    ret = geneObject.locus_id
                return ret
            else:
                node_value = self.pull_gpr(treeObject.left, standard_names)
                node_value = node_value + ' ' + str(treeObject.operator)
                node_value = node_value + ' ' + self.pull_gpr(treeObject.right, standard_names)
                node_value = '( ' + node_value + ' )'
                return node_value



    def pull_rxns(self, model_id, standard_names = True):
        rxns_list = [x for x in self.model_reactions.itervalues() if x.model_id == model_id]
        authorized = [x for x in self.model_reactions.itervalues()]
        cobra_list = []
        self.make_all_metabolites()
        for r in rxns_list:
            try:
                if authorized and r not in authorized:
                    raise Exception
                rxnObject = self.reactions[r.reaction_id]
                rdict = self.pull_rxn_equation(rxnObject.id)
                name = rxnObject.long_name
                rid = rxnObject.name
                cobra_rxn = cobra.core.Reaction(name = name)
                cobra_rxn.id = rid
                cobra_rxn.lower_bound = r.lowerbound
                cobra_rxn.upper_bound = r.upperbound
                cobra_rxn.add_metabolites(rdict)
                gpr = self.pull_gpr(r.gpr_id, standard_names)
                if gpr:
                    cobra_rxn.gene_reaction_rule = gpr
                cobra_list.append(cobra_rxn)
            except Exception as e:
                pass
                #print e
                #print 'Reaction (#%i) - %s not selected in the model'%\
                #    (r.reaction_id, self.reactions[r.reaction_id].name)


        return cobra_list

@timing
def pull_model(version, standard_names = True):

    engine = create_engine(db_connection_string, echo=False)
    Session = sessionmaker(bind=engine)
    session = Session()

    print "Pulling model using the", mode.upper(), "method"
    modelObject = get_model_or_latest(session, version)
    my_model = cobra.core.Model(version)


    #importing the reactions
    pull_request = PullModel(session, modelObject)
    rxns = pull_request.pull_rxns(standard_names)
    mets = pull_request.cobra_metabolites
    my_model.add_reactions(rxns)
    #to add the metabolites that are not in reactions
    my_model.add_metabolites(mets.values())

    return my_model

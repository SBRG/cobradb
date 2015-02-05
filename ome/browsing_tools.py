# -*- coding: utf-8 -*-

"""
@author: P0N3Y (Pierre SALVY)
Submitted by Amyris (Amyris, Inc). for SBRG@UCSD
"""


'''
This module aims at giving pre-written functions that perform tasks on a
kind of optimized way given the organisation of the DB. Really, it's 
just so I don't have to write always the same lines to test stuff. Also,
writing queries with SA is fairly intuitive so you can give it a try.
'''


from sqlalchemy import Table, MetaData, create_engine
from sqlalchemy.orm import sessionmaker

from ome.base import *
from ome.models import *
from settings import db_connection_string


def make_dict(statement, indexing_attributes = [], core_friendly = False):
    '''
    Super cool function that transforms an ORM query into a dict indexed by
    a chosen indexing_attribute. core_friendly returns a dict instead
    of the object, useful to just update data for batch inserts
    '''
    l = statement.all()
    d = {}
    if indexing_attributes:
        if isinstance(indexing_attributes, str):
            for r in l:
                d[getattr(r, indexing_attributes)] = r if not core_friendly else r.__dict__
        else:
            for r in l:
                index = tuple([getattr(r, i) for i in indexing_attributes])
                d[index] = r if not core_friendly else r.__dict__
    else:
        for r in l:
            d[r.id] = r if not core_friendly else r.__dict__

        #Crappy DF solution, a.k.a how to kill a fly with a bazooka
        '''
        d = DataFrame({'object':l})
        if isinstance(indexing_attributes, str):
                indexing_attributes = [indexing_attributes]
            if indexing_attribute:
                for i in indexing_attributes
                    d[i] = d['object'].apply(lambda r:getattr(r, i))
                d.set_index(indexing_attributes, inplace = True)
            else:
                d['id'] = d['object'].apply(lambda r:getattr(r, i))
                d.set_index(['id'], inplace = True)
        '''    
    return d

def make_light_dict(statement, indexing_attribute = ''):
    '''
    Same but using SQLA Core instead of ORM
    '''
    if not indexing_attribute:
        indexing_attribute = 'id'
    l = statement.fetchall()
    d = {}
    for e in l:
        d[e[indexing_attribute]] = e
    return d

def get_the_right_gene(session, sinput, method = 'locus_id'):
    '''
    Pretty much a synonym resolution service that returns the gene object
    corresponding to a name found in litterature.
    '''
    systematic_query = session.query(Gene).filter(getattr(Gene, method) == sinput)
    if systematic_query.count():
        return systematic_query.first()
    else:
        joined_query = session.query(Gene).\
            join(Synonyms, Synonyms.ome_id == Gene.id).\
            filter(Synonyms.synonym == sinput).\
            filter(Synonyms.type == 'gene_standard_name')
        if joined_query.count():
            return joined_query.first()
        else:
            print "No synonym found for :", sinput

class GPRCompiler:
    def __init__(self, genes, gprtrees):
        self.genes = genes
        self.gprtrees = gprtrees

    def compile_gpr(self, gpr_id, method = 'genes'):
        '''
        Returns the GPR Tree as a good-looking string or list depending on the method
        '''
        #tree_query = session.query(GPRTree).filter(GPRTree.id == gpr_id)
        genes = self.genes
        gprtrees = self.gprtrees
        node_value = ""
        if gprtrees.has_key(gpr_id):
            treeObject = gprtrees[gpr_id]
            if treeObject['operator'] == 'leaf':
                geneObject = genes[treeObject['left']]
                if method == 'gpr': return geneObject['locus_id']
                elif method == 'genes': return [geneObject]
            else:
                if method == 'gpr':
                    node_value = self.compile_gpr(treeObject['left'], method = method)
                    node_value = node_value + ' ' + str(treeObject['operator'])
                    node_value = node_value + ' ' + self.compile_gpr(treeObject['right'], method = method)
                    node_value = '( ' + node_value + ' )'
                    return node_value
                elif method == 'genes':
                    gene_list = []
                    gene_list.extend(self.compile_gpr(treeObject['left'], method = method))
                    gene_list.extend(self.compile_gpr(treeObject['right'], method = method))
                    return gene_list


def get_gpr_from_reaction(session, rxn_name, model_id, method = 'genes'):
    '''
    Returns the compiled GPR of a reaction by its name. Like 'exchange_of_stuff', or
    'magical_transhydrogenase', or 'something_synthase_badass_mitochondrial_call_me_maybe'
    '''
    modelObject = session.query(Model).filter(Model.bigg_id == model_id).first()
    cache = cache_data(session, ['genes','gpr']) 
    genes = cache['genes']
    gprtrees = cache['gpr']

    rxn_query = session.query(ModelReaction).\
                join(Reaction, Reaction.id == ModelReaction.reaction_id).\
                join(ReactionMatrix, Reaction.id == ReactionMatrix.reaction_id).\
                filter(ModelReaction.model_id == modelObject.id).\
                filter(Reaction.name == rxn_name)
    if rxn_query.count():
        rxnObject = rxn_query.first()
        compiler = GPRCompiler(genes, gprtrees)
        return compiler.compile_gpr(rxnObject.gpr_id,gprtrees,genes,method)

def get_the_genes_behind_this_metabolite(session, metname, model_id = None, method = 'genes'):
    '''
    OK, this one was an original request. You want to know which genes have an impact on glutathione ?
    just type get_the_genes_behind_this_metabolite(session, 'glutathione', 'my_model') !
    You get, a list of gene objects and some AWESOME !
    Methods :
    (*) genes give you a list of gene objects
    (*) gpr gives you gpr strings
    '''
    model_query = session.query(Model).filter(Model.bigg_id == model_id)
    if not model_id or not model_query.count():
        model_query = session.query(Model).order_by(Model.id.desc())
        modelObject = model_query.first()
    model_reactions = make_dict(session.query(ModelReaction).\
                join(Reaction, Reaction.id == ModelReaction.reaction_id).\
                filter(ModelReaction.model_id == modelObject.id))
    reaction_matrixes = make_dict(session.query(ReactionMatrix).\
                join(Reaction, Reaction.id == ReactionMatrix.reaction_id).\
                filter(Reaction.id.in_(model_reactions)),\
                ['compartmentalized_component_id', 'id'])
    component_query = session.query(Component).filter(Component.name == metname)
    cache = cache_data(session, ['genes','gpr']) 
    genes = cache['genes']
    gprtrees = cache['gpr'] 

    if component_query.count():
        componentObject = component_query.first()
        ccList = session.query(CompartmentalizedComponent).\
            filter(CompartmentalizedComponent.component_id  == componentObject.id).all()
        rxnList = []
        for m in ccList:
            to_add = [model_reactions[x[1]] for x in reaction_matrixes if x[0] == m.id]
            rxnList.extend(to_add)
        glist = []
        compiler = GPRCompiler(genes, gprtrees)
        for r in rxnList:
            if method == 'gpr' :glist.append(compiler.compile_gpr( r.gpr_id, method))
            elif method == 'genes':
                item = compiler.compile_gpr( r.gpr_id, method)
                if item:
                    glist.extend(item)

        return glist


def get_all_genes_to_metabolites_groups(session, model_id = None):
    '''
    Faster implementation of the previous algorithm for all
    the metabolites in the DB.
    '''
    modelObject = get_model_or_latest(session, model_id)
    model_reactions = make_dict(session.query(ModelReaction).\
                join(Reaction, Reaction.id == ModelReaction.reaction_id).\
                filter(ModelReaction.model_id == modelObject.id), 'reaction_id')
    reaction_matrixes = make_dict(session.query(ReactionMatrix).\
                join(Reaction, Reaction.id == ReactionMatrix.reaction_id).\
                filter(Reaction.id.in_(model_reactions)),\
                ['compartmentalized_component_id', 'reaction_id'])
    compartments = make_dict(session.query(Compartment))
    metabolites = make_dict(session.query(Metabolite))
    compartmentalized_components = make_dict(session.query(CompartmentalizedComponent), \
        ['component_id','compartment_id'])
    cache = cache_data(session, ['genes','gpr']) 
    genes = cache['genes']
    gprtrees = cache['gpr']
    d = {}

    compiler = GPRCompiler(genes, gprtrees)

    for cc in compartmentalized_components:
        ccObject = compartmentalized_components[cc]
        rxnList = []
        to_add = [model_reactions[x[1]] for x in reaction_matrixes if x[0] == ccObject.id]
        rxnList.extend(to_add)
        cc_name = metabolites[cc[0]].name+'_'+compartments[cc[1]].symbol
        d[cc_name] = []
        for r in rxnList:
            item = compiler.compile_gpr( r.gpr_id, 'genes')
            if item:
                d[cc_name].extend(item)

    return d


def get_model_or_latest(session, model_id):
    '''
    Fetches the latest model object if the model_id is invalid
    '''
    model_query = session.query(Model).filter(Model.bigg_id == model_id)
    if not model_id or not model_query.count():
        print 'Model not found, fetching latest'
        model_query = session.query(Model).order_by(Model.id.desc())
    modelObject = model_query.first()
    return modelObject




if __name__ == '__main__':
    '''
    Duh.
    '''
    engine = create_engine(db_connection_string, echo=False)
    Session = sessionmaker(bind=engine)
    session = Session()

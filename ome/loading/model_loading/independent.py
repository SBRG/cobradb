# -*- coding: utf-8 -*-

from ome import base
from ome.components import *
from ome.models import *
from ome.loading.model_loading import parse, queries

import logging

def loadGenes(model_list, session):
    for model in model_list:
        for gene in model.genes:
            if not session.query(Gene).filter(Gene.bigg_id == gene.id).count():
                geneObject = Gene(bigg_id = gene.id)
                session.add(geneObject)

def loadModel(session, model, genome_db_id, first_created, pmid):
    modelObject = Model(bigg_id=model.id, first_created=first_created,
                        genome_id=genome_db_id, notes='', description = '')
    session.add(modelObject)
    if session.query(base.Publication).filter(base.Publication.pmid == pmid).count() == 0: 
        p = base.Publication(pmid = pmid)
        session.add(p)
        publication = session.query(base.Publication).filter(base.Publication.pmid == pmid).first()
        pm = base.PublicationModel(publication_id= publication.id, model_id = modelObject.id)
        session.add(pm)
    else:
        publication = session.query(base.Publication).filter(base.Publication.pmid == pmid).first()
        pm = base.PublicationModel(publication_id= publication.id, model_id = modelObject.id)
        session.add(pm)


def loadComponents(session, model_list):
    for model in model_list:
        for component in model.metabolites:
            try:
                met_id = parse.split_compartment(component.id)[0]
            except Exception as e:
                logging.warn("%s. In model %s" % (e, model.id))
                continue

            linkouts = [('KEGGID', 'kegg_id'), 
                        ('CAS_NUMBER', 'cas_number'), 
                        ('SEED', 'seed'),
                        ('METACYC', 'metacyc'),
                        ('CHEBI', 'chebi'),
                        ('BRENDA', 'brenda'),
                        ('UPA', 'upa')]
            def parse_linkout_str(id):
                if id is None:
                    return None
                id_string = str(id)
                for s in ['{', '}', '[', ']', '&apos;', "'",]:
                    id_string = id_string.replace(s, '')
                return id_string

            # if there is no metabolite, add a new one
            metabolite_db = (session
                             .query(Metabolite)
                             .filter(Metabolite.bigg_id == met_id)
                             .first())
            if metabolite_db is None:
                found = {}
                for linkout in linkouts:
                    found[linkout[0]] = parse_linkout_str(component.notes.get(linkout[0]))

                # look for the formula
                if component.notes.get("FORMULA1") is not None:
                    _formula = component.notes.get("FORMULA1")
                else:
                    _formula = component.formula

                metaboliteObject = Metabolite(bigg_id=met_id,
                                              name=component.name,
                                              kegg_id=found['KEGGID'],
                                              cas_number=found['CAS_NUMBER'],
                                              seed=found['SEED'], 
                                              chebi=found['CHEBI'], 
                                              metacyc=found['METACYC'],
                                              upa=found['UPA'], 
                                              brenda=found['BRENDA'],
                                              formula=str(_formula))
                session.add(metaboliteObject)
            else:
                for linkout in linkouts:
                    need_linkout = (component.notes.get(linkout[0]) is not None and
                                    getattr(metabolite_db, linkout[1]) is not None)
                    if need_linkout:
                        setattr(metabolite_db, linkout[1],
                                parse_linkout_str(component.notes.get(linkout[0])))
                            
def loadReactions(session, model_list):
    for model in model_list:
        for reaction in model.reactions:
            reaction_db = queries.get_reaction(session, reaction.id)
            if reaction_db is None:
                new_object = Reaction(bigg_id=reaction.id, name=reaction.name, # fix this, see BiGG2 issue #29
                                      notes='', reaction_hash=parse.hash_reaction(reaction))
                session.add(new_object)

def loadCompartments(session, model_list):
    compartments_all = set()
    for model in model_list:
        for component in model.metabolites:
            if component.id is not None:
                compartments_all.add(parse.split_compartment(component.id)[1])
        for symbol in compartments_all:
            if not session.query(Compartment).filter(Compartment.bigg_id == symbol).count():
                compartmentObject = Compartment(bigg_id = symbol, name = '')
                session.add(compartmentObject)
"""  
                if len(component.id.split('_'))>1:


                    if not session.query(Compartment).filter(Compartment.name == parse.split_compartment(component.id)[1]).count():
                        compartmentObject = Compartment(name = parse.split_compartment(component.id)[1])
                        session.add(compartmentObject)

                else:
                    if not session.query(Compartment).filter(Compartment.name == 'none').count():
                        compartmentObject = Compartment(name = 'none')
                        session.add(compartmentObject)
"""           

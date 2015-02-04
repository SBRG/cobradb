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
                logging.exception("%s. In model %s" % (e, model.id))
                continue

            linkouts = [('KEGGID', 'kegg_id'), 
                        ('CASNUMBER', 'cas_number'), 
                        ('SEED', 'seed'),
                        ('METACYC', 'metacyc'),
                        ('CHEBI', 'chebi'),
                        ('BRENDA', 'brenda'),
                        ('UPA', 'upa'),
                        ('HMDB', 'hmdb'),
                        ('BIOPATH', 'biopath'),
                        ('REACTOME', 'reactome'),
                        ('LIPIDMAPS', 'lipidmaps'),
                        ('CASID', 'casid')]
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
                                              formula=str(_formula))
                
                session.add(metaboliteObject)
                for _key in component.notes.keys():
                    if _key != 'FORMULA1':
                        linkout = LinkOut(external_id = found[parse_linkout_str(_key)], external_source = _key, type = "metabolite", ome_id = metaboliteObject.id)
                        session.add(linkout)
            else:
                for _key in component.notes.keys():
                    if _key !=  'FORMULA1':
                        external_link = (session.query(LinkOut)
                                        .filter(LinkOut.external_source == _key)
                                        .filter(LinkOut.ome_id == metabolite_db.id)
                                        .first())
                        if external_link != None:
                            if external_link.external_id == None or external_link.external_id == "":
                                external_link.external_id = found[parse_linkout_str(_key)]
                        else:
                            linkout = LinkOut(external_id = found[parse_linkout_str(_key)], external_source = _key, type = "metabolite", ome_id = metaboliteObject.id)
                            session.add(linkout)
                            

                            
def loadReactions(session, model_list):
    for model in model_list:
        for reaction in model.reactions:
            reaction_db = queries.get_reaction(session, reaction.id)
            if reaction_db is None:
                new_object = Reaction(bigg_id=reaction.id, name=reaction.name, # fix this, see BiGG2 issue #29
                                      notes='', reaction_hash=parse.hash_reaction(reaction))
                session.add(new_object)
            else:
                if reaction_db.name == '' or (reaction_db.name == reaction_db.bigg_id and (reaction.name != '' or reaction.name != None)):
                    reaction_db.name = reaction.name


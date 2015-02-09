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
                        'CASID']
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
                             
            found = {}
            for linkout in linkouts:
                try:
                    vl = parse_linkout_str(component.notes.get(linkout))
                except:
                    continue
                found[linkout] = vl
            
            if metabolite_db is None:

                # look for the formula
                if component.notes.get("FORMULA1") is not None:
                    _formula = component.notes.get("FORMULA1")
                else:
                    _formula = component.formula

                metabolite_db = Metabolite(bigg_id=met_id,
                                              name=component.name,
                                              formula=str(_formula))
                session.add(metabolite_db)
                session.commit()
    
            #load new linkouts even ones that are pointing to previously created universal
            #metabolites. The only scenario where we don't load a linkout is if the 
            #external id and metabolite is exactly the same as a previous linkout.
            
            for _key in component.notes.keys():
                if _key !=  'FORMULA1' and _key != 'FORMULA':
                    external_id_string  = found[parse_linkout_str(_key)].strip()
                    if external_id_string.lower() != "none":
                        for external_id in [x.strip() for x in external_id_string.split(',')]:
                            if (session
                                .query(LinkOut)
                                .filter(LinkOut.external_id == external_id)
                                .filter(LinkOut.external_source == _key)
                                .filter(LinkOut.type == "metabolite")
                                .filter(LinkOut.ome_id == metabolite_db.id)
                                .count()) == 0:
                                linkout = LinkOut(external_id = external_id, 
                                                    external_source = _key, 
                                                    type = "metabolite", 
                                                    ome_id = metabolite_db.id)
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


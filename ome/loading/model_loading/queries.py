# -*- coding: utf-8 -*-

# A collection of common queries from the model loading scripts

from ome.base import Chromosome, Synonym
from ome.models import *

def get_model(session, bigg_id):
    return (session
            .query(Model)
            .filter(Model.bigg_id == bigg_id)
            .first())
    
def get_reaction(session, bigg_id):
    return (session
            .query(Reaction)
            .filter(Reaction.name == bigg_id)
            .first())

def chromosomes_for_genome(session, genome_id):
    return (session
            .query(Chromosome)
            .filter(Chromosome.genome_id == genome_id)
            .all())

def has_model_gene(session, model_id, gene_id):
    return (session
            .query(ModelGene)
            .join(Gene)
            .filter(ModelGene.model_id == model_id)
            .filter(ModelGene.gene_id == gene_id)
            .count() > 0)

def add_model_gene(session, model_id, gene_id):
    new_object = ModelGene(model_id=model_id, gene_id=gene_id)
    session.add(new_object)
    # session.commit()            # TODO commit here? or later?
                

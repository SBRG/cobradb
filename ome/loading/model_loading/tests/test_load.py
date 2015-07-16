# -*- coding: utf-8 -*-

from ome.loading.model_loading import load_model
from ome.loading.component_loading import load_genome
from ome.loading import AlreadyLoadedError
from ome import base
from ome.models import *
from ome.components import *
from ome import settings

from sqlalchemy.orm import aliased
from sqlalchemy import func
import pytest
import os
from os.path import join

def test_load_model(test_genbank, test_model, test_db, test_prefs, setup_logger):
    session = base.Session()

    # preferences. TODO these would be better as arguments to load_model
    settings.reaction_id_prefs = test_prefs['reaction_id_prefs']
    settings.reaction_hash_prefs = test_prefs['reaction_hash_prefs']
    settings.gene_reaction_rule_prefs = test_prefs['gene_reaction_rule_prefs']

    timestamp = '2014-9-16 14:26:22'
    pmid = 'pmid:25575024'
    # can't load the model without the genome
    with pytest.raises(Exception):
        for x in range(2):
            load_model(test_model[x], test_genbank[x]['genome_id'], timestamp,
                       pmid, session, dump_directory=None,
                       published_directory=None, polished_directory=None)
    
    for x in range(2):
        # load the test genomes
        load_genome(test_genbank[x]['path'], session)
        # load the models
        load_model(test_model[x]['path'], test_genbank[x]['genome_id'],
                   timestamp, pmid, session, dump_directory=None,
                   published_directory=None, polished_directory=None)
        
    # load the third model
    load_model(test_model[2]['path'], test_genbank[0]['genome_id'], timestamp,
               pmid, session, dump_directory=None, published_directory=None,
               polished_directory=None)
    
    # test the model
    assert session.query(Model).count() == 3
    assert session.query(Genome).count() == 2
    assert session.query(Chromosome).count() == 2
    assert session.query(Reaction).count() == 98
    assert session.query(ModelReaction).count() == 286
    assert session.query(CompartmentalizedComponent).count() == 73
    assert session.query(ModelCompartmentalizedComponent).count() == 72 * 3 + 1
    assert session.query(Metabolite).count() == 55
    assert session.query(Gene).count() == 281 # b4151 and b4152 are in genome1 but not in model1
    assert session.query(ModelGene).count() == 414
    
    # test linkouts
    result = (session
              .query(LinkOut.external_source, LinkOut.external_id, LinkOut.ome_id)
              .join(Metabolite, Metabolite.id == LinkOut.ome_id)
              .filter(Metabolite.bigg_id == '13dpg')
              .filter(LinkOut.external_source == 'KEGGID')
              .all())
    assert len(result) == 2

    # TODO test charge

    # check s0001
    assert session.query(Gene).filter(Gene.bigg_id == 's0001').count() == 3
    
    # check alternate transcripts
    # these are in model 2 but not in model 1:
    assert (session
            .query(ModelGene)
            .join(Gene)
            .join(Model)
            .filter(Gene.bigg_id.in_(['b4151', 'b4152']))
            .filter(Model.bigg_id == 'Ecoli_core_model')
            .count()) == 0
    assert (session
            .query(ModelGene)
            .join(Gene)
            .join(Model)
            .filter(Gene.bigg_id.in_(['b4151', 'b4152']))
            .filter(Model.bigg_id == 'Ecoli_core_model_2')
            .count()) == 2

    # these alt. transcripts in model 1:
    GeneSource = aliased(Gene)
    assert (session
            .query(Gene)
            .join(GeneSource,
                  GeneSource.id == Gene.alternative_transcript_of)
            .filter(Gene.bigg_id == '904_AT1')
            .filter(GeneSource.bigg_id == 'b0904')
            .count()) == 1
    assert session.query(Synonym).filter(Synonym.synonym == '904').count() == 3 # 3 in first model, 0 in second model
    assert session.query(Gene).filter(Gene.name == 'focA').count() == 3 # 3 in first model, 0 in second model
    
    assert session.query(Gene).filter(Gene.bigg_id == '904_AT1').count() == 1
    assert session.query(Gene).filter(Gene.bigg_id == '904_AT12').count() == 1
    assert session.query(Gene).filter(Gene.bigg_id == 'b0904').count() == 2
    assert session.query(Gene).filter(Gene.bigg_id == 'frdB').count() == 0
    assert session.query(Gene).filter(Gene.bigg_id == 'frdB_AT1').count() == 1
    assert session.query(ModelGene).join(Gene).filter(Gene.bigg_id == '904_AT1').count() == 1
    assert session.query(ModelGene).join(Gene).filter(Gene.bigg_id == '904_AT12').count() == 1
    assert session.query(ModelGene).join(Gene).filter(Gene.bigg_id == 'gene_with_period_AT22').count() == 1
    
    assert session.query(Synonym).filter(Synonym.ome_id == session.query(Gene).filter(Gene.bigg_id == '904_AT1').first().id).count() == 8
    assert session.query(Synonym).filter(Synonym.ome_id == session.query(Gene).filter(Gene.bigg_id == '904_AT12').first().id).count() == 8

    # make sure the locus tag b4153 is back in this gene_reaction_rule (in place of frdB)
    assert (session
            .query(ModelReaction)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.bigg_id == 'FRD7')
            .first()).gene_reaction_rule == '(904_AT12 and gene_with_period_AT22 and b4153 and b4154)'

    # (2 ny-n) Model 1 has a different ACALD from models 2 and 5. The
    # reaction-hash-prefs file should force the second and third models to have
    # ACALD, and increment the first model.
    assert (session
            .query(ModelReaction)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .join(Model)
            .filter(Reaction.bigg_id == 'ACALD_1')
            .filter(Reaction.pseudoreaction == False)
            .filter(Model.bigg_id == 'Ecoli_core_model')
            .count() == 1)
    assert (session
            .query(ModelReaction)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .join(Model)
            .filter(Reaction.bigg_id == 'ACALD')
            .filter(Model.bigg_id == 'Ecoli_core_model_2')
            .count() == 1)
    assert (session
            .query(ModelReaction)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.bigg_id == 'ACALD_2')
            .count() == 0)
    # (3a ny-n) Matches existing reaction. PFL was renamed in model 2.
    assert (session
            .query(ModelReaction)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.bigg_id == 'PFL')
            .count() == 3)
    # (3a yynn) bigg id and hash both return matches. PDH in model 2 matches ENO from
    # model 1.
    assert (session
            .query(ModelReaction)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.bigg_id == 'ENO')
            .count() == 3)
    assert (session
            .query(ModelReaction)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.bigg_id == 'PDH')
            .count() == 3)

    # pseudoreactions. ATPM should be prefered to ATPM_NGAM based on
    # reaction-id-prefs file. Thus, ATPM should be present 3 times, once with
    # the ATPM(NGAM) synonym, and never as ATPM_1.
    assert (session
            .query(ModelReaction)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.bigg_id == 'ATPM')
            .filter(Reaction.pseudoreaction == True)
            .count() == 3)
    assert (session
            .query(OldIDSynonym)
            .join(Synonym, OldIDSynonym.synonym_id == Synonym.id)
            .filter(Synonym.synonym == 'ATPM(NGAM)')
            .join(ModelReaction, ModelReaction.id == OldIDSynonym.ome_id)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.bigg_id == 'ATPM')
            .count() == 1) 
    assert (session
            .query(ModelReaction)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.bigg_id == 'ATPM_1')
            .filter(Reaction.pseudoreaction == True)
            .count() == 0)
    # NTP1 should not be a pseudoreaction, and should not get incremented
    assert (session
            .query(ModelReaction)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.bigg_id == 'NTP1')
            .filter(Reaction.pseudoreaction == False)
            .count() == 1)
    assert (session
            .query(ModelReaction)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.bigg_id == 'NTP1_1')
            .filter(Reaction.pseudoreaction == False)
            .count() == 0)
    assert (session
            .query(OldIDSynonym)
            .join(Synonym, OldIDSynonym.synonym_id == Synonym.id)
            .filter(Synonym.synonym == 'NTP1')
            .join(ModelReaction, ModelReaction.id == OldIDSynonym.ome_id)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.bigg_id == 'NTP1')
            .count() == 1) 

    # gene reaction matrix
    mr_db = (session
             .query(GeneReactionMatrix, ModelReaction)
             .join(ModelReaction, ModelReaction.id == GeneReactionMatrix.model_reaction_id)
             .join(Model, Model.id == ModelReaction.model_id)
             .join(Reaction, Reaction.id == ModelReaction.reaction_id)
             .filter(Model.bigg_id == 'Ecoli_core_model')
             .filter(Reaction.bigg_id == 'ATPM')
             .filter(ModelReaction.gene_reaction_rule == '')
             .all())
    assert len(mr_db) == 0
    mr_db = (session
             .query(GeneReactionMatrix, ModelReaction)
             .join(ModelReaction, ModelReaction.id == GeneReactionMatrix.model_reaction_id)
             .join(Model, Model.id == ModelReaction.model_id)
             .join(Reaction, Reaction.id == ModelReaction.reaction_id)
             .filter(Model.bigg_id == 'Ecoli_core_model')
             .filter(Reaction.bigg_id == 'NTP1')
             .filter(ModelReaction.gene_reaction_rule == '(b0650) or (b4161)')
             .all())
    assert len(mr_db) == 2

    # gene_reaction_rule_prefs to fix reaction rules
    mr_db = (session
             .query(ModelReaction)
             .join(Reaction, Reaction.id == ModelReaction.reaction_id)
             .filter(Model.bigg_id == 'Ecoli_core_model')
             .filter(Reaction.bigg_id == 'ACKr')
             .first())
    assert mr_db.gene_reaction_rule == '(b1849 or b2296 or b3115)'

    # old ids
    assert (session
            .query(OldIDSynonym)
            .join(Synonym, OldIDSynonym.synonym_id == Synonym.id)
            .filter(Synonym.synonym == 'ACALD')
            .join(ModelReaction, ModelReaction.id == OldIDSynonym.ome_id)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.bigg_id == 'ACALD_1')
            .count() == 1)

    assert (session
            .query(OldIDSynonym)
            .join(Synonym, OldIDSynonym.synonym_id == Synonym.id)
            .filter(Synonym.synonym == 'gln_L_c')
            .join(ModelCompartmentalizedComponent,
                  ModelCompartmentalizedComponent.id == OldIDSynonym.ome_id)
            .join(CompartmentalizedComponent,
                  CompartmentalizedComponent.id == ModelCompartmentalizedComponent.compartmentalized_component_id)
            .join(Metabolite,
                  Metabolite.id == CompartmentalizedComponent.component_id)
            .filter(Metabolite.bigg_id == 'gln__L')
            .count() == 3)

    # remove leading underscores (_13dpg in Model 1)
    assert (session
            .query(Metabolite)
            .filter(Metabolite.bigg_id == '_13dpg')
            .first()) is None

    # linkouts
    assert (session
            .query(LinkOut)
            .join(Metabolite, Metabolite.id == LinkOut.ome_id)
            .filter(Metabolite.bigg_id == '13dpg')
            .filter(LinkOut.external_source == 'KEGGID')
            .count()) == 2
    assert (session
            .query(LinkOut)
            .join(Metabolite, Metabolite.id == LinkOut.ome_id)
            .filter(Metabolite.bigg_id == '13dpg')
            .filter(LinkOut.external_source == 'BIOPATH')
            .count()) == 1
    assert (session
            .query(LinkOut)
            .join(Metabolite, Metabolite.id == LinkOut.ome_id)
            .filter(Metabolite.bigg_id == '13dpg')
            .filter(LinkOut.external_source == 'CHEBI')
            .count()) == 5

    # formulas added in second model
    assert (session
            .query(Metabolite)
            .filter(Metabolite.bigg_id == 'atp')
            .first()).formula == 'C10H12N5O13P3'
    
    # test reaction attributes
    r_db =  (session.query(ModelReaction)
             .join(Reaction)
             .filter(Reaction.bigg_id == 'GAPD')
             .first())
    assert r_db.objective_coefficient == 0
    assert r_db.upper_bound == 1000
    assert r_db.lower_bound == -1000

    # can't load the same model twice
    with pytest.raises(AlreadyLoadedError):
        load_model(test_model[0]['path'], test_genbank[0]['genome_id'],
                   timestamp, pmid, session, dump_directory=None,
                   published_directory=None, polished_directory=None)

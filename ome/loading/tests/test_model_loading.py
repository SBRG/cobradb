# -*- coding: utf-8 -*-

from ome.loading.model_loading import load_model, GenbankNotFound
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


@pytest.mark.usefixtures('load_models')
class TestsWithModels:
    def test_cannot_load_same_model_twice(self, session, test_model_files):
        model_details = test_model_files[0]
        with pytest.raises(AlreadyLoadedError):
            load_model(model_details['path'], model_details['pmid'],
                       model_details['genome_ref'], session)

    def test_counts(self, session):
        # test the model
        assert session.query(Model).count() == 3
        assert session.query(Genome).count() == 2
        assert session.query(Chromosome).count() == 2
        assert session.query(Reaction).count() == 98
        assert session.query(ModelReaction).count() == 288
        assert session.query(CompartmentalizedComponent).count() == 73
        assert session.query(ModelCompartmentalizedComponent).count() == 72 * 3 + 1
        assert session.query(Metabolite).count() == 55
        assert session.query(Gene).count() == 286
        assert session.query(ModelGene).count() == 415

    def test_no_charge_in_linkouts(self, session):
        assert (session
                .query(Synonym)
                .join(DataSource)
                .filter(DataSource.name == 'CHARGE')
                .count()) == 0

    def test_s0001(self, session):
        assert session.query(Gene).filter(Gene.bigg_id == 's0001').count() == 3

    def test_name_scrubbing(self, session):
        assert session.query(Reaction).filter(Reaction.bigg_id == 'ACALD').first().name == 'Acetaldehyde dehydrogenase (acetylating)'

    def test_name_filtering_met(self, session):
        # Filter for higher quality descriptive names
        assert session.query(Metabolite).filter(Metabolite.name == 'E4P c').first() is None
        assert session.query(Metabolite).filter(Metabolite.name == 'D-Erythrose-4-phosphate').first() is not None

    def test_name_filtering_rxn(self, session):
        # Filter for higher quality descriptive names
        assert session.query(Reaction).filter(Reaction.name == 'Atps4   z').first() is None
        assert session.query(Reaction).filter(Reaction.name == 'ATP synthase (four protons for one ATP)').first() is not None

    # alternative transcripts
    def test_alternative_transcripts_counts(self, session):
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

    def test_alternative_transcripts(self, session):
        # these alt. transcripts in model 1:
        GeneSource = aliased(Gene)
        res_db = (session
                  .query(Gene)
                  .join(GeneSource,
                        GeneSource.id == Gene.alternative_transcript_of)
                  .filter(Gene.bigg_id == '904_AT1')
                  .filter(GeneSource.bigg_id == 'b0904'))
        assert res_db.count() == 1
        assert res_db.first().name == 'focA'

    def test_gene_reaction_rule_handling(self, session):
        # make sure the locus tag b4153 is back in this gene_reaction_rule (in
        # place of frdB).
        # NOTE: COBRApy now reformats these without extra parens.
        res_db = (session
                  .query(ModelReaction.gene_reaction_rule)
                  .join(Reaction, Reaction.id == ModelReaction.reaction_id)
                  .filter(Reaction.bigg_id == 'FRD7')
                  .all())
        assert '904_AT12 and gene_with_period_AT22 and b4153 and b4154' in [x[0] for x in res_db]

    def tests_reaction_collisions(self, session):
        # (2 ny-n) Model 1 has a different ACALD from models 2 and 3. The
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

    def test_repeated_reaction_different_bounds(self, session):
        # exact copies should get merged
        assert (session
                .query(Reaction)
                .join(ModelReaction)
                .join(Model)
                .filter(Model.bigg_id == 'Ecoli_core_model')
                .filter(Reaction.bigg_id == 'EX_glu__L_e')
                .count()) == 1
        # copies with different bounds should be separated
        res_db = (session
                  .query(Reaction, ModelReaction)
                  .join(ModelReaction)
                  .join(Model)
                  .filter(Model.bigg_id == 'Ecoli_core_model')
                  .filter(Reaction.bigg_id == 'EX_gln__L_e'))
        assert res_db.count() == 2
        assert {(x.lower_bound, x.upper_bound) for x in (y[1] for y in res_db)} == {(0, 50), (0, 1000)}

    def tests_pseudoreactions(self, session):
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
                .filter(OldIDSynonym.type == 'model_reaction')
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
                .filter(OldIDSynonym.type == 'model_reaction')
                .join(Synonym, OldIDSynonym.synonym_id == Synonym.id)
                .filter(Synonym.synonym == 'NTP1')
                .join(ModelReaction, ModelReaction.id == OldIDSynonym.ome_id)
                .join(Reaction, Reaction.id == ModelReaction.reaction_id)
                .filter(Reaction.bigg_id == 'NTP1')
                .count() == 1)

    def test_gene_reaction_matrix(self, session):
        len_mr_db = (session
                     .query(GeneReactionMatrix, ModelReaction)
                     .join(ModelReaction, ModelReaction.id == GeneReactionMatrix.model_reaction_id)
                     .join(Model, Model.id == ModelReaction.model_id)
                     .join(Reaction, Reaction.id == ModelReaction.reaction_id)
                     .filter(Model.bigg_id == 'Ecoli_core_model')
                     .filter(Reaction.bigg_id == 'ATPM')
                     .filter(ModelReaction.gene_reaction_rule == '')
                     .count())
        assert len_mr_db == 0
        len_mr_db = (session
                     .query(GeneReactionMatrix, ModelReaction)
                     .join(ModelReaction, ModelReaction.id == GeneReactionMatrix.model_reaction_id)
                     .join(Model, Model.id == ModelReaction.model_id)
                     .join(Reaction, Reaction.id == ModelReaction.reaction_id)
                     .filter(Model.bigg_id == 'Ecoli_core_model')
                     .filter(Reaction.bigg_id == 'NTP1')
                     .filter(ModelReaction.gene_reaction_rule == 'b0650 or b4161')
                     .count())
        assert len_mr_db == 2

    def test_gene_reaction_matrix_multiple_reaction_copies(self, session):
        res_db = (session
                  .query(GeneReactionMatrix, Gene, ModelReaction)
                  .join(ModelGene)
                  .join(Gene)
                  .join(ModelReaction)
                  .join(Reaction)
                  .join(Model, Model.id == ModelReaction.model_id)
                  .filter(Reaction.bigg_id == 'ADK1')
                  .filter(Model.bigg_id == 'Ecoli_core_model')
                  .all())
        # two genes
        assert len(res_db) == 2
        assert {x[1].bigg_id for x in res_db} == {'b0474', 'b0474_test_copy'}
        # should be separate entries for each ModelReaction
        assert res_db[0][2].id != res_db[1][2].id

    def test_gene_reaction_rule_prefs(self, session):
        # gene_reaction_rule_prefs to fix reaction rules
        mr_db = (session
                .query(ModelReaction)
                .join(Reaction, Reaction.id == ModelReaction.reaction_id)
                .filter(Model.bigg_id == 'Ecoli_core_model')
                .filter(Reaction.bigg_id == 'ACKr')
                .first())
        assert mr_db.gene_reaction_rule == 'b1849 or b2296 or b3115'

    def test_multiple_reaction_copies(self, session):
        # make sure both copies of ADK1 are here
        res = (session
               .query(ModelReaction)
               .join(Reaction)
               .join(Model)
               .filter(Reaction.bigg_id == 'ADK1')
               .filter(Model.bigg_id == 'Ecoli_core_model')
               .order_by(ModelReaction.copy_number)
               .all())
        assert len(res) == 2
        assert res[0].copy_number == 1
        assert res[0].lower_bound == -1000
        assert res[1].copy_number == 2
        assert res[1].lower_bound == 0

    def test_multiple_metabolite_copies(self, session):
        # e.g. ID collision
        res = (session
               .query(Synonym.synonym, Component, Compartment)
               .join(OldIDSynonym)
                .filter(OldIDSynonym.type == 'model_compartmentalized_component')
               .filter(Synonym.type == 'compartmentalized_component')
               .join(ModelCompartmentalizedComponent, ModelCompartmentalizedComponent.id == OldIDSynonym.ome_id)
               .join(CompartmentalizedComponent)
               .join(Component)
               .join(Compartment)
               .join(Model)
               .filter(Component.bigg_id == 'glc__D')
               .filter(Compartment.bigg_id == 'e')
               .filter(Model.bigg_id == 'Ecoli_core_model'))
        assert res.count() == 2

    def test_multiple_metabolite_copies_2(self, session):
        # e.g. ID collision
        res = (session
               .query(Synonym.synonym, Component, Compartment)
               .join(OldIDSynonym)
                .filter(OldIDSynonym.type == 'model_compartmentalized_component')
               .filter(Synonym.type == 'component')
               .join(ModelCompartmentalizedComponent, ModelCompartmentalizedComponent.id == OldIDSynonym.ome_id)
               .join(CompartmentalizedComponent)
               .join(Component)
               .join(Compartment)
               .join(Model)
               .filter(Component.bigg_id == 'glc__D')
               .filter(Compartment.bigg_id == 'e')
               .filter(Model.bigg_id == 'Ecoli_core_model'))
        assert res.count() == 1

    def test_multiple_gene_copies(self, session):
        # for T. maritima, the genes TM0846 and TM_0846 were the same, so
        # multiple OldIdSynonyms
        res_db = (session
                  .query(Synonym.synonym)
                  .join(OldIDSynonym)
                  .filter(OldIDSynonym.type == 'model_gene')
                  .join(ModelGene, ModelGene.id == OldIDSynonym.ome_id)
                  .join(Gene)
                  .join(Model)
                  .filter(Gene.bigg_id == 'b3528')
                  .filter(Model.bigg_id == 'Ecoli_core_model'))
        assert [x[0] for x in res_db] == ['b3528', 'b_3528']

    def test_old_reaction_id(self, session):
        assert (session
                .query(OldIDSynonym)
                .filter(OldIDSynonym.type == 'model_reaction')
                .join(Synonym, OldIDSynonym.synonym_id == Synonym.id)
                .filter(Synonym.synonym == 'ACALD')
                .join(ModelReaction, ModelReaction.id == OldIDSynonym.ome_id)
                .join(Reaction, Reaction.id == ModelReaction.reaction_id)
                .filter(Reaction.bigg_id == 'ACALD_1')
                .count()) == 1

    def test_old_metabolite_id(self, session):
        assert (session
                .query(OldIDSynonym)
                .filter(OldIDSynonym.type == 'model_compartmentalized_component')
                .join(Synonym, OldIDSynonym.synonym_id == Synonym.id)
                .filter(Synonym.synonym == 'gln_L_c')
                .join(ModelCompartmentalizedComponent,
                    ModelCompartmentalizedComponent.id == OldIDSynonym.ome_id)
                .join(CompartmentalizedComponent,
                    CompartmentalizedComponent.id == ModelCompartmentalizedComponent.compartmentalized_component_id)
                .join(Metabolite,
                    Metabolite.id == CompartmentalizedComponent.component_id)
                .filter(Metabolite.bigg_id == 'gln__L')
                .count()) == 3

    def test_old_gene_id(self, session):
        assert (session
                .query(OldIDSynonym)
                .filter(OldIDSynonym.type == 'model_gene')
                .join(Synonym, OldIDSynonym.synonym_id == Synonym.id)
                .filter(Synonym.synonym == 'gene_with_period.22')
                .join(ModelGene, ModelGene.id == OldIDSynonym.ome_id)
                .join(Gene, Gene.id == ModelGene.gene_id)
                .filter(Gene.bigg_id == 'gene_with_period_AT22')
                .count()) == 1

    def test_leading_underscores(self, session):
        # remove leading underscores (_13dpg in Model 1)
        assert (session
                .query(Metabolite)
                .filter(Metabolite.bigg_id == '_13dpg')
                .first()) is None

    def test_linkout_old_bigg_id(self, session):
        res_db = (session
                  .query(Metabolite, Synonym, DataSource)
                  .filter(Metabolite.bigg_id == '13dpg')
                  .join(CompartmentalizedComponent)
                  .join(Synonym, Synonym.ome_id == CompartmentalizedComponent.id)
                  .filter(Synonym.type == 'compartmentalized_component')
                  .join(DataSource)
                  .filter(DataSource.bigg_id == 'old_bigg_id')
                  .all())
        assert len(res_db) == 2
        assert set(x[1].synonym for x in res_db) == {'_13dpg_c', '13dpg_c'}

    def tests_reaction_attributes(self, session):
        r_db =  (session.query(ModelReaction)
                .join(Reaction)
                .filter(Reaction.bigg_id == 'GAPD')
                .first())
        assert r_db.objective_coefficient == 0
        assert r_db.upper_bound == 1000
        assert r_db.lower_bound == -1000
        assert r_db.subsystem == 'Glycolysis/Gluconeogenesis'

    def test_metabolite_attributes(self, session):
        res_db = (session
                  .query(ModelCompartmentalizedComponent)
                  .join(CompartmentalizedComponent)
                  .join(Component)
                  .join(Model)
                  .filter(Component.bigg_id == '13dpg')
                  .filter(Model.bigg_id == 'Ecoli_core_model')
                  .first())
        assert res_db.formula == 'C3H4O10P2'
        assert res_db.charge == -4

    def test_reaction_direction_hash_1(self, session):
        # Use PGI direction from the first model
        rm_db =  (session
                  .query(ReactionMatrix)
                  .join(Reaction)
                  .join(CompartmentalizedComponent)
                  .join(Component)
                  .filter(Reaction.bigg_id == 'PGI')
                  .filter(Component.bigg_id == 'g6p')
                  .first())
        assert rm_db.stoichiometry == 1

    def test_reaction_direction_hash_2(self, session):
        # No incremented PGI, even though PGI in the two models are in different directions
        assert session.query(Reaction).filter(Reaction.bigg_id == 'PGI_1').count() == 0

    def test_reaction_direction_hash_3(self, session):
        # The PGI in model 2 should have its upper and lower bounds reversed
        r_db = (session
                .query(ModelReaction)
                .join(Model)
                .join(Reaction)
                .filter(Model.bigg_id == 'Ecoli_core_model_2')
                .filter(Reaction.bigg_id == 'PGI')
                .first())
        assert r_db.lower_bound == -1000
        assert r_db.upper_bound == 0

    def test_reaction_direction_hash_4(self, session):
        # ignore reverse hash mapping when reactions are in the same model
        assert session.query(Reaction).filter(Reaction.bigg_id == 'FRD7').count() == 1
        assert session.query(Reaction).filter(Reaction.bigg_id == 'SUCDi').count() == 1

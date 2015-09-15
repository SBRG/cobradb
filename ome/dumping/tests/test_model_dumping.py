# -*- coding: utf-8 -*-

from ome.dumping.model_dumping import dump_model
from ome.loading import load_model, load_genome
from ome.base import *
from ome.models import *

import pytest
import logging
import cobra.io
from os.path import join, exists
import shutil


# Dumping
@pytest.fixture(scope='session')
def dumped_model(load_models, session):
    bigg_id = session.query(Model).first().bigg_id
    model = dump_model(bigg_id)
    return model


# TIP: always use the fixtures to avoid trouble
def test_cannot_dump_unknown_model(dumped_model, session):
    with pytest.raises(Exception):
        dump_model('C3PO', session)


# Model content
def test_dumped_model(dumped_model):
    # COBRApy uses the description as the ID sometimes. See https://github.com/opencobra/cobrapy/pull/152
    assert dumped_model.description == dumped_model.id

    assert len(dumped_model.reactions) == 96
    assert len(dumped_model.metabolites) == 73
    assert len(dumped_model.genes) == 140
    assert dumped_model.genes.get_by_id('b0114').name == 'aceE'
    assert dumped_model.genes.get_by_id('b3528').name == 'dctA'

    # check reaction
    assert 'GAPD' in dumped_model.reactions
    assert dumped_model.reactions.get_by_id('GAPD').name == 'glyceraldehyde-3-phosphate dehydrogenase'

    # check metabolite
    assert 'g3p_c' in dumped_model.metabolites
    assert dumped_model.metabolites.get_by_id('g3p_c').name == 'Glyceraldehyde-3-phosphate'


def test_solve_model(dumped_model):
    # test solve
    dumped_model.reactions.get_by_id('EX_glc__D_e').lower_bound = -10
    sol = dumped_model.optimize()
    assert sol.f > 0.5 and sol.f < 1.0


def test_sbml_dump(dumped_model, tmpdir):
    # temp sbml dump
    m_file = join(str(tmpdir), 'test_model.sbml')
    cobra.io.write_sbml_model(dumped_model, m_file)
    loaded_model = cobra.io.read_sbml_model(m_file)
    # check dumped_model id
    assert loaded_model.id == dumped_model.id
    assert loaded_model.description == dumped_model.description

    # Check for buggy metabolite. This was being saved without an ID in the SBML
    # dumped_model for some reason.
    assert 'formmfr_b_c' in [x.id for x in dumped_model.metabolites]
    assert 'formmfr_b_c' in [x.id for x in loaded_model.metabolites]
    assert '' not in [x.id for x in loaded_model.metabolites]

    # delete test directory
    shutil.rmtree(str(tmpdir))


def test_dumped_pseudoreactions(dumped_model):
    # make sure ATPM and NTP1 are represented
    r1 = dumped_model.reactions.get_by_id('NTP1')
    r2 = dumped_model.reactions.get_by_id('ATPM')
    assert r1.lower_bound == 3.15
    assert r2.lower_bound == 8.39
    assert r1.notes['original_bigg_id'] == 'NTP1'
    assert r2.notes['original_bigg_id'] == 'ATPM(NGAM)'


# Old IDs
def test_reaction_notes(dumped_model):
    assert dumped_model.reactions.get_by_id('ATPM').notes['original_bigg_id'] == 'ATPM(NGAM)'


def test_metabolite_notes(dumped_model):
    assert dumped_model.metabolites.get_by_id('13dpg_c').notes['original_bigg_id'] == '_13dpg_c'


def test_gene_notes(dumped_model):
    assert dumped_model.genes.get_by_id('b0114').notes['original_bigg_id'] == 'b0114'

# -*- coding: utf-8 -*-

from cobradb import settings
from cobradb import base
from cobradb.loading.component_loading import load_genome
from cobradb.loading.model_loading import load_model

import pytest
from sqlalchemy import create_engine
import sys
import os
from os.path import join, realpath, dirname
import cobra.io
import logging


test_data_dir = realpath(join(dirname(__file__), 'test_data'))


@pytest.fixture(scope='session')
def session(request):
    """Make a session"""
    def teardown():
        base.Session.close_all()
    request.addfinalizer(teardown)

    return base.Session()


@pytest.fixture(scope='session')
def setup_logger():
    logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)


@pytest.fixture(scope='session')
def test_genbank_files():
    return [(('primary_type', 'core'), join(test_data_dir, 'core.gb')),
            (('primary_type', 'core_2'), join(test_data_dir, 'core_2.gb'))]


@pytest.fixture(scope='session')
def test_model_files():
    return [{'path': join(test_data_dir, 'ecoli_core_model.xml'),
             'genome_ref': ('primary_type', 'core'),
             'pmid': ('pmid', '25575024')},
            {'path': join(test_data_dir, 'ecoli_core_model_2.xml'),
             'genome_ref': ('primary_type', 'core_2'),
             'pmid': ('pmid', '25575024')},
            {'path': join(test_data_dir, 'ecoli_core_model_3.xml'),
             'genome_ref': ('primary_type', 'core'),
             'pmid': ('pmid', '25575024')}]


@pytest.fixture(scope='session')
def test_prefs():
    return {'reaction_id_prefs': join(test_data_dir, 'reaction-id-prefs.txt'),
            'reaction_hash_prefs': join(test_data_dir, 'reaction-hash-prefs.txt'),
            'data_source_preferences': join(test_data_dir, 'data-source-prefs.txt'),
            'gene_reaction_rule_prefs': join(test_data_dir, 'gene-reaction-rule-prefs.txt')}


@pytest.fixture(scope='session')
def test_db_create(setup_logger):
    user = settings.postgres_user
    test_db = settings.postgres_test_database
    # make sure the test database is clean
    os.system('dropdb %s' % test_db)
    os.system('createdb %s -U %s' % (test_db, user))
    logging.info('Dropped and created database %s' % test_db)


@pytest.fixture(scope='session')
def test_db(request, test_db_create):
    user = settings.postgres_user
    test_db = settings.postgres_test_database
    engine = create_engine("postgresql://%s:%s@%s/%s" % (user,
                                                         settings.postgres_password,
                                                         settings.postgres_host,
                                                         test_db))
    base.Base.metadata.create_all(engine)
    base.Session.configure(bind=engine)
    logging.info('Loaded database schema')

    def teardown():
        # close all sessions. Comment this line out to see if cobradb functions are
        # closing their sessions properly.
        base.Session.close_all()
        # clear the db for the next test
        base.Base.metadata.drop_all(engine)
        logging.info('Dropped database schema')
    request.addfinalizer(teardown)


@pytest.fixture(scope='session')
def load_genomes(test_db, test_genbank_files, test_prefs, session):
    # preferences
    settings.reaction_id_prefs = test_prefs['reaction_id_prefs']
    settings.reaction_hash_prefs = test_prefs['reaction_hash_prefs']
    settings.gene_reaction_rule_prefs = test_prefs['gene_reaction_rule_prefs']
    settings.data_source_preferences = test_prefs['data_source_preferences']

    # load the test genomes
    for genome_ref, gb in test_genbank_files:
        load_genome(genome_ref, [gb], session)


@pytest.fixture(scope='session')
def load_models(load_genomes, test_model_files, session):
    # fixture load_genomes will have loaded 2 genomes

    out = []
    for model_details in test_model_files:
        # load the models
        out.append(load_model(model_details['path'],
                              model_details['pmid'],
                              model_details['genome_ref'],
                              session))
    assert out == ['Ecoli_core_model', 'Ecoli_core_model_2', 'Ecoli_core_model_3']


# Session
#  - test_db_create
#     -> dropbd
#     -> createdb
#
#  - test_db
#     = create_all
#     -> test_function_1
#     = drop_all
#
#  - test_db
#     = create_all
#     -> test_function_2
#     = drop_all

# py.test --capture=no # gives you the logs
# py.test --capture=no -k load # only runs tests with 'load' in the name

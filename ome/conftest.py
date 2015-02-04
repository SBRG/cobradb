# -*- coding: utf-8 -*-

from ome import settings
from ome import base

import pytest
from sqlalchemy import create_engine
import sys
import os
from os.path import join, realpath, dirname
import cobra.io
import logging

@pytest.fixture(scope='session')
def setup_logger():
    logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)

@pytest.fixture(scope='session')
def test_genbank():
    return [{ 'genome_id': 'PRJNA57779-core',
             'path': realpath(join(dirname(__file__), 'test_data', 'core.gb')) },
             { 'genome_id': 'PRJNA57779-core-2',
             'path': realpath(join(dirname(__file__), 'test_data', 'core_2.gb'))}]

@pytest.fixture(scope='session')
def test_model():
    return [{ 'id': 'ecoli_core_model',
             'bigg_id': 'Ecoli_core_model', # TODO bigg_id and id should not be different
             'dir': realpath(join(dirname(__file__), 'test_data')) },
             { 'id': 'ecoli_core_model_2',
             'bigg_id': 'Ecoli_core_model_2', # TODO bigg_id and id should not be different
             'dir': realpath(join(dirname(__file__), 'test_data')) }]

@pytest.fixture(scope='session')
def test_db_create(setup_logger):
    user = settings.postgres_user
    test_db = settings.postgres_test_database
    # make sure the test database is clean
    os.system('dropdb %s' % test_db)
    os.system('createdb %s -U %s' % (test_db, user))
    logging.info('Dropped and created database %s' % test_db)
    
@pytest.fixture(scope='function')
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
        # close all sessions. Comment this line out to see if ome functions are
        # closing their sessions properly.
        base.Session.close_all()
        # clear the db for the next test
        base.Base.metadata.drop_all(engine)
        logging.info('Dropped database schema')
    request.addfinalizer(teardown)

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



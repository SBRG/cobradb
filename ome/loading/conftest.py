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

@pytest.fixture
def test_genbank():
    return { 'genome_id': 'PRJNA57779-core',
             'path': realpath(join(dirname(__file__), 'test_data', 'core.gb')) }

@pytest.fixture
def test_model():
    return { 'id': 'ecoli_core_model',
             'dir': realpath(join(dirname(__file__), 'test_data')) }

@pytest.fixture(scope='session')
def test_db(request):
    user = settings.postgres_user
    test_db = settings.postgres_test_database
    os.system('dropdb %s' % test_db)
    os.system('createdb %s -U %s' % (test_db, user))
    engine = create_engine("postgresql://%s:%s@%s/%s" % (user,
                                                         settings.postgres_password,
                                                         settings.postgres_host,
                                                         test_db))
    base.Base.metadata.create_all(engine)
    base.Session.configure(bind=engine)

    def teardown():
        os.system('dropdb %s' % test_db)
    request.addfinalizer(teardown)

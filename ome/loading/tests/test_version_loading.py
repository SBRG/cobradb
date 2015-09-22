# -*- coding: utf-8 -*-

from ome.base import DatabaseVersion
from ome.loading.version_loading import *

def test_load_version_date(test_db, session):
    load_version_date(session)
    res_db_1 = session.query(DatabaseVersion).one()
    t1 = res_db_1.date_time
    load_version_date(session)
    res_db_2 = session.query(DatabaseVersion).one()
    assert t1 != res_db_2.date_time

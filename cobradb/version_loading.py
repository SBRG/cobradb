# -*- coding: utf-8 -*-

from cobradb.models import DatabaseVersion

from datetime import datetime

def load_version_date(session):
    vers_db = (session
               .query(DatabaseVersion)
               .first())
    time = datetime.now()
    if vers_db is None:
        vers_db = DatabaseVersion(time)
        session.add(vers_db)
    else:
        vers_db.date_time = time
    session.commit()

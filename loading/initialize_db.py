from ome.lib import settings
import os

def set_schema():
    """(re)set the schema for the entire database (wipes existing data)"""
    psql = settings.psql_full
    schema_base = settings.ome_directory + "/schemas/ome"
    #  "2>&1" will redirect stderr into stdout so it also gets logged
    os.system("%s %s < %s_base.sql > log1.log 2>&1" % (psql, settings.dev_database, schema_base))
    os.system("%s %s < %s_components.sql > log2.log 2>&1" % (psql, settings.dev_database, schema_base))
    os.system("%s %s < %s_data.sql > log3.log 2>&1" % (psql, settings.dev_database, schema_base))

    
    
if __name__ == "__main__":
    # creating the schema will drop the entire database
    set_schema()

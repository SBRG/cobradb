"""retrive local user settings"""

from ConfigParser import SafeConfigParser, NoOptionError
import os as os
from os.path import join, split, abspath, isfile, expanduser, dirname
from sys import modules

self = modules[__name__]

# define various filepaths

config = SafeConfigParser()

# overwrite defaults settings with settings from the file
filepath = abspath(join(dirname(__file__), '..', 'settings.ini'))
if isfile(filepath):
    config.read(filepath)
else:
    raise Exception('No settings files at path: %s' % filepath)

# save options as variables
self.postgres_user = config.get('DATABASE', 'postgres_user')
self.postgres_password = config.get('DATABASE', 'postgres_password')
self.postgres_database = config.get('DATABASE', 'postgres_database')
self.postgres_host = config.get('DATABASE', 'postgres_host')
self.postgres_port = config.get('DATABASE', 'postgres_port')
self.postgres_test_database = config.get('DATABASE', 'postgres_test_database')

if self.postgres_host == '' and self.postgres_password == '' \
        and self.postgres_user == '':
    self.db_connection_string = 'postgresql:///%s' % self.postgres_database
else:
    self.db_connection_string = 'postgresql://%s:%s@%s/%s' % \
        (self.postgres_user, self.postgres_password,
            self.postgres_host, self.postgres_database)

self.java = config.get('EXECUTABLES', 'java')

# these are required
try:
    self.model_directory = expanduser(config.get('DATA', 'model_directory'))
except NoOptionError:
    raise Exception('model_directory was not supplied in settings.ini')

try:
    self.refseq_directory = expanduser(config.get('DATA', 'refseq_directory'))
except NoOptionError:
    raise Exception('refseq_directory was not supplied in settings.ini')
try:
    self.model_genome = expanduser(config.get('DATA', 'model_genome'))
except NoOptionError:
    raise Exception('model_genome path was not supplied in settings.ini')

# these are optional
for data_pref in ['compartment_names', 'reaction_id_prefs',
                    'reaction_hash_prefs', 'gene_reaction_rule_prefs',
                    'data_source_preferences']:
    try:
        setattr(self, data_pref, expanduser(config.get('DATA', data_pref)))
    except NoOptionError:
        setattr(self, data_pref, None)

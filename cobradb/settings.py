# -*- coding: utf-8 -*-

"""Retrive local user settings"""

from configparser import SafeConfigParser, NoOptionError
import os
from os.path import join, split, abspath, isfile, expanduser, dirname
from sys import modules
import six

self = modules[__name__]

# define various filepaths

config = SafeConfigParser()

# overwrite defaults settings with settings from the file
filepath = abspath(join(dirname(__file__), '..', 'settings.ini'))
if isfile(filepath):
    config.read(filepath)
else:
    raise Exception('No settings files at path: %s' % filepath)

# prefer environment variables for database settings
env_names = {
    'postgres_host': 'COBRADB_POSTGRES_HOST',
    'postgres_port': 'COBRADB_POSTGRES_PORT',
    'postgres_user': 'COBRADB_POSTGRES_USER',
    'postgres_password': 'COBRADB_POSTGRES_PASSWORD',
    'postgres_database': 'COBRADB_POSTGRES_DATABASE',
    'postgres_test_database': 'COBRADB_POSTGRES_TEST_DATABASE',
}
for setting_name, env_name in six.iteritems(env_names):
    if env_name in os.environ:
        print('Setting %s with environment variable %s' % (setting_name,
                                                           env_name))
        setattr(self, setting_name, os.environ[env_name])
    else:
        setattr(self, setting_name, config.get('DATABASE', setting_name))

# set up the database connection string
if self.postgres_host == '' and self.postgres_password == '' \
        and self.postgres_user == '':
    self.db_connection_string = 'postgresql:///%s' % self.postgres_database
else:
    self.db_connection_string = 'postgresql://%s:%s@%s/%s' % \
        (self.postgres_user, self.postgres_password,
            self.postgres_host, self.postgres_database)

# get the java executable (optional, for running Model Polisher)
if config.has_option('EXECUTABLES', 'java'):
    self.java = config.get('EXECUTABLES', 'java')
else:
    print('No Java executable provided.')

if not config.has_section('DATA'):
    raise Exception('DATA section was not found in settings.ini')

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

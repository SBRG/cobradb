"""retrive local user settings"""

from ConfigParser import SafeConfigParser, NoOptionError
import os as os
from os.path import join, split, abspath, isfile, expanduser
from sys import modules

self = modules[__name__]

# define various filepaths
omelib_directory = join(split(abspath(__file__))[0], "")
ome_directory = join(abspath(join(omelib_directory, "..")), "")

def which(program):
    """returns path to an executable if it is found in the path"""
    fpath, fname = split(program)
    if fpath:
        if isfile(program) and os.access(program, os.X_OK):
            return program
    else:
        paths_to_search = os.environ["PATH"].split(os.pathsep)
        paths_to_search.extend((omelib_directory, ome_directory))
        for path in paths_to_search:
            exe_file = join(path, program)
            if isfile(exe_file) and os.access(exe_file, os.X_OK):
                return exe_file
    if os.name == "nt" and not program.endswith(".exe"):
        return which(program + ".exe")
    return None

def _escape_space(program):
    """escape spaces in for windows"""
    if os.name == "nt" and ' ' in program:
        return '"' + program + '"'
    else:
        return program


config = SafeConfigParser()

# set the default settings
config.add_section("DATABASE")
config.set("DATABASE", "postgres_host", "localhost")
config.set("DATABASE", "postgres_port", "5432")
config.set("DATABASE", "postgres_database", "ome_stage_2")
config.set("DATABASE", "postgres_user", "dbuser")
config.set("DATABASE", "postgres_password", "")
config.set("DATABASE", "postgres_test_database", "ome_test")

config.add_section("DATA")

config.add_section("EXECUTABLES")

# overwrite defaults settings with settings from the file
def load_settings_from_file(filepath='settings.ini', in_ome_dir=True):
    """Reload settings from a different settings file.

    Arguments
    ---------

    filepath: The path to the settings file to use.

    in_ome_dir: Whether or not the path given is a relative path from the ome
    directory.

    """
    if in_ome_dir:
        filepath = join(ome_directory, filepath)
    config.read(filepath)

    # attempt to intellegently determine more difficult settings
    if not config.has_option("DATABASE", "user"):
        if "USERNAME" in os.environ:  # windows
            user = os.environ["USERNAME"]
        elif "USER" in os.environ:  # unix
            user = os.environ["USER"]
        config.set("DATABASE", "user", user)

    # executables
    if not config.has_option("EXECUTABLES", "psql"):
        psql = which("psql91")
        if psql is None:
            psql = which("psql")
        config.set("EXECUTABLES", "psql", psql)
    if not config.has_option("EXECUTABLES", "R"):
        R = which("R")
        config.set("EXECUTABLES", "R", R)
    if not config.has_option("EXECUTABLES", "Rscript"):
        Rscript = which("Rscript")
        config.set("EXECUTABLES", "Rscript", Rscript)
    if not config.has_option("EXECUTABLES", "primer3"):
        primer3 = which("primer3_core")
        config.set("EXECUTABLES", "primer3", primer3)
    if not config.has_option("EXECUTABLES", "cufflinks"):
        cufflinks = which("cufflinks")
        config.set("EXECUTABLES", "cufflinks", cufflinks)
    if not config.has_option("EXECUTABLES", "java"):
        java = which("java")
        config.set("EXECUTABLES", "java", java)

    # save options as variables
    self.postgres_user = config.get("DATABASE", "postgres_user")
    self.postgres_password = config.get("DATABASE", "postgres_password")
    if len(self.postgres_password) > 0:
        os.environ["PGPASSWORD"] = self.postgres_password
    self.postgres_database = config.get("DATABASE", "postgres_database")
    self.postgres_host = config.get("DATABASE", "postgres_host")
    self.postgres_port = config.get("DATABASE", "postgres_port")
    self.postgres_test_database = config.get("DATABASE", "postgres_test_database")
    self.psql = _escape_space(config.get("EXECUTABLES", "psql"))
    self.R = _escape_space(config.get("EXECUTABLES", "R"))
    self.Rscript = _escape_space(config.get("EXECUTABLES", "Rscript"))
    self.primer3 = _escape_space(config.get("EXECUTABLES", "primer3"))
    self.cufflinks = config.get("EXECUTABLES", "cufflinks")
    self.java = config.get("EXECUTABLES", "java")

    # make a psql string with the database options included
    self.psql_full = "%s --host=%s --username=%s --port=%s " % \
                     (self.psql, self.postgres_host, self.postgres_user, self.postgres_port)

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

load_settings_from_file()


del SafeConfigParser, modules

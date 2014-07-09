"""retrive local user settings"""

from ConfigParser import SafeConfigParser
import os as __os
from os.path import split as __split, join as __join, abspath as __abspath, \
    isfile as __isfile
from sys import modules

self = modules[__name__]

# define various filepaths
omlib_directory = __join(__split(__abspath(__file__))[0], "")
om_directory = __join(__abspath(__join(omlib_directory, "..")), "")

def which(program):
    """returns path to an executable if it is found in the path"""
    fpath, fname = __split(program)
    if fpath:
        if __isfile(program) and __os.access(program, __os.X_OK):
            return program
    else:
        paths_to_search = __os.environ["PATH"].split(__os.pathsep)
        paths_to_search.extend((omlib_directory, om_directory))
        for path in paths_to_search:
            exe_file = __join(path, program)
            if __isfile(exe_file) and __os.access(exe_file, __os.X_OK):
                return exe_file
    if __os.name == "nt" and not program.endswith(".exe"):
        return which(program + ".exe")
    return None

def _escape_space(program):
    """escape spaces in for windows"""
    if __os.name == "nt" and ' ' in program:
        return '"' + program + '"'
    else:
        return program


config = SafeConfigParser()
# set the default settings
config.add_section("DATABASE")
config.set("DATABASE", "postgres_host", "localhost:5432")
config.set("DATABASE", "postgres_database", "om")
config.set("DATABASE", "password", "")

config.add_section("MISC")
config.set("MISC", "entrez_email", "SET_ENTREZ_EMAIL")
try:
    default_home = __os.environ["HOME"]
except:
    try:
        default_home = __os.environ["USERPROFILE"]
    except:
        default_home = ""

default_data_dir = __join(default_home, "om_data", "")
config.set("MISC", "home_directory", default_home)
config.set("MISC", "data_directory", default_data_dir)
del default_home, default_data_dir
config.add_section("EXECUTABLES")

# overwrite defaults settings with settings from the file
def load_settings_from_file(filepath="settings.ini", in_omlib=True):
    """reload settings from a different settings file

    filepath: The path to the settings file to use

    in_omlib: Whether or not the path given is a relative path from the omlib
        directory"""
    if in_omlib:
        filepath = __join(omlib_directory, filepath)
    config.read(filepath)

    # attempt to intellegently determine more difficult settings
    if not config.has_option("DATABASE", "user"):
        if "USERNAME" in __os.environ:  # windows
            user = __os.environ["USERNAME"]
        elif "USER" in __os.environ:  # unix
            user = __os.environ["USER"]
        config.set("DATABASE", "user", user)
    if not config.has_option("EXECUTABLES", "psql"):
        psql = which("psql91")
        if psql is None:
            psql = which("psql")
        if psql is None:
            psql = "SET_PATH_TO_PSQL"
        config.set("EXECUTABLES", "psql", psql)
    if not config.has_option("EXECUTABLES", "R"):
        R = which("R")
        if R is None:
            R = "SET_PATH_TO_R"
        config.set("EXECUTABLES", "R", R)
    if not config.has_option("EXECUTABLES", "Rscript"):
        Rscript = which("Rscript")
        if Rscript is None:
            Rscript = "SET_PATH_TO_Rscript"
        config.set("EXECUTABLES", "Rscript", Rscript)
    if not config.has_option("EXECUTABLES", "primer3"):
        primer3 = which("primer3_core")
        if primer3 is None:
            primer3 = "./primer3_core"
        config.set("EXECUTABLES", "primer3", primer3)

    # write the options back to the file
    with open(filepath, "w") as outfile:
        config.write(outfile)

    # save options as variables
    self.postgres_user = config.get("DATABASE", "postgres_user")
    self.postgres_password = config.get("DATABASE", "postgres_password")
    if len(self.postgres_password) > 0:
        __os.environ["PGPASSWORD"] = self.postgres_password
    self.postgres_database = config.get("DATABASE", "postgres_database")
    self.postgres_host = config.get("DATABASE", "postgres_host")
    self.psql = _escape_space(config.get("EXECUTABLES", "psql"))
    self.R = _escape_space(config.get("EXECUTABLES", "R"))
    self.Rscript = _escape_space(config.get("EXECUTABLES", "Rscript"))
    self.primer3 = _escape_space(config.get("EXECUTABLES", "primer3"))
    self.cufflinks = config.get("EXECUTABLES", "cufflinks")

    # make a psql string with the database options included
    self.hostname, self.port = postgres_host.split(":")
    self.psql_full = "%s --host=%s --username=%s --port=%s " % \
        (self.psql, self.hostname, self.postgres_user, self.port)
    self.entrez_email = config.get("MISC", "entrez_email")
    if self.entrez_email == "SET_ENTREZ_EMAIL": self.entrez_email = None
    #set home directory
    self.home_directory = config.get("MISC", "home_directory")
    self.data_directory = config.get("MISC", "data_directory")


load_settings_from_file()
del SafeConfigParser, modules

_base_site_file = \
"""WSGIScriptAlias /om OMLIB_DIRserver.py
<Directory OMLIB_DIR>
    WSGIScriptReloading On
    Order allow,deny
    Allow from 132.239
    Allow from 128.54
    Allow from 127
    Allow from 137.110
</Directory>
"""

def write_apache_site_file():
    with open(om_directory + "om.site", "w") as outfile:
        outfile.write(_base_site_file.replace("OMLIB_DIR", omlib_directory))

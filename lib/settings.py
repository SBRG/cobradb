"""retrive local user settings"""

from ConfigParser import SafeConfigParser
import os as __os
from os.path import split as __split, join as __join, abspath as __abspath, \
    isfile as __isfile
from sys import modules

self = modules[__name__]

# define various filepaths
trnlib_directory = __join(__split(__abspath(__file__))[0], "")
trn_directory = __join(__abspath(__join(trnlib_directory, "..")), "")

def which(program):
    """returns path to an executable if it is found in the path"""
    fpath, fname = __split(program)
    if fpath:
        if __isfile(program) and __os.access(program, __os.X_OK):
            return program
    else:
        paths_to_search = __os.environ["PATH"].split(__os.pathsep)
        paths_to_search.extend((trnlib_directory, trn_directory))
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
config.set("DATABASE", "host", "localhost:5432")
config.set("DATABASE", "database", "ome")
config.set("DATABASE", "password", "")
config.set("DATABASE", "schema", "ecoli")

config.add_section("MISC")
config.set("MISC", "entrez_email", "SET_ENTREZ_EMAIL")
try:
    default_home = __os.environ["HOME"]
except:
    try:
        default_home = __os.environ["USERPROFILE"]
    except:
        default_home = ""
default_dropbox = __join(default_home, "Dropbox", "")
# the data_directory is by default in dropbox
default_data_dir = __join(default_dropbox, "trn_data", "")
config.set("MISC", "home_directory", default_home)
config.set("MISC", "dropbox_directory", default_dropbox)
config.set("MISC", "data_directory", default_data_dir)
del default_home, default_dropbox, default_data_dir
config.add_section("EXECUTABLES")

# overwrite defaults settings with settings from the file
def load_settings_from_file(filepath="settings.ini", in_trnlib=True):
    """reload settings from a different settings file
    
    filepath: The path to the settings file to use
    
    in_trnlib: Whether or not the path given is a relative path from the trnlib
        directory"""
    if in_trnlib:
        filepath = __join(trnlib_directory, filepath)
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
    self.user = config.get("DATABASE", "user")
    self.password = config.get("DATABASE", "password")
    if len(self.password) > 0:
        __os.environ["PGPASSWORD"] = self.password
    self.host = config.get("DATABASE", "host")
    self.database = config.get("DATABASE", "postgres_database")
    self.psql = _escape_space(config.get("EXECUTABLES", "psql"))
    self.R = _escape_space(config.get("EXECUTABLES", "R"))
    self.Rscript = _escape_space(config.get("EXECUTABLES", "Rscript"))
    self.primer3 = _escape_space(config.get("EXECUTABLES", "primer3"))
    # make a psql string with the database options included
    self.hostname, self.port = host.split(":")
    self.psql_full = "%s --host=%s --username=%s --port=%s " % \
        (self.psql, self.hostname, self.user, self.port)
    self.entrez_email = config.get("MISC", "entrez_email")
    if self.entrez_email == "SET_ENTREZ_EMAIL": self.entrez_email = None
    #set home directory
    self.home_directory = config.get("MISC", "home_directory")
    self.data_directory = config.get("MISC", "data_directory")
    
load_settings_from_file()
del SafeConfigParser, modules

_base_site_file = \
"""WSGIScriptAlias /trn TRNLIB_DIRserver.py
<Directory TRNLIB_DIR>
    WSGIScriptReloading On
    Order allow,deny
    Allow from 132.239
    Allow from 128.54
    Allow from 127
    Allow from 137.110
</Directory>
"""

def write_apache_site_file():
    with open(trn_directory + "trn.site", "w") as outfile:
        outfile.write(_base_site_file.replace("TRNLIB_DIR", trnlib_directory))

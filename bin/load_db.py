#! /usr/bin/env python

# configure the logger before imports so other packages do not override this
# setup
import logging
import time
def configure_logger(log_file=None, level=logging.INFO, overwrite_log=True,
                     format=logging.BASIC_FORMAT):
    # console and file
    if log_file is None:
        logging.basicConfig(stream=sys.stdout, level=level, format=format)
    else:
        logging.basicConfig(filename=log_file, level=level,
                            filemode=('w' if overwrite_log else 'a'),
                            format=format)
        console = logging.StreamHandler()
        console.setLevel(level)
        console.setFormatter(logging.Formatter(format))
        logging.getLogger("").addHandler(console)
configure_logger('%s OME load_db.log' % time.strftime('%Y-%m-%d %H:%M:%S'),
                 level=logging.INFO)


from ome import base, settings, components, datasets, models, timing, util
from ome.loading import AlreadyLoadedError
from ome.loading import dataset_loading
from ome.loading import component_loading
from ome.loading.component_loading import BadGenomeError
from ome.loading import model_loading
from ome.loading import map_loading

from sqlalchemy.schema import Sequence, CreateSequence
from warnings import warn
import sys
import os
from os.path import join
import argparse

try:
    from pymongo import ASCENDING
    MONGO_INSTALLED = True
except ImportError:
    logging.warn('pymongo not installed')
    MONGO_INSTALLED = False

parser = argparse.ArgumentParser()
parser.add_argument('--drop-all', help='Empty database and reload data. NOTE: Does not drop types (e.g. enum categories)', action='store_true')
parser.add_argument('--drop-models', help='Empty model and map data', action='store_true')
parser.add_argument('--drop-maps', help='Empty map data', action='store_true')
parser.add_argument('--skip-genomes', help='Skip genome loading', action='store_true')
parser.add_argument('--skip-models', help='Skip model loading', action='store_true')

args = parser.parse_args()

def drop_all_tables(engine, enums_to_drop=None):
    """Drops all tables and, optionally, user enums from a postgres database.

    Adapted from: http://www.siafoo.net/snippet/85


    NOTE: To list all the user defined enums, use something like this:

        SELECT DISTINCT t.typname
        FROM pg_catalog.pg_type t
        JOIN pg_catalog.pg_enum e ON t.oid = e.enumtypid;

    """

    from sqlalchemy.sql.expression import text

    table_sql = ("SELECT table_name FROM information_schema.tables "
                 "WHERE table_schema='public' AND table_name NOT LIKE 'pg_%%'")

    for table in [name for (name, ) in engine.execute(text(table_sql))]:
        engine.execute(text('DROP TABLE %s CASCADE' % table))

    # drop the enum types
    if enums_to_drop is not None:
        for enum in enums_to_drop:
            engine.execute(text('DROP TYPE IF EXISTS %s CASCADE' % enum))

if __name__ == "__main__":
    if args.drop_all:
        logging.info("Dropping everything from the database")
        drop_all_tables(base.engine, base.custom_enums.keys())

        try:
            base.omics_database.genome_data.drop()
        except:
            pass

    logging.info("Building the database models")
    base.Base.metadata.create_all()

    if args.drop_models:
        logging.info('Dropping rows from models')
        connection = base.engine.connect()
        trans = connection.begin()
        try:
            connection.execute('TRUNCATE model, reaction, component, compartment CASCADE;')
            trans.commit()
        except:
            trans.rollback()

    # make the session
    session = base.Session()

    # get the models and genomes
    model_dir = join(settings.data_directory, 'models')
    model_genome_path = settings.model_genome
    logging.info('Loading models and genomes using %s' % model_genome_path)
    lines = util.load_tsv(model_genome_path, required_column_num=3)
    models_list = []; found_genomes = {}
    for line in lines:
        val_nones = [(x if x.strip() != '' else None) for x in line]
        model_filename, bioproject_id, pub_ref = val_nones
        if bioproject_id is not None:
            found_genomes[bioproject_id] = False
        models_list.append((model_filename, bioproject_id, pub_ref))


    if not args.skip_genomes:
        logging.info('Finding matching GenBank files')
        genbank_dir = join(settings.data_directory, 'annotation', 'genbank')
        # get the genbank files with matching bioproject ids
        genome_files = []
        for filename in os.listdir(genbank_dir):
            logging.debug('Looking for BioProject ID in %s' % filename)
            try:
                bioproject_id, _ = component_loading.get_bioproject_id(join(genbank_dir, filename),
                                                                       fast=True)
            except BadGenomeError as e:
                logging.warn(e.message)
                continue
            if bioproject_id in found_genomes:
                found_genomes[bioproject_id] = True
                genome_files.append(filename)
        not_found = [k for k, v in found_genomes.iteritems() if v is False]
        if len(not_found) > 0:
            logging.warn('No genbank file for for %s' % not_found)

        logging.info('Loading genomes')
        n = len(genome_files)
        for i, genbank_file in enumerate(genome_files):
            logging.info('Loading genome from genbank file (%d of %d) %s' % (i + 1, n, genbank_file))
            try:
                component_loading.load_genome(join(genbank_dir, genbank_file), session)
            except Exception as e:
                logging.exception(e)

        # chromosome gffs
        data_genomes = (session
                        .query(base.Genome)
                        .filter(base.Genome.bioproject_id.in_(['PRJNA57779']))
                        .all())
        raw_flag = False
        normalize_flag = False
        for genome in data_genomes:
            for chromosome in genome.chromosomes:
                component_loading.write_chromosome_annotation_gff(base, components,
                                                                  chromosome)

    if not args.skip_models:
        logging.info("Loading models")
        n = len(models_list)
        for i, (model_filename, bioproject_id, pub_ref) in enumerate(models_list):
            logging.info('Loading model (%d of %d) %s' % (i + 1, n, model_filename))
            try:
                model_loading.load_model(join(model_dir, model_filename),
                                         bioproject_id, pub_ref, session)
            except AlreadyLoadedError as e:
                logging.info(e.message)
            except Exception as e:
                logging.error('Could not load model %s.' % model_filename)
                logging.exception(e)

    logging.info("Loading Escher maps")
    map_loading.load_maps_from_server(session, drop_maps=(args.drop_models or
                                                          args.drop_maps))

    if MONGO_INSTALLED and base.omics_database is not None:
        genome_data = base.omics_database.genome_data
        genome_data.create_index([("data_set_id", ASCENDING),
                                  ("leftpos", ASCENDING)])

    session.close()
    base.Session.close_all()

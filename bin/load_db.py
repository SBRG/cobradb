#! /usr/bin/env python

# configure the logger before imports so other packages do not override this
# setup
import logging
import time
import sys

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
from ome.loading.component_loading import BadGenomeError, get_genbank_accessions
from ome.loading import model_loading
from ome.loading import map_loading

import os
from os import listdir
from os.path import join
import argparse
from collections import defaultdict


parser = argparse.ArgumentParser()
parser.add_argument('--drop-all', help='Empty database and reload data. NOTE: Does not drop types (e.g. enum categories)', action='store_true')
parser.add_argument('--drop-models', help='Empty model and map data', action='store_true')
parser.add_argument('--drop-maps', help='Empty map data', action='store_true')
parser.add_argument('--skip-genomes', help='Skip genome loading', action='store_true')
parser.add_argument('--skip-models', help='Skip model loading', action='store_true')
parser.add_argument('--skip-maps', help='Skip map loading', action='store_true')

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
    model_dir = settings.model_directory
    model_genome_path = settings.model_genome
    logging.info('Loading models and genomes using %s' % model_genome_path)
    lines = util.load_tsv(model_genome_path, required_column_num=3)
    models_list = []
    for line in lines:
        model_filename, pub_ref_string, genome_ref_string = line
        if genome_ref_string is None:
            genome_ref = None
        else:
            genome_ref = util.ref_str_to_tuple(genome_ref_string)
        if pub_ref_string is None:
            pub_ref = None
        else:
            pub_ref = util.ref_str_to_tuple(pub_ref_string)
        models_list.append({'model_filename': model_filename,
                            'pub_ref': pub_ref,
                            'genome_ref': genome_ref})

    genomes_for_models = {}
    if not args.skip_genomes:
        # find the accessions for the genbank files
        logging.info('Finding GenBank files')
        refseq_dir = settings.refseq_directory
        # unique refs
        genome_refs = {r['genome_ref'] for r in models_list}
        # loop through all the files
        genome_file_locations = defaultdict(list)
        for refseq_filename in listdir(refseq_dir):
            if refseq_filename.startswith('.'):
                continue
            refseq_filepath = join(refseq_dir, refseq_filename)
            # check both accession and assembly for a match
            ids = get_genbank_accessions(refseq_filepath, fast=True)
            # if the ids couldn't be found
            if all(x is None for x in ids.values()):
                logging.warn('Could not find accessions for genbank file %s' % refseq_filepath)
                continue
            # look for matching ids from the model-genome file
            for genome_ref in ids.iteritems():
                if genome_ref in genome_refs:
                    genome_file_locations[genome_ref].append(refseq_filepath)

        # load the genomes
        n = len(genome_refs)
        for i, genome_ref in enumerate(genome_refs):
            logging.info('Loading genome ({} of {}) {}'.format(i + 1, n, genome_ref))
            file_paths = genome_file_locations[genome_ref]
            try:
                component_loading.load_genome(genome_ref, file_paths, session)
            except Exception as e:
                logging.exception(e)

    if not args.skip_models:
        logging.info("Loading models")
        n = len(models_list)
        model_dir = settings.model_directory
        for i, model_dict in enumerate(models_list):
            logging.info('Loading model ({} of {}) {}'
                         .format(i + 1, n, model_dict['model_filename']))
            try:
                model_loading.load_model(join(model_dir, model_dict['model_filename']),
                                         model_dict['pub_ref'],
                                         model_dict['genome_ref'],
                                         session)
            except AlreadyLoadedError as e:
                logging.info(str(e))
            except Exception as e:
                logging.error('Could not load model %s.' % model_filename)
                logging.exception(e)

    if not args.skip_maps:
        logging.info("Loading Escher maps")
        map_loading.load_maps_from_server(session, drop_maps=(args.drop_models or
                                                              args.drop_maps))

    session.close()
    base.Session.close_all()

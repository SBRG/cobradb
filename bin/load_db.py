#! /usr/bin/python

from ome import base, settings, components, datasets, models, timing
from ome.loading import dataset_loading
from ome.loading import component_loading
from ome.loading import model_loading

from sqlalchemy.schema import Sequence, CreateSequence
from warnings import warn
import sys
import os
from os.path import join
import argparse
import logging

# configure the logger
logging.basicConfig(level=logging.INFO, stream=sys.stdout)

try:
    from pymongo import ASCENDING
    MONGO_INSTALLED = True
except ImportError:
    warn('pymongo not installed')
    MONGO_INSTALLED = False

parser = argparse.ArgumentParser()
parser.add_argument('--drop-all', help='Empty database and reload data', action='store_true')
parser.add_argument('--drop-models', help='Empty model data', action='store_true')
parser.add_argument('--skip-genomes', help='Skip genome loading', action='store_true')
parser.add_argument('--skip-models', help='Skip model loading', action='store_true')

args = parser.parse_args()

def drop_all_tables(engine):
    """Drops all tables from a postgres database.

    Adapted from: http://www.siafoo.net/snippet/85

    """
    
    from sqlalchemy.sql.expression import text
     
    table_sql = ("SELECT table_name FROM information_schema.tables "
                 "WHERE table_schema='public' AND table_name NOT LIKE 'pg_%%'")

    for table in [name for (name, ) in engine.execute(text(table_sql))]:
        engine.execute(text('DROP TABLE %s CASCADE' % table))

if __name__ == "__main__":
    if args.drop_all:
        logging.info("Dropping everything from the database")
        drop_all_tables(base.engine)
        logging.info("Building the database models")
        base.Base.metadata.create_all()

        try:
            base.omics_database.genome_data.drop()
        except: 
            pass
    
    if args.drop_models:
        logging.info('Dropping rows from models')
        connection = base.engine.connect()
        trans = connection.begin()
        try:
            connection.execute('TRUNCATE model, reaction, component, compartment CASCADE;')
            trans.commit()
        except:
            trans.rollback()
                        
    if not args.skip_genomes:
        logging.info('Loading genomes')
        genbank_dir = join(settings.data_directory, 'annotation', 'genbank')
        dirs = os.listdir(genbank_dir)
        n = len(dirs)
        for i, genbank_file in enumerate(dirs):
            logging.info('Loading genome from genbank file (%d of %d) %s' % (i + 1, n, genbank_file))
            try:
                if genbank_file != '.DS_Store':
                    component_loading.load_genome(join(genbank_dir, genbank_file))
            except Exception as e:
                logging.error(str(e))

        session = base.Session()

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
        model_dir = join(settings.data_directory, 'models')
        model_genome_file = join(settings.data_directory,
                                 'annotation',
                                 'model-genome.txt')
        with open(model_genome_file, 'r') as f:
            lines = f.readlines()
            n = len(lines)
            for i, line in enumerate(lines):
                model_id, genome_id, timestamp, pmid = line.rstrip('\n').split(',')
                logging.info('Loading model (%d of %d) %s' % (i + 1, n, model_id))
                try:
                    model_loading.load_model(model_id, model_dir, genome_id,
                                             timestamp, pmid)
                except Exception as e:
                    logging.error('Could not load model %s. %s' % (model_id, e))

    genome_data = base.omics_database.genome_data

    if MONGO_INSTALLED:
        genome_data.create_index([("data_set_id", ASCENDING),
                                  ("leftpos", ASCENDING)])

    session.close()

"""
        dataset_loading.load_raw_files(settings.data_directory+'/chip_experiment/bam/crp', group_name='crp', normalize=normalize_flag, raw=raw_flag)
        dataset_loading.load_raw_files(settings.data_directory+'/chip_experiment/bam/yome', group_name='yome', normalize=normalize_flag, raw=raw_flag)

        dataset_loading.load_raw_files(settings.data_directory+'/chip_experiment/gff', group_name='trn', normalize=False, raw=raw_flag)

        dataset_loading.load_raw_files(settings.data_directory+'/rnaseq_experiment/fastq/crp', group_name='crp', normalize=False, raw=False)
        dataset_loading.load_raw_files(settings.data_directory+'/rnaseq_experiment/fastq/yome', group_name='yome', normalize=False, raw=False)
        #dataset_loading.load_raw_files(settings.data_directory+'/rnaseq_experiment/bam', normalize=True)
        #dataset_loading.load_raw_files(settings.data_directory+'/chip_experiment/bam', normalize=False)
        dataset_loading.load_raw_files(settings.data_directory+'/microarray/asv2', group_name='asv2', raw=False)
        dataset_loading.load_raw_files(settings.data_directory+'/microarray/ec2', group_name='ec2', raw=False)


        experiment_sets = dataset_loading.query_experiment_sets()
        dataset_loading.load_experiment_sets(experiment_sets)


        for chromosome in genome.chromosomes:
            component_loading.load_metacyc_proteins(base, components, chromosome)
            component_loading.load_metacyc_bindsites(base, components, chromosome)
            component_loading.load_metacyc_transcription_units(base, components, chromosome)

        old_gff_file = settings.data_directory+'/annotation/NC_000913.2_old.gff'

        #dataset_loading.run_cuffquant(base, datasets, genome, group_name='crp')
        #dataset_loading.run_cuffnorm(base, datasets, genome, group_name='crp', gff_file=old_gff_file, overwrite=True)
        #dataset_loading.run_cuffnorm(base, datasets, genome, group_name='yome', overwrite=True)
        #dataset_loading.run_cuffdiff(base, datasets, genome, group_name='crp', gff_file=old_gff_file, overwrite=True)
        #dataset_loading.run_cuffdiff(base, datasets, genome, group_name='yome', overwrite=True)
        #dataset_loading.run_gem(base, datasets, genome, debug=True)


        dataset_loading.load_gem(session.query(ChIPPeakAnalysis).all(), base, datasets, genome)
        dataset_loading.load_gff_chip_peaks(session.query(ChIPPeakAnalysis).all(), base, datasets, genome, group_name='gff-BK')

        dataset_loading.load_extra_analyses(base, datasets, genome, settings.data_directory+'/ChIP_peaks/gps-curated-HL28Aug14', group_name='gps-curated-HL28Aug14')
        dataset_loading.load_gff_chip_peaks(session.query(ChIPPeakAnalysis).all(), base, datasets, genome, group_name='gps-curated-HL28Aug14')

        component_loading.load_kegg_pathways(base, components)

        dataset_loading.load_cuffnorm(base, datasets, group_name='crp')
        dataset_loading.load_cuffnorm(base, datasets, group_name='yome')
        dataset_loading.load_cuffdiff(group_name='crp')
        dataset_loading.load_cuffdiff(group_name='yome')

        dataset_loading.load_arraydata(settings.data_directory+'/microarray/formatted_asv2.txt', group_name='asv2')
        dataset_loading.load_arraydata(settings.data_directory+'/microarray/formatted_ec2.txt', group_name='ec2')

        dataset_loading.run_array_ttests(base, datasets, genome, group_name='asv2')
        dataset_loading.run_array_ttests(base, datasets, genome, group_name='ec2')

        dataset_loading.make_genome_region_map(base, datasets, genome)
        """

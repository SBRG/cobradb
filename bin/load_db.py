#! /usr/bin/python

from ome import base,settings,components,datasets,models,timing

from ome.loading import dataset_loading
from ome.loading import component_loading
from ome.loading import model_loading

from sqlalchemy.schema import Sequence,CreateSequence
from warnings import warn
import sys
import os
import argparse

try:
    from pymongo import ASCENDING
    MONGO_INSTALLED = True
except ImportError:
    warn('pymongo not installed')
    MONGO_INSTALLED = False

parser = argparse.ArgumentParser()
parser.add_argument("--dropall", help="will empty database and reload data", action="store_true")
parser.add_argument("--dropmodels", help="will empty model data", action="store_true")

args = parser.parse_args()

if __name__ == "__main__":
    """
    #if not dataset_loading.query_yes_no('This will drop the ENTIRE database and load from scratch, ' + \
                        'are you sure you want to do this?'): sys.exit()
    """
    
    
    if args.dropall:
        base.Base.metadata.drop_all()
        base.Base.metadata.create_all()

        try: base.omics_database.genome_data.drop()
        except: None
        print "dropping all"
        #base.engine.execute(CreateSequence(Sequence('wids')))
    
    if args.dropmodels:
        print "dropping rows from models"
        connection = base.engine.connect()
        trans = connection.begin()
        try:
            connection.execute('TRUNCATE model,reaction,component,compartment CASCADE;')
            trans.commit()
        except:
            trans.rollback()
                        
    for genbank_file in os.listdir(settings.data_directory+'annotation/genbank'):
        #if genbank_file not in ['NC_000913.2.gb']: continue

        component_loading.load_genome(genbank_file, base, components, debug=False)


    session = base.Session()

    data_genomes = session.query(base.Genome).filter(base.Genome.bioproject_id.in_(['PRJNA57779'])).all()


    raw_flag = False
    normalize_flag = False

    for genome in data_genomes:

        for chromosome in genome.chromosomes:
            component_loading.write_chromosome_annotation_gff(base, components, chromosome)
    
    with open(settings.data_directory+'/annotation/model-genome.txt') as file:
        for line in file:
            model_id,genome_id,model_creation_timestamp,pmid = line.rstrip('\n').split(',')
            try:
                model_loading.load_model(model_id, genome_id, model_creation_timestamp, pmid)
            except Exception as e:
                warn('Could not load model %s. %s' % (model_id, e))

    genome_data = base.omics_database.genome_data

    if MONGO_INSTALLED:
        genome_data.create_index([("data_set_id", ASCENDING), ("leftpos", ASCENDING)])

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

        #dataset_loading.run_cuffquant(base, datasets, genome, group_name='crp', debug=False)
        #dataset_loading.run_cuffnorm(base, datasets, genome, group_name='crp', gff_file=old_gff_file, debug=False, overwrite=True)
        #dataset_loading.run_cuffnorm(base, datasets, genome, group_name='yome', debug=False, overwrite=True)
        #dataset_loading.run_cuffdiff(base, datasets, genome, group_name='crp', gff_file=old_gff_file, debug=False, overwrite=True)
        #dataset_loading.run_cuffdiff(base, datasets, genome, group_name='yome', debug=False, overwrite=True)
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

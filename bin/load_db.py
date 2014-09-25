from om import base,settings,components,data,timing
from om.data import *
from om.model import *
from om.loading import data_loading
from om.loading import component_loading

from sqlalchemy import func
import simplejson as json
import numpy as np
import os,sys,math,shutil,subprocess

from IPython import embed



@timing
def load_raw_files(directory_path, group_name='default', normalize=True, overwrite=False, raw=True):
    """This will load raw .bam files, .gff files, and affymetrix .CEL files into the associated
    genome_data tables.  A genome object must also be supplied for proper mappping. The raw signal
    will go into a mongoDB genome_data table and the rest will go into corresponding SQL tables.
    This function will by default assume name based loading, however, it can also be configured
    to read a metadata file with corresponding experiment details.
    """
    session = base.Session()

    parseable_file_extensions = ['bam','sorted.bam','fastq','fastq.gz','R1.fastq.gz','R2.fastq.gz','CEL','gff']
    experiments = []

    fastq_experiment_paths = {}

    for file_name in os.listdir(directory_path):
        if file_name[0:4] == 'chip': file_name = 'ChIP'+file_name[4:] #scrub misnamed files

        file_prefix = file_name.split('.')[0]
        file_suffix = '.'.join(file_name.split('.')[1:])


        """if an inner directory is found, load that into a separate experiment group"""
        if os.path.isdir(directory_path+'/'+file_name):

            load_raw_files(directory_path+'/'+file_name, group_name=file_name, normalize=False)


        if file_suffix not in parseable_file_extensions: continue

        """if the experiment was split across multiple sequencing runs the fastq files will be split
           and they will be preceded by an integer e.g. 1_exp_name, 2_exp_name, etc...
        """

        if file_name[0] in ['1','2','3']:
          file_prefix,fastq_paths = get_split_fastq_files(file_name, directory_path)
          fastq_experiment_paths[file_prefix] = fastq_paths

        elif file_suffix in ['fastq','fastq.gz','R1.fastq.gz','R2.fastq.gz']:
          #print file_prefix
          try: fastq_experiment_paths[file_prefix].append(file_name)
          except: fastq_experiment_paths[file_prefix] = [file_name]

        experiment = data_loading.create_name_based_experiment(session, file_prefix, group_name)

        experiments.append(experiment)


    for experiment in set(experiments):
      if experiment.name in fastq_experiment_paths.keys():
          data_loading.run_bowtie2(experiment, fastq_experiment_paths[experiment.name], overwrite=overwrite, debug=True)


    if normalize:
        normalization_factors = data_loading.calculate_normalization_factors(experiments[0].type, group_name)
    else:
        normalization_factors = {exp.name: 1. for exp in experiments}


    for experiment in set(experiments):
        if not raw: continue

        if experiment.type == 'chip_experiment':
            norm_factor = normalization_factors[experiment.name]
            if experiment.protocol_type == 'ChIPExo':
                data_loading.load_raw_experiment_data(experiment, loading_cutoff=10, flip=False, five_prime=True, norm_factor=norm_factor)
            elif experiment.protocol_type == 'ChIPchip':
                data_loading.load_raw_gff_to_db(experiment)

        elif experiment.type == 'rnaseq_experiment':
            norm_factor = normalization_factors[experiment.name]
            data_loading.load_raw_experiment_data(experiment, loading_cutoff=10., flip=True, five_prime=False, norm_factor=norm_factor)

    session.close()


def get_split_fastq_files(file_name, directory_path):
    fastq_files = []

    file_prefix = file_name.split('.')[0]
    file_suffix = '.'.join(file_name.split('.')[1:])

    #suffix = '_'.join(file_name.split('_')[1:])

    for file_name2 in os.listdir(directory_path):
        file_prefix2 = file_name2.split('.')[0]
        if file_name2[0] in ['1','2','3'] and file_prefix[2:] == file_prefix2[2:]:
            fastq_files.append(file_name2)
    return file_prefix[2:],fastq_files


def query_experiment_sets():
    """This queries all the replicates of an experiment and groups them for further processing"""
    experiment_sets = {}

    ome = base.Session()

    experiment_sets['RNAseq'] = ome.query(func.array_agg(RNASeqExperiment.name),RNASeqExperiment.group_name).\
                                          group_by(RNASeqExperiment.group_name, RNASeqExperiment.strain_id,
                                                   RNASeqExperiment.environment_id, RNASeqExperiment.machine_id,
                                                   RNASeqExperiment.sequencing_type).all()

    experiment_sets['array'] = ome.query(func.array_agg(ArrayExperiment.name)).\
                                         group_by(ArrayExperiment.strain_id, ArrayExperiment.environment_id,\
                                         ArrayExperiment.platform).all()

    experiment_sets['ChIP'] = ome.query(func.array_agg(ChIPExperiment.name)).\
                                        group_by(ChIPExperiment.strain_id, ChIPExperiment.environment_id,\
                                        ChIPExperiment.antibody, ChIPExperiment.protocol_type,\
                                        ChIPExperiment.target).all()
    ome.close()
    return experiment_sets


def load_experiment_sets(experiment_sets):
    """This will create the database entries for the grouped experiments"""

    ome = base.Session()

    for exp_group in experiment_sets['RNAseq']:
        vals = exp_group[0][0].split('_')
        if len(vals) > 6: exp_group_name = '_'.join(vals[0:5]+vals[-1:])
        else: exp_group_name = '_'.join(vals[0:5])

        exp = ome.query(RNASeqExperiment).filter_by(name=exp_group[0][0], group_name=exp_group[1]).one()
        exp_analysis = ome.get_or_create(NormalizedExpression, replicate=1, name=exp_group_name, environment_id=exp.environment.id,\
                                         strain_id=exp.strain.id, group_name=exp.group_name, expression_type='rnaseq_experiment')
        for exp_name in exp_group[0]:
            expt = ome.query(RNASeqExperiment).filter_by(name=exp_name, group_name=exp.group_name).one()
            ome.get_or_create(AnalysisComposition, analysis_id = exp_analysis.id, data_set_id = expt.id)

    for exp_group in experiment_sets['array']:
        exp_group_name = '_'.join(exp_group[0][0].split('_')[:-2])
        exp = ome.query(ArrayExperiment).filter_by(name=exp_group[0][0]).one()
        exp_analysis = ome.get_or_create(NormalizedExpression, replicate=1, name=exp_group_name, environment_id=exp.environment.id,\
                                         strain_id=exp.strain.id, group_name=exp.group_name, expression_type='array_experiment')
        for exp_name in exp_group[0]:
            exp = ome.query(ArrayExperiment).filter_by(name=exp_name).one()
            ome.get_or_create(AnalysisComposition, analysis_id = exp_analysis.id, data_set_id = exp.id)

    default_parameters = {'mrc':20, 'smooth':3, 'nrf':'', 'outNP':'', 'nf':'', 'k_min': 4, 'k_max': 22, 'k_win':150}

    for exp_group in experiment_sets['ChIP']:
        parameters = {'mrc':20, 'smooth':3, 'outNP':'', 'nrf':'', 'nf':'','k_min': 4, 'k_max': 22, 'k_win':150}
        if not set(parameters.items()) - set(default_parameters.items()):
            parameter_name = 'default'
        else:
            parameter_name = '_'.join([y+'-'+str(z) for y,z in dict(set(parameters.items()) - set(default_parameters.items())).iteritems()])


        vals = exp_group[0][0].split('_')
        exp_group_name = '_'.join(vals[0:5]+vals[6:]+[parameter_name,'peaks'])


        exp = ome.query(ChIPExperiment).filter_by(name=exp_group[0][0]).one()


        exp_analysis = ome.get_or_create(ChIPPeakAnalysis, name=exp_group_name, environment_id=exp.environment.id, strain_id=exp.strain.id,\
                                                           parameters=json.dumps(parameters), replicate=1, group_name=exp.group_name)
        for exp_name in exp_group[0]:
            expt = ome.query(ChIPExperiment).filter_by(name=exp_name, group_name=exp.group_name).one()
            ome.get_or_create(AnalysisComposition, analysis_id = exp_analysis.id, data_set_id = expt.id)

    ome.close()


if __name__ == "__main__":

    #if not query_yes_no('This will drop the ENTIRE database and load from scratch, ' + \
    #                    'are you sure you want to do this?'): sys.exit()


    base.Base.metadata.drop_all()
    base.omics_database.genome_data.drop()
    base.Base.metadata.create_all()

    component_loading.load_genomes(base, components)


    session = base.Session()

    data_genomes = session.query(base.Genome).filter(base.Genome.ncbi_id.in_(['NC_000913.2'])).all()


    raw_flag = True
    normalize_flag = True

    for genome in data_genomes:

        component_loading.write_genome_annotation_gff(base, components, genome)

        load_raw_files(settings.data_directory+'/chip_experiment/bam/crp', group_name='crp', normalize=normalize_flag, raw=raw_flag)
        load_raw_files(settings.data_directory+'/chip_experiment/bam/yome', group_name='yome', normalize=normalize_flag, raw=raw_flag)

        load_raw_files(settings.data_directory+'/chip_experiment/gff', group_name='trn', normalize=False, raw=raw_flag)

        load_raw_files(settings.data_directory+'/rnaseq_experiment/fastq/crp', group_name='crp', normalize=False, raw=False)
        load_raw_files(settings.data_directory+'/rnaseq_experiment/fastq/yome', group_name='yome', normalize=False, raw=False)
        #load_raw_files(settings.data_directory+'/rnaseq_experiment/bam', normalize=True)
        #load_raw_files(settings.data_directory+'/chip_experiment/bam', normalize=False)
        load_raw_files(settings.data_directory+'/microarray/asv2', group_name='asv2', raw=False)
        load_raw_files(settings.data_directory+'/microarray/ec2', group_name='ec2', raw=False)


        experiment_sets = query_experiment_sets()
        load_experiment_sets(experiment_sets)


        component_loading.load_metacyc_proteins(base, components, genome)
        component_loading.load_metacyc_bindsites(base, components, genome)
        component_loading.load_metacyc_transcription_units(base, components, genome)

        old_gff_file = settings.data_directory+'/annotation/NC_000913.2_old.gff'

        #data_loading.run_cuffquant(base, data, genome, group_name='crp', debug=False)
        #data_loading.run_cuffnorm(base, data, genome, group_name='crp', gff_file=old_gff_file, debug=False, overwrite=True)
        #data_loading.run_cuffnorm(base, data, genome, group_name='yome', debug=False, overwrite=True)
        #data_loading.run_cuffdiff(base, data, genome, group_name='crp', gff_file=old_gff_file, debug=False, overwrite=True)
        #data_loading.run_cuffdiff(base, data, genome, group_name='yome', debug=False, overwrite=True)
        #data_loading.run_gem(base, data, genome, debug=True)


        data_loading.load_gem(session.query(ChIPPeakAnalysis).all(), base, data, genome)
        data_loading.load_gff_chip_peaks(session.query(ChIPPeakAnalysis).all(), base, data, genome, group_name='gff-BK')

        data_loading.load_extra_analyses(base, data, genome, settings.data_directory+'/ChIP_peaks/gps-curated-HL28Aug14', group_name='gps-curated-HL28Aug14')
        data_loading.load_gff_chip_peaks(session.query(ChIPPeakAnalysis).all(), base, data, genome, group_name='gps-curated-HL28Aug14')

        component_loading.load_kegg_pathways(base, components)

        data_loading.load_cuffnorm(base, data, group_name='crp')
        data_loading.load_cuffnorm(base, data, group_name='yome')
        data_loading.load_cuffdiff(group_name='crp')
        data_loading.load_cuffdiff(group_name='yome')

        data_loading.load_arraydata(settings.data_directory+'/microarray/formatted_asv2.txt', group_name='asv2')
        data_loading.load_arraydata(settings.data_directory+'/microarray/formatted_ec2.txt', group_name='ec2')

        data_loading.run_array_ttests(base, data, genome, group_name='asv2')
        data_loading.run_array_ttests(base, data, genome, group_name='ec2')

        data_loading.make_genome_region_map(base, data, genome)


    genome_data = base.omics_database.genome_data

    genome_data.create_index([("data_set_id",ASCENDING), ("leftpos", ASCENDING)])

    session.close()

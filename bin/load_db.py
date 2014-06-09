from om import base,settings,components,timing
from om.data import *
from om.loading import data_loading
from om.loading import component_loading

from sqlalchemy import func
import simplejson as json
import numpy as np
import os,sys,math,shutil,subprocess



def calculate_normalization_factors():
    exp_types = ['ChIP','RNAseq']
    mapped_read_norm_factor = {}
    for exp_type in exp_types:
        mapped_reads = {}
        directory = settings.dropbox_directory+'/crp/data/'+exp_type+'/bam'
        for file_name in os.listdir(directory):
            if file_name[-3:] != 'bam': continue
            #mapped_reads[file_name] = int(subprocess.check_output(['samtools', 'flagstat', directory+'/'+file_name]).split('\n')[2].split()[0])
            mapped_reads[file_name] = 30000000
        mean_read_count = np.array(mapped_reads.values()).mean()
        for exp,value in mapped_reads.iteritems():
            #mapped_read_norm_factor[exp] = mean_read_count/value
            mapped_read_norm_factor[exp] = .9
    return mapped_read_norm_factor

@timing
def load_raw_files(file_names):
    """This will load raw .bam files and affymetrix .CEL files into the associated
    genome_data tables.  The raw signal will go into a mongoDB genome_data table and
    the rest will go into corresponding SQL tables
    """
    mapped_read_norm_factor = calculate_normalization_factors()
    for file_name in file_names:
        if file_name[-3:] != 'bam' and file_name[-3:] != 'CEL': continue
        try:
            norm_factor = mapped_read_norm_factor[file_name]
            data_loading.name_based_experiment_loading(file_name, bulk_file_load=True, norm_factor=norm_factor)
            data_loading.name_based_experiment_loading(file_name, bulk_file_load=True, norm_factor=1.)
        except:
            data_loading.name_based_experiment_loading(file_name, bulk_file_load=True)


def query_experiment_sets():
    """This queries all the replicates of an experiment and groups them for further processing"""
    experiment_sets = {}

    ome = base.Session()

    experiment_sets['RNAseq'] = ome.query(func.array_agg(RNASeqExperiment.name),func.array_agg(RNASeqExperiment.file_name)).\
                                          filter(RNASeqExperiment.normalization_factor == 1.).\
                                          group_by(RNASeqExperiment.strain_id, RNASeqExperiment.environment_id,\
                                          RNASeqExperiment.machine_id,RNASeqExperiment.sequencing_type).all()

    experiment_sets['array'] = ome.query(func.array_agg(ArrayExperiment.name),func.array_agg(ArrayExperiment.file_name)).\
                                         group_by(ArrayExperiment.strain_id, ArrayExperiment.environment_id,\
                                         ArrayExperiment.platform).all()

    experiment_sets['ChIP'] = ome.query(func.array_agg(ChIPExperiment.name),func.array_agg(ChIPExperiment.file_name)).\
                                        filter(ChIPExperiment.normalization_factor == 1.).\
                                        group_by(ChIPExperiment.strain_id, ChIPExperiment.environment_id,\
                                        ChIPExperiment.antibody, ChIPExperiment.protocol_type,\
                                        ChIPExperiment.target).all()
    ome.close()
    return experiment_sets


def load_experiment_sets(experiment_sets):
    """This will create the database entries for the grouped experiments"""

    ome = base.Session()

    for exp_group in experiment_sets['RNAseq']:
        exp_group_name = '_'.join(exp_group[0][0].split('_')[:-2])
        exp = ome.query(RNASeqExperiment).filter_by(name=exp_group[0][0]).one()
        exp_analysis = ome.get_or_create(NormalizedExpression, name=exp_group_name, environment_id=exp.environment.id,\
                                         strain_id=exp.strain.id)
        for exp_name in exp_group[0]:
            exp = ome.query(RNASeqExperiment).filter_by(name=exp_name).one()
            ome.get_or_create(AnalysisComposition, analysis_id = exp_analysis.id, data_set_id = exp.id)

    for exp_group in experiment_sets['array']:
        exp_group_name = '_'.join(exp_group[0][0].split('_')[:-2])
        exp = ome.query(ArrayExperiment).filter_by(name=exp_group[0][0]).one()
        exp_analysis = ome.get_or_create(NormalizedExpression, name=exp_group_name, environment_id=exp.environment.id,\
                                         strain_id=exp.strain.id)
        for exp_name in exp_group[0]:
            exp = ome.query(ArrayExperiment).filter_by(name=exp_name).one()
            ome.get_or_create(AnalysisComposition, analysis_id = exp_analysis.id, data_set_id = exp.id)

    default_parameters = {'mrc':20, 'smooth':3, 'nrf':'', 'outNP':'', 'nf':'', 'k_min': 4, 'k_max': 22}

    for exp_group in experiment_sets['ChIP']:
        parameters = {'mrc':20, 'smooth':3, 'outNP':'', 'nrf':'', 'nf':'','k_min': 4, 'k_max': 22}
        if not set(parameters.items()) - set(default_parameters.items()):
            parameter_name = 'default'
        else:
            parameter_name = '_'.join([y+'-'+str(z) for y,z in dict(set(parameters.items()) - set(default_parameters.items())).iteritems()])


        exp_group_name = '_'.join(exp_group[0][0].split('_')[0:5]+exp_group[0][0].split('_')[6:7]+[parameter_name,'peaks'])
        exp = ome.query(ChIPExperiment).filter_by(name=exp_group[0][0]).one()

        exp_analysis = ome.get_or_create(ChIPPeakAnalysis, name=exp_group_name, environment_id=exp.environment.id, strain_id=exp.strain.id,\
                                     parameters=json.dumps(parameters))
        for exp_name in exp_group[0]:
            if exp_name.split('_')[-1] != '1.0': continue
            exp = ome.query(ChIPExperiment).filter_by(name=exp_name).one()
            ome.get_or_create(AnalysisComposition, analysis_id = exp_analysis.id, data_set_id = exp.id)

    ome.close()


def run_parallel_cuffquant():
    from IPython.parallel import Client

    c = Client()
    c[:].execute('from PrototypeDB.loading import load_data')
    v = c[:]
    d = v.map(lambda x:load_data.run_cuffquant(x), [exp for exp in ome.query(RNASeqExperiment).all()])


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")


if __name__ == "__main__":

    #if not query_yes_no('This will drop the ENTIRE database and load from scratch, ' + \
    #                    'are you sure you want to do this?'): sys.exit()
    """
    base.Base.metadata.drop_all()
    base.omics_database.genome_data.drop()
    base.Base.metadata.create_all()


    file_names = os.listdir(settings.dropbox_directory+'/crp/data/ChIP/bam') + \
                 os.listdir(settings.dropbox_directory+'/crp/data/RNAseq/bam') + \
                 os.listdir(settings.dropbox_directory+'/om_data/Microarray/asv2') + \
                 os.listdir(settings.dropbox_directory+'/om_data/Microarray/ec2')

    load_raw_files(file_names)

    experiment_sets = query_experiment_sets()

    load_experiment_sets(experiment_sets)
    """
    component_loading.load_genes(base, components)
    #component_loading.load_proteins(base, components)
    #component_loading.load_bindsites(base, components)
    #component_loading.load_transcription_units(base, components)

    session = base.Session()

    #data_loading.run_parallel_cuffquant()
    #data_loading.run_cuffnorm(experiment_sets['RNAseq'])
    """query all RNASeqExperiments grouped across replicates"""
    rna_seq_exp_sets = session.query(NormalizedExpression).join(AnalysisComposition, NormalizedExpression.id == AnalysisComposition.analysis_id).\
                                                           join(RNASeqExperiment, RNASeqExperiment.id == AnalysisComposition.data_set_id).all()
    #data_loading.run_cuffdiff(rna_seq_exp_sets, debug=False)
    #data_loading.run_gem(session.query(ChIPPeakAnalysis).filter_by(id=335).all(),debug=False)



    #data_loading.load_gem(session.query(ChIPPeakAnalysis).all())
    #data_loading.load_cuffnorm()
    data_loading.load_cuffdiff()
    #data_loading.load_arraydata(settings.dropbox_directory+'/om_data/Microarray/formatted_asv2.txt', type='asv2')
    #data_loading.load_arraydata(settings.dropbox_directory+'/om_data/Microarray/formatted_ec2.txt', type='ec2')


    genome_data = base.omics_database.genome_data
    genome_data.create_index([("data_set_id",ASCENDING), ("leftpos", ASCENDING)])
    session.close()

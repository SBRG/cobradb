from PrototypeDB.orm import base
from PrototypeDB.orm import data
from PrototypeDB.orm.data import *
from PrototypeDB.orm import components
from PrototypeDB.orm.components import *
from PrototypeDB.loading.load_data import *
from PrototypeDB.lib import settings
from sqlalchemy import func
import simplejson as json
import numpy as np
import os
import math
import shutil
import subprocess


ome = base.Session()


files_to_load = os.listdir(settings.dropbox_directory+'/crp/data/ChIP/bam') + \
                os.listdir(settings.dropbox_directory+'/crp/data/RNAseq/bam') + \
                os.listdir(settings.dropbox_directory+'/ome/data/Microarray/asv2') + \
                os.listdir(settings.dropbox_directory+'/ome/data/Microarray/ec2')

for file_name in files_to_load:
    if file_name[-3:] != 'bam' and file_name[-3:] != 'CEL': continue
    name_based_experiment_loading(file_name, bulk_file_load=True)
    print file_name

#genome_data.create_index([("data_set_id",ASCENDING), ("leftpos", ASCENDING)])

rnaseq_exp_sets = ome.query(func.array_agg(RNASeqExperiment.name),func.array_agg(RNASeqExperiment.file_name)).\
                            group_by(RNASeqExperiment.strain_id, RNASeqExperiment.environment_id,\
                                     RNASeqExperiment.machine_ID,RNASeqExperiment.sequencing_type).all()

array_exp_sets = ome.query(func.array_agg(ArrayExperiment.name),func.array_agg(ArrayExperiment.file_name)).\
                            group_by(ArrayExperiment.strain_id, ArrayExperiment.environment_id,\
                                     ArrayExperiment.platform).all()
    
chip_exp_sets = ome.query(func.array_agg(ChIPExperiment.name),func.array_agg(ChIPExperiment.file_name)).\
                          group_by(ChIPExperiment.strain_id, ChIPExperiment.environment_id,\
                                   ChIPExperiment.antibody, ChIPExperiment.protocol_type,\
                                   ChIPExperiment.target).all()

                                                 
                               
for exp_group in rnaseq_exp_sets: 
    exp_group_name = '_'.join(exp_group[0][0].split('_')[:-1])
    exp = ome.query(RNASeqExperiment).filter_by(name=exp_group[0][0]).one()
    exp_analysis = ome.get_or_create(NormalizedExpression, name=exp_group_name, environment=exp.environment, strain=exp.strain)
    for exp_name in exp_group[0]:
        exp = ome.query(RNASeqExperiment).filter_by(name=exp_name).one()
        ome.get_or_create(AnalysisComposition, analysis_id = exp_analysis.id, data_set_id = exp.id)

for exp_group in array_exp_sets: 
    exp_group_name = '_'.join(exp_group[0][0].split('_')[:-1])
    exp = ome.query(ArrayExperiment).filter_by(name=exp_group[0][0]).one()
    exp_analysis = ome.get_or_create(NormalizedExpression, name=exp_group_name, environment=exp.environment, strain=exp.strain)
    for exp_name in exp_group[0]:
        exp = ome.query(ArrayExperiment).filter_by(name=exp_name).one()
        ome.get_or_create(AnalysisComposition, analysis_id = exp_analysis.id, data_set_id = exp.id)
    
for exp_group in chip_exp_sets:
    exp_group_name = '_'.join(exp_group[0][0].split('_')[0:5]+exp_group[0][0].split('_')[6:7]) 
    exp = ome.query(ChIPExperiment).filter_by(name=exp_group[0][0]).one()
    exp_analysis = ome.get_or_create(ChIPPeak, name=exp_group_name, environment=exp.environment, strain=exp.strain)
    for exp_name in exp_group[0]:
        exp = ome.query(ChIPExperiment).filter_by(name=exp_name).one()
        ome.get_or_create(AnalysisComposition, analysis_id = exp_analysis.id, data_set_id = exp.id)

#load_gem()
load_cuffnorm()
#load_cuffdiff()
load_arraydata(settings.dropbox_directory+'/ome/data/Microarray/formatted_asv2.txt')
load_arraydata(settings.dropbox_directory+'/ome/data/Microarray/formatted_ec2.txt')
    
    
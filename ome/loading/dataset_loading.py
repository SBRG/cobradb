####Many parts of this code are derived from sequtil written by aebrahim####
#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
from ome import base, datasets, components, settings, timing
from os.path import split
from math import log
from itertools import combinations

import os,subprocess,sys,math

from scipy.stats import ttest_ind
from numpy import zeros, roll, array, mean, genfromtxt, zeros_like, arange
from sqlalchemy import func, or_, and_
from IPython import embed

import pysam
import simplejson as json


def count_coverage(samfile, flip=False, include_insert=False):
    """counts coverage per base in a strand-specific manner

    include_insert: If the insert between paired end reads should be
        included in the counts.

    flip: Whether or not the strands should be flipped.
    This should be true for RNA-seq, and false for ChIP-exo
"""

    all_counts = {}
    plus_strands = []
    minus_strands = []
    chromosome_sizes = dict(zip(samfile.references, samfile.lengths))
    for i in chromosome_sizes:
        chromosome_sizes[i] += 2  # allows to roll later, extra 0's never hurt
    for reference in samfile.references:  # create an array for each reference
        plus_strands.append(zeros((chromosome_sizes[reference],)))
        minus_strands.append(zeros((chromosome_sizes[reference],)))
    # iterate through each mapped read
    for i, read in enumerate(samfile):
        #if i > 1e3: break
        if read.is_unmapped:
            continue
        # for paired and data get entire insert only from read 1
        if include_insert and read.is_proper_pair:
            if read.is_read2:
                continue  # will get handled with read 1
            if read.is_reverse:
                minus_strands[read.tid][read.pnext:read.aend] += 1
            else:
                plus_strands[read.tid][read.pos:read.pos + read.isize] += 1
        else:
            # Truth table for where reads are mapped
            # read2 is flipped

            # is_read1  is_reverse      outcome
            # ---------------------------------
            # True      False           +
            # True      True            -
            # False     False           -
            # False     True            +

            # therefore read1 == is_reverse --> negative
            #           read1 != is_reverse --> positive
            if read.is_reverse == read.is_read1:
                minus_strands[read.tid][read.pos:read.aend] += 1
            else:
                plus_strands[read.tid][read.pos:read.aend] += 1
    # store the results per reference
    for i, reference in enumerate(samfile.references):
        all_counts[reference] = {}
        # roll shifts by 1, so the first base position (at index 0) is now at
        # index 1
        print reference
        if flip:
            all_counts[reference]["-"] = roll(plus_strands[i], 1)
            all_counts[reference]["+"] = roll(minus_strands[i], 1)
        else:
            all_counts[reference]["+"] = roll(plus_strands[i], 1)
            all_counts[reference]["-"] = roll(minus_strands[i], 1)
    return all_counts


def count_coverage_5prime(samfile, flip=False):
    """counts the coverage of 5' ends per base in a strand-specific manner

    flip: Whether or not the strands should be flipped.
    This should be true for RNA-seq, and false for ChIP-exo
"""

    all_counts = {}
    plus_strands = []
    minus_strands = []
    chromosome_sizes = dict(zip(samfile.references, samfile.lengths))
    for i in chromosome_sizes:
        chromosome_sizes[i] += 2  # allows to roll later, extra 0's never hurt
    for reference in samfile.references:  # create an array for each reference
        plus_strands.append(zeros((chromosome_sizes[reference],)))
        minus_strands.append(zeros((chromosome_sizes[reference],)))
    # iterate through each mapped read
    for i, read in enumerate(samfile):
        #if i > 1e5: break
        if read.is_unmapped:
            continue
        if read.is_reverse:
            minus_strands[read.tid][read.aend - 1] += 1
        else:
            plus_strands[read.tid][read.pos] += 1
    # store the results per reference
    for i, reference in enumerate(samfile.references):
        all_counts[reference] = {}
        # roll shifts by 1, so the first base position (at index 0) is now at
        # index 1
        if flip:
            all_counts[reference]["-"] = roll(plus_strands[i], 1)
            all_counts[reference]["+"] = roll(minus_strands[i], 1)
        else:
            all_counts[reference]["+"] = roll(plus_strands[i], 1)
            all_counts[reference]["-"] = roll(minus_strands[i], 1)
    return all_counts


def gff_variance(gff_file_path):

    n = 0
    mean = 0
    M2 = 0

    for line in open(gff_file_path, 'r').readlines():
        vals = line.split('\t')

        x = abs(float(vals[5]))
        n = n + 1
        delta = x - mean
        mean = mean + delta/n
        M2 = M2 + delta*(x - mean)

    variance = M2/(n - 1)
    return variance


def write_samfile_to_gff(sam_filename, out_filename, flip=False, log2=False,
        separate_strand=False, include_insert=False, five_prime=False,
        track=None):
    """
    write samfile object to an output object in a gff format

    flip: Whether or not the strands should be flipped.
    This should be true for RNA-seq, and false for ChIP-exo

    separate_strand: Whether the forward and reverse strands should be made
    into separate tracks (True) or the negative strand should be rendered
    as negative values (False)

    log2: Whether intensities should be reported as log2.
    """

    samfile = pysam.Samfile(sam_filename)
    if five_prime:
        all_counts = count_coverage_5prime(samfile, flip=flip)
    else:
        all_counts = count_coverage(samfile,
            include_insert=include_insert, flip=flip)
    if track is None:
        name = split(samfile.filename)[1]
    else:
        name = track
    gff_base = "%s\t\t%s\t%d\t%d\t%s\t%s\t.\t.\n"
    if log2:
        str_func = lambda x, s: "%.2f" % (log(x, 2) * s)
    else:
        str_func = lambda x, s: "%d" % (x * s)
    output = open(out_filename, "w")
    for reference in all_counts:
        for strand in all_counts[reference]:
            factor = 1 if strand == "+" else -1
            track_name = "%s_(%s)" % (name, strand) if separate_strand else name
            counts = all_counts[reference][strand]
            for i in counts.nonzero()[0]:
                output.write(gff_base % (reference, track_name, i, i,
                                         str_func(counts[i], factor), strand))
    output.close()
    samfile.close()


def load_samfile_to_db(sam_filepath, dataset_id, loading_cutoff=0, bulk_file_load=False,
                       flip=False, log2=False, separate_strand=False, include_insert=False,
                       five_prime=False, track=None, norm_factor=1.):
    """
    write samfile object to an output object in a gff format

    flip: Whether or not the strands should be flipped.
    This should be true for RNA-seq, and false for ChIP-exo

    separate_strand: Whether the forward and reverse strands should be made
    into separate tracks (True) or the negative strand should be rendered
    as negative values (False)

    log2: Whether intensities should be reported as log2.
    """

    samfile = pysam.Samfile(sam_filepath)
    if five_prime:
        all_counts = count_coverage_5prime(samfile, flip=flip)
    else:
        all_counts = count_coverage(samfile,
            include_insert=include_insert, flip=flip)
    #connection to mongo data store
    genome_data = base.omics_database.genome_data

    if log2:
        str_func = lambda x, s: "%.2f" % (log(x, 2) * s)
    else:
        str_func = lambda x, s: "%d" % (x * s)

    for reference in all_counts:
        for strand in all_counts[reference]:
            factor = 1 if strand == "+" else -1
            counts = all_counts[reference][strand]
            entries = []
            for i in counts.nonzero()[0]:
                if abs(float(counts[i])) < loading_cutoff: continue

                entries.append({
                                "leftpos": int(i),
                                "rightpos": int(i),
                                "value": float(counts[i])*norm_factor,
                                "strand": strand,
                                "dataset_id": dataset_id})

                if i%50000 == 0:
                    genome_data.insert(entries)
                    entries = []
            genome_data.insert(entries)
    if not bulk_file_load:
        genome_data.create_index([("dataset_id", ASCENDING), ("leftpos", ASCENDING)])

    samfile.close()


@timing
def load_raw_gff_to_db(experiment):
    filepath = settings.data_directory+'/chip_experiment/gff/'+experiment.name+'.gff'
    assert filepath.endswith(".gff")
    genome_data = base.omics_database.genome_data
    entries = []
    with open(filepath) as infile:
        for i,line in enumerate(infile):
            if line[0] == '#': continue
            data = line.split("\t")
            if float(data[5]) < .1: continue

            entries.append({
                "leftpos": int(data[3]),
                "rightpos": int(data[4]),
                "value": float(data[5]),
                "strand": data[6],
                "dataset_id": experiment.id})

            if i%50000 == 0:
                genome_data.insert(entries)
                entries = []
        genome_data.insert(entries)



def calculate_normalization_factors(experiment_type, group_name):
    if experiment_type not in ['rnaseq_experiment','chip_experiment']: return
    if group_name == 'default':
        directory_path = os.path.join(settings.data_directory, experiment_type, 'bam')
    else:
        directory_path = os.path.join(settings.data_directory, experiment_type, 'bam', group_name)

    mapped_read_norm_factor = {}
    mapped_reads = {}
    for file_name in os.listdir(directory_path):
        if file_name[-3:] != 'bam': continue
        prefix = file_name.split('.')[0]

        mapped_reads[prefix] = int(subprocess.check_output(['samtools', 'flagstat', directory_path+'/'+file_name]).split('\n')[2].split()[0])

        #mapped_reads[prefix] = 3e6
    mean_read_count = array(mapped_reads.values()).mean()
    for exp,value in mapped_reads.iteritems():
        mapped_read_norm_factor[exp] = mean_read_count/value
    return mapped_read_norm_factor


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

        experiment = create_name_based_experiment(session, file_prefix, group_name)

        experiments.append(experiment)


    for experiment in set(experiments):
      if experiment.name in fastq_experiment_paths.keys():
          run_bowtie2(experiment, fastq_experiment_paths[experiment.name], overwrite=overwrite, debug=True)


    if normalize:
        normalization_factors = calculate_normalization_factors(experiments[0].type, group_name)
    else:
        normalization_factors = {exp.name: 1. for exp in experiments}


    for experiment in set(experiments):
        if not raw: continue

        if experiment.type == 'chip_experiment':
            norm_factor = normalization_factors[experiment.name]
            if experiment.protocol_type == 'ChIPExo':
                load_raw_experiment_data(experiment, loading_cutoff=10, flip=False, five_prime=True, norm_factor=norm_factor)
            elif experiment.protocol_type == 'ChIPchip':
                load_raw_gff_to_db(experiment)

        elif experiment.type == 'rnaseq_experiment':
            norm_factor = normalization_factors[experiment.name]
            load_raw_experiment_data(experiment, loading_cutoff=10., flip=True, five_prime=False, norm_factor=norm_factor)

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

    RNASeqExperiment = datasets.RNASeqExperiment
    ArrayExperiment = datasets.ArrayExperiment
    ChIPExperiment = datasets.ChIPExperiment

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

        exp = ome.query(datasets.RNASeqExperiment).filter_by(name=exp_group[0][0], group_name=exp_group[1]).one()

        exp_analysis = ome.get_or_create(datasets.NormalizedExpression, replicate=1, name=exp_group_name, environment_id=exp.environment.id,\
                                         strain_id=exp.strain.id, group_name=exp.group_name, expression_type='rnaseq_experiment')

        for exp_name in exp_group[0]:
            expt = ome.query(datasets.RNASeqExperiment).filter_by(name=exp_name, group_name=exp.group_name).one()
            ome.get_or_create(datasets.AnalysisComposition, analysis_id = exp_analysis.id, dataset_id = expt.id)

    for exp_group in experiment_sets['array']:
        exp_group_name = '_'.join(exp_group[0][0].split('_')[:-2])
        exp = ome.query(datasets.ArrayExperiment).filter_by(name=exp_group[0][0]).one()
        exp_analysis = ome.get_or_create(datasets.NormalizedExpression, replicate=1, name=exp_group_name, environment_id=exp.environment.id,\
                                         strain_id=exp.strain.id, group_name=exp.group_name, expression_type='array_experiment')
        for exp_name in exp_group[0]:
            exp = ome.query(datasets.ArrayExperiment).filter_by(name=exp_name).one()
            ome.get_or_create(datasets.AnalysisComposition, analysis_id = exp_analysis.id, dataset_id = exp.id)

    default_parameters = {'mrc':20, 'smooth':3, 'nrf':'', 'outNP':'', 'nf':'', 'k_min': 4, 'k_max': 22, 'k_win':150}

    for exp_group in experiment_sets['ChIP']:
        parameters = {'mrc':20, 'smooth':3, 'outNP':'', 'nrf':'', 'nf':'','k_min': 4, 'k_max': 22, 'k_win':150}
        if not set(parameters.items()) - set(default_parameters.items()):
            parameter_name = 'default'
        else:
            parameter_name = '_'.join([y+'-'+str(z) for y,z in dict(set(parameters.items()) - set(default_parameters.items())).iteritems()])


        vals = exp_group[0][0].split('_')
        exp_group_name = '_'.join(vals[0:5]+vals[6:]+[parameter_name,'peaks'])


        exp = ome.query(datasets.ChIPExperiment).filter_by(name=exp_group[0][0]).one()


        exp_analysis = ome.get_or_create(datasets.ChIPPeakAnalysis, name=exp_group_name, environment_id=exp.environment.id, strain_id=exp.strain.id,\
                                                                parameters=json.dumps(parameters), replicate=1, group_name=exp.group_name)
        for exp_name in exp_group[0]:
            expt = ome.query(datasets.ChIPExperiment).filter_by(name=exp_name, group_name=exp.group_name).one()
            ome.get_or_create(datasets.AnalysisComposition, analysis_id = exp_analysis.id, dataset_id = expt.id)

    ome.close()




def create_name_based_experiment(session, exp_name, group_name, lab='palsson', institution='UCSD'):
    vals = exp_name.split('_')

    if len(vals) < 6: return  #required to have experiment-type_strain_carbon-source_nitrogen-source_electron-acceptor_replicate

    exp_type = vals[0].split('-')

    try: supplements = vals[7]
    except:
        try:
            if vals[6][0:4] == 'anti' or vals[6] in ['asv2','ec2']:
                supplements = ''
            else:
                supplements = vals[6]
        except: supplements = ''

    strain = session.get_or_create(datasets.Strain, name=vals[1])

    data_source = session.get_or_create(base.DataSource, name=vals[0], lab=lab, institution=institution)

    environment = session.get_or_create(datasets.InVivoEnvironment, name='_'.join(vals[2:5]+[supplements]), carbon_source=vals[2],\
                                        nitrogen_source=vals[3], electron_acceptor=vals[4], temperature=37,\
                                        supplements=supplements)


    if exp_type[0][0:4] == 'ChIP':

        experiment = session.get_or_create(datasets.ChIPExperiment, name=exp_name, replicate=vals[5],\
                                           strain_id=strain.id, data_source_id=data_source.id, environment_id=environment.id,\
                                           protocol_type=exp_type[0], antibody=vals[6], target=exp_type[1], group_name=group_name)


    elif exp_type[0][0:6] == 'RNAseq':
        experiment = session.get_or_create(datasets.RNASeqExperiment, name=exp_name, replicate=vals[5],\
                                           strain_id=strain.id, data_source_id=data_source.id, environment_id=environment.id,\
                                           machine_id='miseq', sequencing_type='unpaired', group_name=group_name)


    elif exp_type[0][0:7] == 'affyexp':

        experiment = session.get_or_create(datasets.ArrayExperiment, name=exp_name, replicate=vals[5],\
                                           strain_id=strain.id, data_source_id=data_source.id, environment_id=environment.id,\
                                           platform=vals[6], group_name=group_name)

    session.commit()
    return experiment


@timing
def load_raw_experiment_data(experiment, loading_cutoff=5., bulk_file_load=True, five_prime=False, flip=True, norm_factor=1.):

    if experiment.group_name == 'default':
        file_path = os.path.join(settings.data_directory, experiment.type, 'bam', experiment.name+'.bam')
    else:
        file_path = os.path.join(settings.data_directory, experiment.type, 'bam', experiment.group_name, experiment.name+'.bam')


    load_samfile_to_db(file_path, experiment.id, loading_cutoff=loading_cutoff,\
                       bulk_file_load=bulk_file_load, five_prime=five_prime, flip=flip,
                       norm_factor=norm_factor)


def run_bowtie2(experiment, fastq_paths, overwrite=False, debug=False):
    """This function runs bowtie2 on the raw fastq files associated with an experiment object.

       It first checks to see if a bam file exists and by default prompts the user on whether or
       not to overwrite an existing bam file.
    """

    if experiment.group_name == 'default':
        bam_file_path = settings.data_directory+'/'+experiment.type+'/bam/'+experiment.name
        fastq_dir_path = settings.data_directory+'/'+experiment.type+'/fastq/'
    else:
        bam_file_path = settings.data_directory+'/'+experiment.type+'/bam/'+experiment.group_name+'/'+experiment.name
        fastq_dir_path = settings.data_directory+'/'+experiment.type+'/fastq/'+experiment.group_name+'/'

    if os.path.isfile(bam_file_path+'.bam') and overwrite:
        os.remove(bam_file_path+'.bam')
    elif os.path.isfile(bam_file_path+'.bam'): return

    """Check for proper bowtie2 index and if not build it"""
    genbank_dir = settings.data_directory+'/annotation/genbank/'
    os.chdir(genbank_dir)

    if not os.path.isfile(genbank_dir+'/NC_000913.2.1.bt2'):
        os.system("bowtie2-build %s %s" % (genbank_dir+'/NC_000913.2.fna', 'NC_000913.2'))

    R1_list = []
    R2_list = []
    for path in fastq_paths:
        suffix = '.'.join(path.split('.')[1:])
        if suffix[0:2] == 'R1': R1_list.append(path)
        elif suffix[0:2] == 'R2': R2_list.append(path)

    #print fastq_paths
    if len(R1_list+R2_list) == 0: #unpaired
        bowtie_string = "bowtie2 -N %d -p %d -x %s -U %s | samtools view -bS - | samtools sort - %s" % \
                                   (1, 8, 'NC_000913.2', ','.join([fastq_dir_path+x for x in fastq_paths]), bam_file_path)
    else: #paired
        bowtie_string =  "bowtie2 -N %d -p %d -x %s -1 %s -2 %s | samtools view -bS - | samtools sort - %s" % \
                                   (1, 8, 'NC_000913.2', ','.join([fastq_dir_path+x for x in R1_list]),
                                                         ','.join([fastq_dir_path+x for x in R2_list]), bam_file_path)

    if debug:
      print bowtie_string
    else:
      print bowtie_string
      os.system(bowtie_string)


@timing
def run_cuffquant(base, datasets, genome, group_name=None, overwrite=False, debug=False):


    gff_file = settings.data_directory+'/annotation/'+genome.ncbi_id+'_'+group_name+'.gff'
    cxb_dir = settings.data_directory+'/rnaseq_experiment/cxb/'+group_name

    if not os.path.exists(cxb_dir):
        os.mkdir(cxb_dir)

    session = base.Session()

    for experiment in session.query(datasets.RNASeqExperiment).filter(datasets.RNASeqExperiment.group_name == group_name).all():

        out_path = cxb_dir+'/'+experiment.name

        exp_file = settings.data_directory+'/rnaseq_experiment/bam/'+experiment.group_name+'/'+experiment.name+'.bam'

        cuffquant_string = '%s -p %d -v --library-type fr-firststrand %s %s' % (settings.cufflinks+'/cuffquant', 8, gff_file, exp_file)

        if debug:
            print cuffquant_string
            continue

        if os.path.exists(out_path+'/abundances.cxb') and overwrite:
            os.system('rm -r '+out_path)
            os.mkdir(out_path)
            os.chdir(out_path)
            os.system(cuffquant_string)

        elif not os.path.exists(out_path+'/abundances.cxb'):
            if not os.path.exists(out_path): os.mkdir(out_path)
            os.chdir(out_path)
            os.system(cuffquant_string)





@timing
def run_cuffnorm(base, datasets, genome, group_name, gff_file=None, debug=False, overwrite=False):

    if not gff_file:
        gff_file = settings.data_directory+'/annotation/'+genome.ncbi_id+'.gff'

    cxb_dir = settings.data_directory+'/rnaseq_experiment/cxb/'+group_name
    out_path = settings.data_directory+'/rnaseq_experiment/cuffnorm/'+group_name

    session = base.Session()
    experiments = session.query(func.array_agg(datasets.RNASeqExperiment.name)).\
                                      filter(datasets.RNASeqExperiment.group_name == group_name).\
                                      group_by(datasets.RNASeqExperiment.strain_id, datasets.RNASeqExperiment.environment_id,\
                                               datasets.RNASeqExperiment.machine_id, datasets.RNASeqExperiment.sequencing_type).all()


    cuffnorm_string = '%s -p %d --library-type fr-firststrand -L %s %s %s' % \
                                    (settings.cufflinks+'/cuffnorm', 24,
                                     ','.join([exps[0][0] for exps in experiments]),
                                     gff_file,
                                     ' '.join([','.join([cxb_dir+'/'+exp_name+'/abundances.cxb' for exp_name in exps[0]]) for exps in experiments]))

    session.close()

    if debug:
        print cuffnorm_string
        return

    if os.path.exists(out_path) and overwrite:
        os.system('rm -r '+out_path)
        os.mkdir(out_path)
        os.chdir(out_path)
        os.system(cuffnorm_string)

    elif not os.path.exists(out_path):
        os.mkdir(out_path)
        os.chdir(out_path)
        os.system(cuffnorm_string)




def find_single_factor_pairwise_contrasts(dataset_list):
    dataset_contrasts = []
    dataset_conditions = {}
    for dataset in dataset_list:
        dataset_conditions[(dataset.strain, dataset.environment)] = dataset
    for c in combinations(dataset_conditions.keys(), 2):
        e1, e2 = sorted(c, key=str)
        s1, c1 = e1
        s2, c2 = e2
        if s1 == s2:  # if the strains are the same
           # single shift only
            differences = (c1.carbon_source != c2.carbon_source) + \
                          (c1.nitrogen_source != c2.nitrogen_source) + \
                          (c1.electron_acceptor != c2.electron_acceptor) + \
                          (c1.supplements != c2.supplements)
            if differences != 1:
                continue
        else:  # if the strains are different
            if not (s1.name == "wt" or s2.name == "wt"):
                continue
            if c1 != c2:  # make sure the conditions are the same
                continue
        exp_1 = dataset_conditions[(s1, c1)]
        exp_2 = dataset_conditions[(s2, c2)]
        dataset_contrasts.append([exp_1,exp_2])

    return dataset_contrasts


def generate_cuffdiff_contrasts(normalized_expression_objects, group_name):
    contrasts = find_single_factor_pairwise_contrasts(normalized_expression_objects)
    with open(settings.data_directory+'/rnaseq_experiment/cuffdiff/'+group_name+'/contrasts.txt', 'wb') as contrast_file:
        contrast_file.write('condition_A\tcondition_B\n')
        for contrast in contrasts:
            contrast_file.write(contrast[0].name+'\t'+contrast[1].name+'\n')


@timing
def run_cuffdiff(base, datasets, genome, group_name, gff_file=None, debug=False, overwrite=False):

    if not gff_file:
        gff_file = settings.data_directory+'/annotation/'+genome.ncbi_id+'.gff'

    cxb_dir = settings.data_directory+'/rnaseq_experiment/cxb/'+group_name
    out_path = settings.data_directory+'/rnaseq_experiment/cuffdiff/'+group_name


    session = base.Session()
    exp_objects = session.query(datasets.NormalizedExpression).\
                                   join(datasets.AnalysisComposition, datasets.NormalizedExpression.id == datasets.AnalysisComposition.analysis_id).\
                                   join(datasets.RNASeqExperiment, datasets.RNASeqExperiment.id == datasets.AnalysisComposition.dataset_id).\
                                   filter(datasets.RNASeqExperiment.group_name == group_name).all()


    cuffdiff_string = '%s -v -p %d --library-type fr-firststrand --FDR 0.05 -C %s -L %s %s %s' % \
                          (settings.cufflinks+'/cuffdiff', 24, out_path+'/contrasts.txt',
                           ','.join([x.name for x in exp_objects]), gff_file,
                           ' '.join([','.join([cxb_dir+'/'+x.name+'/abundances.cxb' for x in exp.children]) for exp in exp_objects]))

    if debug:
        print cuffdiff_string
        return

    if os.path.exists(out_path) and overwrite:
        os.system('rm -r '+out_path)
        os.mkdir(out_path)
        os.chdir(out_path)
        generate_cuffdiff_contrasts(exp_objects, group_name)
        os.system(cuffdiff_string)

    elif not os.path.exists(out_path):
        os.mkdir(out_path)
        os.chdir(out_path)
        generate_cuffdiff_contrasts(exp_objects, group_name)
        os.system(cuffdiff_string)



def calculate_differential_expression(base, datasets, experiment1, experiment2):
    """calculate differential expression (fold change and q)"""
    session = base.Session()
    platform = experiment1.children[0].platform

    genes = array([g[0] for g in session.query(func.distinct(datasets.GenomeData.genome_region_id)).\
                                         filter(datasets.GenomeData.dataset_id.in_([x.id for x in experiment1.children+experiment2.children])).all()])


    n_genes = len(genes)
    # query data
    data_vals1 = session.query(datasets.GenomeData.value).filter(datasets.GenomeData.dataset_id.in_([x.id for x in experiment1.children])).\
                                                        order_by(datasets.GenomeData.genome_region_id).all()

    data_vals2 = session.query(datasets.GenomeData.value).filter(datasets.GenomeData.dataset_id.in_([x.id for x in experiment2.children])).\
                                                        order_by(datasets.GenomeData.genome_region_id).all()

    if n_genes == 0 or len(data_vals1) % n_genes != 0 or len(data_vals1) % n_genes != 0 or len(data_vals1) == 0 or len(data_vals2) == 0:
        return array(0),array(0),array(0)

    d1 = array(data_vals1).reshape(n_genes, -1)
    d2 = array(data_vals2).reshape(n_genes, -1)

    # filter the data, using the average of all values in the platform file as a cutoff
    cutoff = 3.5 #mean(genfromtxt("%s/microarray/%s/IG_formatted_%s.tab" % (settings.data_directory, platform, platform))[:, 1:])
    selection = (d1.max(axis=1) > cutoff) + (d2.max(axis=1) > cutoff)  # at least one of the samples must have one value over the cutoff
    d1 = d1[selection, :]
    d2 = d2[selection, :]
    genes = genes[selection]
    # calculate fold change (difference of means)
    fold_change = (d1.mean(axis=1) - d2.mean(axis=1))
    # calculate t-test statistics
    t, p = ttest_ind(d1, d2, axis=1)  # independent t-test
    # perform Benjamini-Hochberg correction (similar to p.adjust(p_values, method="BH") in R)
    n_total = len(p)
    ranks = zeros_like(p)
    ranks[p.argsort()] = arange(n_total) + 1.0  # ranks must be floats starting with 1
    q = p * n_total / ranks  # each entry is scaled by n_total / it's rank
    q[q > 1] = 1.0  # maximum value is 1
    if not q[0] >= 0:
        from IPython import embed; embed()
    session.close()
    return genes, fold_change, q




@timing
def run_array_ttests(base, datasets, genome, group_name, debug=False, overwrite=False):
    session = base.Session()
    exp_objects = session.query(datasets.NormalizedExpression).filter(datasets.NormalizedExpression.group_name == group_name).\
                                   join(datasets.AnalysisComposition, datasets.NormalizedExpression.id == data.AnalysisComposition.analysis_id).\
                                   join(datasets.ArrayExperiment, data.ArrayExperiment.id == datasets.AnalysisComposition.dataset_id).all()

    contrasts = find_single_factor_pairwise_contrasts(exp_objects)

    for experiment1,experiment2 in contrasts:

        x = experiment1.name.split('_')
        y = experiment2.name.split('_')
        exp_name = ''
        if len(x) != len(y):
            if len(x) > len(y):
                exp_name = '_'.join(y)+'_supp/'+x[-1]
            else:
                exp_name = '_'.join(x)+'_supp/'+y[-1]
        else:
            for i in range(len(x)):
                if x[i] == y[i]:
                    exp_name += x[i]+'_'
                else: exp_name += x[i]+'/'+y[i]+'_'
            exp_name = exp_name.rstrip('_')

        #print experiment1,experiment2
        diff_exp = session.get_or_create(datasets.DifferentialExpression, name=exp_name, replicate=1, norm_method='gcrma',fdr=.05, group_name=group_name)

        session.get_or_create(datasets.AnalysisComposition, analysis_id = diff_exp.id, dataset_id = experiment1.id)
        session.get_or_create(datasets.AnalysisComposition, analysis_id = diff_exp.id, dataset_id = experiment2.id)


        genes, fold_change, q = calculate_differential_expression(base, datasets, experiment1, experiment2)

        if genes.all() == 0: continue

        for i,gene_id in enumerate(genes):

            if math.isnan(fold_change[i]): value = 0.
            else: value = fold_change[i]

            if math.isnan(q[i]): pval = 0.
            else: pval = q[i]

            gene_data = datasets.DiffExpData(dataset_id = diff_exp.id,\
                                             genome_region_id = gene_id,\
                                             value=value, pval=pval)
            session.add(gene_data)

    session.flush()
    session.commit()
    session.close()


@timing
def run_gem(base, datasets, genome, debug=False, overwrite=False, with_control=False):
    default_parameters = {'mrc':20, 'smooth':3, 'nrf':'', 'outNP':''}
    gem_path = settings.home_directory+'/libraries/gem'
    bam_dir = settings.data_directory+'/chip_experiment/bam/'

    session = base.Session()

    for chip_peak_analysis in session.query(datasets.ChIPPeakAnalysis).all():

        outdir = chip_peak_analysis.name
        out_path = settings.data_directory+'/chip_peaks/gem/'+outdir


        input_files = ' '.join(['--expt'+x.name+' '+bam_dir+x.group_name+'/'+x.name+'.bam' for x in chip_peak_analysis.children])


        #control_peak_analysis = session.query(ChIPPeakAnalysis).join(AnalysisComposition, ChIPPeakAnalysis.id == AnalysisComposition.analysis_id).\
        #                                                        join(ChIPExperiment, ChIPExperiment.id == AnalysisComposition.dataset_id).\
        #                                                        join(Strain).\
        #                                                        filter(and_(Strain.name == 'delta-crp',
        #                                                                    ChIPExperiment.antibody == 'anti-crp')).one()

        if with_control:
            from random import randint
            """since every replicate needs a corresponding control replicate for GEMs and since we might have less or more
               control replicates than experimental replicates, we are going to randomly sample from the set of control
               replicates for each experimental replicate.  This should be fine because control replicates should be random
               noise anyways. However, if specific control replicates are generated this function needs to be changed.
            """
            control_input_files = ' '.join(['--ctrl'+x.name+' '+bam_dir+control_chip_peak_analysis.children[randint(0,len(control_chip_peak_analysis.children))].file_name for i,x in enumerate(chip_peak_analysis.children)])
        else:
            control_input_files = ''


        params = json.loads(chip_peak_analysis.parameters)
        parameter_string = ' '.join(['--'+y+' '+str(z) for y,z in params.iteritems()])


        gem_string = "java -Xmx5G -jar %s/gem.jar --d %s/Read_Distribution_ChIP-exo.txt --g %s --genome %s %s %s --f SAM %s" %\
                     (gem_path, gem_path, settings.data_directory+'/annotation/ec_mg1655.sizes', settings.data_directory+'/annotation',\
                     input_files, control_input_files, parameter_string)

        if debug:
            print gem_string
            continue


        if os.path.exists(out_path) and overwrite:
            os.system('rm -r '+out_path)
            os.mkdir(out_path)
            os.chdir(out_path)
            os.system(gem_string)

        elif not os.path.exists(out_path):
            os.mkdir(out_path)
            os.chdir(out_path)
            os.system(gem_string)

            gem_peak_file = open(out_path+'/out_GPS_events.narrowPeak','r')
            with open(out_path+'/'+chip_peak_analysis.name+'_gps.gff', 'wb') as peaks_gff_file:

                for line in gem_peak_file.readlines():
                    vals = line.split('\t')

                    position = int(vals[3].split(':')[1])

                    peaks_gff_file.write('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n' %\
                                         ('NC_000913','.', chip_peak_analysis.name+'_gps', position-1, position+1, float(vals[6])*3.2, '+', '.','.'))
                    peaks_gff_file.write('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n' %\
                                         ('NC_000913','.', chip_peak_analysis.name+'_gps', vals[1], vals[2], float(vals[6]), '+', '.','.'))

            try:
                gem_peak_file = open(out_path+'/out_GEM_events.narrowPeak','r')
                with open(out_path+'/'+chip_peak_analysis.name+'_gem.gff', 'wb') as peaks_gff_file:

                    for line in gem_peak_file.readlines():
                        vals = line.split('\t')

                        position = int(vals[3].split(':')[1])

                        peaks_gff_file.write('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n' %\
                                             ('NC_000913','.', chip_peak_analysis.name+'_gem', position-1, position+1, float(vals[6])*3.2, '+', '.','.'))
                        peaks_gff_file.write('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n' %\
                                             ('NC_000913','.', chip_peak_analysis.name+'_gem', vals[1], vals[2], float(vals[6]), '+', '.','.'))
            except: None


    session.close()


@timing
def load_cuffnorm(base, datasets, group_name):

    session = base.Session()

    cuffnorm_output = open(settings.data_directory+'/rnaseq_experiment/cuffnorm/'+group_name+'/isoforms.fpkm_table','r')
    header = cuffnorm_output.readline().rstrip('\n').split('\t')

    exp_id_map = {}
    for name in header[1:]:

        vals = name.split('_')
        if len(vals) > 7:
            exp_name = '_'.join(vals[0:5]+[str(int(vals[5])+int(vals[7]))]+[vals[6]])
        else:
            exp_name = '_'.join(vals[0:5]+[str(int(vals[5])+int(vals[6]))])

        exp_id_map[name] = session.query(datasets.RNASeqExperiment).filter_by(name=exp_name, group_name=group_name).one().id



    for line in cuffnorm_output.readlines():
        vals = line.split('\t')

        try: gene = session.query(components.Gene).filter_by(locus_id=vals[0]).one()
        except: continue

        for i,val in enumerate(vals[1:]):

            try: value = float(val)
            except: continue

            genome_data = datasets.GenomeData(dataset_id=exp_id_map[header[i+1]],
                                              genome_region_id = gene.id,
                                              value=value)
            session.add(genome_data)

    session.flush()
    session.commit()
    session.close()



@timing
def load_cuffdiff(group_name):
    cuffdiff_output = open(settings.data_directory+'/rnaseq_experiment/cuffdiff/'+group_name+'/gene_exp.diff','r')
    header = cuffdiff_output.readline()
    diff_exps = {}

    session = base.Session()
    for line in cuffdiff_output.readlines():
        vals = line.split('\t')

        try:
            value = float(vals[9])
            pvalue = float(vals[11])
        except: continue

        if pvalue > .25: continue

        try: gene = session.query(components.Gene).filter_by(locus_id=vals[0]).one()
        except: continue

        if str(vals[4:6]) not in diff_exps.keys():
            x = vals[4].split('_')
            y = vals[5].split('_')
            exp_name = ''
            if len(x) != len(y):
                if len(x) > len(y):
                    exp_name = '_'.join(y)+'_supp/'+x[-1]
                else:
                    exp_name = '_'.join(x)+'_supp/'+y[-1]
            else:
                for i in range(len(x)):
                    if x[i] == y[i]:
                        exp_name += x[i]+'_'
                    else: exp_name += x[i]+'/'+y[i]+'_'
                exp_name = exp_name.rstrip('_')

            diff_exp = session.get_or_create(datasets.DifferentialExpression, name=exp_name, replicate=1, norm_method='classic-fpkm',fdr=.05, group_name=group_name)


            exp1 = session.query(datasets.Analysis).filter_by(name=vals[4], group_name=group_name).one()
            exp2 = session.query(datasets.Analysis).filter_by(name=vals[5], group_name=group_name).one()
            session.get_or_create(datasets.AnalysisComposition, analysis_id = diff_exp.id, dataset_id = exp1.id)
            session.get_or_create(datasets.AnalysisComposition, analysis_id = diff_exp.id, dataset_id = exp2.id)
            diff_exps[str(vals[4:6])] = diff_exp.id



        diff_exp_data = datasets.DiffExpData(dataset_id=diff_exps[str(vals[4:6])],\
                                         genome_region_id = gene.id,\
                                         value=value, pval=pvalue)
        session.add(diff_exp_data)

    session.flush()
    session.commit()
    session.close()



@timing
def load_gem(chip_peak_analyses, base, datasets, genome):
    gem_path = settings.data_directory+'/chip_peaks/gem/'
    session = base.Session()
    for chip_peak_analysis in chip_peak_analyses:
        if chip_peak_analysis.children[0].protocol_type == 'ChIPchip': continue

        try: gem_peak_file = open(gem_path+chip_peak_analysis.name+'/out_GPS_events.narrowPeak','r')
        except: continue
        for line in gem_peak_file.readlines():
            vals = line.split('\t')

            position = int(vals[3].split(':')[1])

            peak_region = session.get_or_create(base.GenomeRegion, leftpos=vals[1], rightpos=vals[2], strand='+', genome_id=genome.id)

            peak_data = session.get_or_create(datasets.ChIPPeakData, dataset_id=chip_peak_analysis.id, genome_region_id=peak_region.id,\
                                                value=vals[6], eventpos=position, pval=vals[8])

    session.close()



@timing
def load_extra_analyses(base, datasets, genome, analyses_folder, group_name=None):
    session = base.Session()

    analysis_types = [x[0] for x in session.query(func.distinct(datasets.Analysis.type)).all()]

    for file_name in os.listdir(analyses_folder):
        prefix = file_name.split('.')[0]
        vals = prefix.split('_')
        if len(vals) < 6: continue

        strain = session.query(datasets.Strain).filter_by(name=vals[1]).one()

        environment = session.query(datasets.InVivoEnvironment).filter_by(carbon_source=vals[2],
                                                                      nitrogen_source=vals[3],
                                                                      electron_acceptor=vals[4],
                                                                      temperature=37,
                                                                      supplements='').one()


        if vals[7] == 'peaks':
            ce = datasets.ChIPExperiment
            experiments = session.query(ce).filter(and_(ce.strain_id == strain.id,
                                                        ce.environment_id == environment.id,
                                                        ce.antibody == vals[5])).all()

            peak_analysis = session.get_or_create(datasets.ChIPPeakAnalysis,
                                                  name=file_name.split('.')[0],
                                                  environment_id=environment.id,
                                                  strain_id=strain.id,
                                                  parameters=vals[6],
                                                  replicate=1,
                                                  group_name=group_name
                                                 )

            for experiment in experiments:
                session.get_or_create(datasets.AnalysisComposition, analysis_id = peak_analysis.id, dataset_id = experiment.id)

    session.flush()
    session.commit()
    session.close()

@timing
def load_gff_chip_peaks(chip_peak_analyses, base, datasets, genome, group_name):
    gff_path = settings.data_directory+'/chip_peaks/'+group_name
    session = base.Session()
    for chip_peak_analysis in chip_peak_analyses:

        try: gff_peak_file = open(gff_path+'/'+chip_peak_analysis.name+'.gff','r')
        except: continue

        for line in gff_peak_file.readlines():
            if line[0] == '#': continue
            vals = line.split('\t')

            position = (int(vals[3])+int(vals[4]))/2

            peak_region = session.get_or_create(base.GenomeRegion, leftpos=vals[3], rightpos=vals[4], strand='+', genome_id=genome.id)

            peak_data = session.get_or_create(datasets.ChIPPeakData, dataset_id=chip_peak_analysis.id, genome_region_id=peak_region.id,\
                                                value=vals[5], eventpos=position, pval=0.)


    session.flush()
    session.commit()
    session.close()


@timing
def load_arraydata(file_path, datasets, group_name='ec2'):
    array_data_file = open(file_path)

    header = array_data_file.readline().rstrip('\n').split('\t')

    session = base.Session()

    exp_id_map = {}
    for i,name in enumerate(header[2:]):

        try:
            exp_id_map[i] = session.query(datasets.ArrayExperiment).filter(datasets.ArrayExperiment.group_name == group_name).\
                                    filter(func.lower(datasets.ArrayExperiment.name) == str(name[:-4]+'_'+group_name).lower()).one().id
        except: print name

    for line in array_data_file.readlines():
        vals = line.split('\t')


        try: gene = session.query(components.Gene).filter(or_(components.Gene.name == vals[0],\
                                                              components.Gene.locus_id == vals[0])).one()
        except: continue


        for i,val in enumerate(vals[2:]):

            try:
                value = float(val)
                exp_id_map[i]
            except: continue

            array_data = data.GenomeData(dataset_id = exp_id_map[i],
                                         genome_region_id = gene.id,
                                         value=value)
            session.add(array_data)

    session.flush()
    session.commit()
    session.close()


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


@timing
def make_genome_region_map(base, datasets, genome):
    session = base.Session()

    base.GenomeRegionMap.__table__.drop()
    base.GenomeRegionMap.__table__.create()


    #genome_regions = session.query(data.GenomeRegion).all()
    duplicate_set = set()
    tus = session.query(components.TU).all()

    for peak in session.query(datasets.ChIPPeakData).all():
        genome_region_1 = peak.genome_region
        print genome_region_1

        for tu in tus:
            genome_region_2 = tu.genome_region

            if (genome_region_1.id, genome_region_2.id) in duplicate_set: continue
            else: duplicate_set.add((genome_region_1.id, genome_region_2.id))


            right_left_distance = abs(genome_region_1.rightpos - genome_region_2.leftpos)
            left_right_distance = abs(genome_region_1.leftpos - genome_region_2.rightpos)
            if right_left_distance < 1000 and genome_region_2.strand == '+':
            	session.add(base.GenomeRegionMap(genome_region_1.id, genome_region_2.id, right_left_distance))
            elif left_right_distance < 1000 and genome_region_2.strand == '-':
                session.add(base.GenomeRegionMap(genome_region_1.id, genome_region_2.id, left_right_distance))

    session.flush()
    session.commit()
    session.close()

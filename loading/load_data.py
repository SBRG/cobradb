####This code is just a minor modification of sequtil/make_gff.py written by aebrahim####
#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
from PrototypeDB.orm import base
from PrototypeDB.orm import data
from PrototypeDB.lib import settings
from os.path import split
from math import log
import os

from numpy import zeros, roll

import pysam

session = base.Session()

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
    """write samfile object to an output object in a gff format

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


def load_samfile_to_db(sam_filepath, data_set_id, loading_cutoff=0, bulk_file_load=False,
                       flip=False, log2=False, separate_strand=False, include_insert=False,
                       five_prime=False, track=None):
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
    if track is None:
        name = split(samfile.filename)[1]
    else:
        name = track
    #connection to mongo data store
    genome_data = base.omics_database.genome_data

    if log2:
        str_func = lambda x, s: "%.2f" % (log(x, 2) * s)
    else:
        str_func = lambda x, s: "%d" % (x * s)

    for reference in all_counts:
        for strand in all_counts[reference]:
            factor = 1 if strand == "+" else -1
            track_name = "%s_(%s)" % (name, strand) if separate_strand else name
            counts = all_counts[reference][strand]
            entries = []
            for i in counts.nonzero()[0]:    
                if float(counts[i]) < loading_cutoff: continue
            
                entries.append({
                                "position": int(i),
                                "value": float(counts[i]),
                                "strand": strand,
                                "data_set_id": data_set_id})
            
                if i%10000 == 0:
                    genome_data.insert(entries)
                    entries = []
    
    if not bulk_file_load: 
        genome_data.create_index([("data_set_id", ASCENDING), ("position", ASCENDING)])                
    
    samfile.close()


def name_based_experiment_loading(exp_name, lab='palsson', institution='UCSD', bulk_file_load=False):
    vals = exp_name.split('_')
    strain = session.get_or_create(data.Strain, name=vals[1])
    data_source = session.get_or_create(data.DataSource, name=vals[0], lab=lab, institution=institution)
    environment = session.get_or_create(data.InVivoEnvironment, name='_'.join(vals[2:5]), carbon_source=vals[2],\
                                    nitrogen_source=vals[3], electron_acceptor=vals[4], temperature=37)
    exp_type = vals[0].split('-')
    try: vals[6] = vals[6].split('.')[0]
    except: vals[5] = vals[5].split('.')[0]
    if exp_type[0][0:4] == 'chip': 
        experiment = session.get_or_create(data.ChIPExperiment, name='_'.join(vals[0:7]), replicate=vals[5],\
                                           strain=strain, data_source=data_source, environment=environment,\
                                           protocol_type=exp_type[0], antibody=vals[6],\
                                           target=exp_type[1], file_name=exp_name)
        #variance = gff_variance(settings.dropbox_directory+'/crp/data/ChIP/gff/raw/'+exp_name)
        #data.load_genome_data(settings.dropbox_directory+'/crp/data/ChIP/gff/raw/'+exp_name, experiment.id, bulk_file_load,\
        #                      loading_cutoff=math.sqrt(variance)/4.)
    
    elif exp_type[0][0:6] == 'RNAseq':
        experiment = session.get_or_create(data.RNASeqExperiment, name='_'.join(vals[0:7]), replicate=vals[5],\
                                           strain=strain, data_source=data_source, environment=environment,\
                                           machine_id='miseq', sequencing_type='unpaired',\
                                           file_name=exp_name)

        #load_samfile_to_db(settings.dropbox_directory+'/crp/data/RNAseq/bam/'+exp_name, experiment.id, loading_cutoff=10, bulk_file_load=bulk_file_load)

def run_cuffquant(exp):

    os.chdir(settings.dropbox_directory+'/crp/data/RNAseq/cbx')
    gtf_file = settings.dropbox_directory+'/crp/data/annotation/e_coli_notRNA_rRNA.gtf'
    
    out_path = settings.dropbox_directory+'/crp/data/RNAseq/cbx/'+exp.name
    os.system('rm -r '+out_path)
    os.mkdir(out_path)
    os.chdir(out_path)
    exp_file = settings.dropbox_directory+'/crp/data/RNAseq/bam/'+exp.file_name
    
    os.system('%s -p %d %s %s' % ('cuffquant', 8, gtf_file, exp_file))
    
    
def run_cuffnorm(exps):
    gtf_file = settings.dropbox_directory+'/crp/data/annotation/e_coli_notRNA_rRNA.gtf'
    
    out_path = settings.dropbox_directory+'/crp/data/RNAseq/cuffnorm'
    os.system('rm -r '+out_path)
    os.mkdir(out_path)
    os.chdir(out_path)
    #bam_dir = settings.dropbox_directory+'/crp/data/RNAseq/bam/'
    cxb_dir = settings.dropbox_directory+'/crp/data/RNAseq/cxb/'

    os.system('cuffnorm -p %d --library-type fr-firststrand -L %s %s %s' % (24,\
             ','.join([x[1][0][:-2] for x in exps]), gtf_file,\
             ' '.join([','.join([cxb_dir+x.split('.')[0]+'/abundances.cxb' for x in exp[0]]) for exp in exps]))) 
    
    
def run_cuffdiff(exps):
    gtf_file = settings.dropbox_directory+'/crp/data/annotation/e_coli_notRNA_rRNA.gtf'
    
    out_path = settings.dropbox_directory+'/crp/data/RNAseq/cuffdiff'
    os.system('rm -r '+out_path)
    os.mkdir(out_path)
    os.chdir(out_path)
    cxb_dir = settings.dropbox_directory+'/crp/data/RNAseq/cxb/'
    
    os.system('cuffdiff -p %d --library-type fr-firststrand --upper-quartile-norm --FDR 0.05 -L %s %s %s' % (24,\
             ','.join([x[1][0][:-2] for x in exps]), gtf_file,\
             ' '.join([','.join([cxb_dir+x.split('.')[0]+'/abundances.cxb' for x in exp[0]]) for exp in exps])))



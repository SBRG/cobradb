####This code is just a minor modification of sequtil/make_gff.py written by aebrahim####
#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
from PrototypeDB.orm import base
from PrototypeDB.orm import data
from PrototypeDB.orm import components
from PrototypeDB.lib import settings
from os.path import split
from math import log
import os

from numpy import zeros, roll
from sqlalchemy import func

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
        experiment = session.get_or_create(data.RNASeqExperiment, name='_'.join(vals[0:6]), replicate=vals[5],\
                                           strain=strain, data_source=data_source, environment=environment,\
                                           machine_id='miseq', sequencing_type='unpaired',\
                                           file_name=exp_name)
        
    elif exp_type[0][0:7] == 'affyexp':
        experiment = session.get_or_create(data.ArrayExperiment, name='_'.join(vals[0:6]), replicate=vals[5],\
                                           strain=strain, data_source=data_source, environment=environment,\
                                           platform=vals[6], file_name=exp_name)


        #load_samfile_to_db(settings.dropbox_directory+'/crp/data/RNAseq/bam/'+exp_name, experiment.id, loading_cutoff=10, bulk_file_load=bulk_file_load)

def run_cuffquant(exp):

    os.chdir(settings.dropbox_directory+'/crp/data/RNAseq/cxb')
    gtf_file = settings.dropbox_directory+'/crp/data/annotation/e_coli_notRNA_rRNA.gtf'
    
    out_path = settings.dropbox_directory+'/crp/data/RNAseq/cxb/'+exp.name
    os.system('rm -r '+out_path)
    os.mkdir(out_path)
    os.chdir(out_path)
    exp_file = settings.dropbox_directory+'/crp/data/RNAseq/bam/'+exp.file_name
    
    os.system('%s -p %d %s %s' % ('cuffquant', 8, gtf_file, exp_file))
    
    
def run_cuffnorm(exp_sets):
    gtf_file = settings.dropbox_directory+'/crp/data/annotation/e_coli_notRNA_rRNA.gtf'
    cxb_dir = settings.dropbox_directory+'/crp/data/RNAseq/cxb/'
    out_path = settings.dropbox_directory+'/crp/data/RNAseq/cuffnorm'
    
    os.system('rm -r '+out_path)
    os.mkdir(out_path)
    os.chdir(out_path)

    os.system('cuffnorm -p %d --library-type fr-firststrand --library-norm-method classic-fpkm -L %s %s %s' % (24,\
             ','.join([x[0][0][:-2] for x in exp_sets]), gtf_file,\
             ' '.join([','.join([cxb_dir+x.split('.')[0]+'/abundances.cxb' for x in exp[0]]) for exp in exp_sets]))) 
    
    
def run_cuffdiff(exp_sets):
    gtf_file = settings.dropbox_directory+'/crp/data/annotation/e_coli_notRNA_rRNA.gtf'
    cxb_dir = settings.dropbox_directory+'/crp/data/RNAseq/cxb/'
    out_path = settings.dropbox_directory+'/crp/data/RNAseq/cuffdiff'
    
    
    os.system('rm -r '+out_path)
    os.mkdir(out_path)
    os.chdir(out_path)
    
    os.system('cuffdiff -v -p %d --library-type fr-firststrand --upper-quartile-norm --FDR 0.05 -L %s %s %s' % (24,\
             ','.join([x[0][0][:-2] for x in exp_sets]), gtf_file,\
             ' '.join([','.join([cxb_dir+x.split('.')[0]+'/abundances.cxb' for x in exp[0]]) for exp in exp_sets])))
    
    #print 'cuffdiff -p %d --library-type fr-firststrand --upper-quartile-norm --FDR 0.05 -L %s -C %s %s %s' % (24,\
    #         ','.join([x[0][0][:-2] for x in exp_sets]), out_path+'/contrasts2.txt', gtf_file,\
    #         ' '.join([','.join([cxb_dir+x.split('.')[0]+'/abundances.cxb' for x in exp[0]]) for exp in exp_sets]))
    

def run_gem(exp_sets):
    gem_path = settings.home_directory+'/libraries/gem'
    bam_dir = settings.dropbox_directory+'/crp/data/ChIP/bam/'
    """ This could easily be parallelized """
    for exp in exp_sets:
        outdir = '_'.join(exp[0][0].split('_')[0:5]+exp[0][0].split('_')[6:7])    
        out_path = settings.dropbox_directory+'/crp/data/ChIP_peaks/gem/'+outdir
        os.system('rm -r '+out_path)
        os.mkdir(out_path)
        os.chdir(out_path)
    
        input_files = ' '.join(['--expt'+exp[0][i]+' '+bam_dir+exp[1][i] \
                                for i in range(len(exp[0]))])

        os.system("java -Xmx5G -jar %s/gem.jar --d %s/Read_Distribution_ChIP-exo.txt --g %s --genome %s %s --f SAM --k_min %d --smooth 3 --outNP" %\
                  (gem_path, gem_path, settings.dropbox_directory+'/crp/data/annotation/ec_mg1655.sizes', settings.dropbox_directory+'/crp/data/annotation',\
                   input_files, 1))
            
                                    
def load_cuffnorm():
    cuffnorm_genes = open(settings.dropbox_directory+'/crp/data/RNAseq/cuffnorm/isoforms.attr_table','r')
    cuffnorm_genes.readline()
    cuffnorm_gene_dict = {}
    for line in cuffnorm_genes.readlines():
        vals = line.split('\t')
        cuffnorm_gene_dict[vals[0]] = {'name':vals[4],'leftpos':vals[6].split(':')[1].split('-')[0],\
                                                  'rightpos':vals[6].split(':')[1].split('-')[1]}
        
    cuffnorm_output = open(settings.dropbox_directory+'/crp/data/RNAseq/cuffnorm/isoforms.fpkm_table','r')
    header = cuffnorm_output.readline().rstrip('\n').split('\t')
    print header
    exp_id_map = {i:session.query(data.RNASeqExperiment).filter_by(name=name[:-2]+'_'+str(int(name[-1])+1)).one().id\
                  for i,name in enumerate(header[1:])} 

    for line in cuffnorm_output.readlines():
        vals = line.split('\t')
        
        try: gene = session.query(components.Gene).filter_by(locus_id=vals[0]).one()
        except: 
            x = cuffnorm_gene_dict[vals[0]]
            gene = session.get_or_create(base.GenomeRegion, name=x['name'], leftpos=x['leftpos'],\
                                                            rightpos=x['rightpos'], strand='+')

        for i,val in enumerate(vals[1:]):
        
            try: value = float(val)
            except: continue
        
            session.get_or_create(data.GenomeData, data_set_id=exp_id_map[i-1],\
                                                   genome_region_id = gene.id,\
                                                   value=value)
    
    
def load_cuffdiff():
    cuffdiff_output = open(settings.dropbox_directory+'/crp/data/RNAseq/cuffdiff/gene_exp.diff','r')
    header = cuffdiff_output.readline()
    diff_exps = {}
    for line in cuffdiff_output.readlines():
        vals = line.split('\t')
        
        try: 
            value = float(vals[9])
            pvalue = float(vals[11])
        except: continue
        
        if pvalue > .25: continue
        
        if str(vals[4:6]) not in diff_exps.keys():
            diff_exp = session.get_or_create(data.DifferentialExpression, name=vals[4]+'\\'+vals[5], norm_method='classic-fpkm',fdr=.05)
            print vals[4]
            exp1 = session.query(data.Analysis).filter_by(name=vals[4]).one()
            exp2 = session.query(data.Analysis).filter_by(name=vals[5]).one()
            session.get_or_create(data.AnalysisComposition, analysis_id = diff_exp.id, data_set_id = exp1.id)
            session.get_or_create(data.AnalysisComposition, analysis_id = diff_exp.id, data_set_id = exp2.id)
            diff_exps[str(vals[4:6])] = diff_exp.id

        try: gene = session.query(components.Gene).filter_by(locus_id=vals[0]).one()
        except: 
            x = {'name':vals[2], 'leftpos':vals[3].split(':')[1].split('-')[0],\
                                 'rightpos':vals[3].split(':')[1].split('-')[1]}
                                                                             
            gene = session.get_or_create(base.GenomeRegion, name=x['name'], leftpos=x['leftpos'],\
                                                            rightpos=x['rightpos'], strand='+')

        
        
    
        session.get_or_create(data.DiffExpData, data_set_id=diff_exps[str(vals[4:6])],\
                                               genome_region_id = gene.id,\
                                               value=value, pval=pvalue)


def load_gem():
    gem_path = settings.dropbox_directory+'/crp/data/ChIP_peaks/gem/'
    
    for exp_name,exp_id in session.query(data.ChIPPeak.name,data.ChIPPeak.id).all():
        gem_peak_file = open(gem_path+exp_name+'/out_GPS_events.narrowPeak','r')
        
        for line in gem_peak_file.readlines():
            vals = line.split('\t')

            position = int(vals[3].split(':')[1])
        
            with open(gem_path+exp_name+'/'+exp_name+'_peaks.gff', 'wb') as peaks_gff_file:
                peaks_gff_file.write('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n' %\
                                     ('NC_000913','.',exp_name+'_peaks', position-1, position+1, float(vals[6])+20, '+', '.','.'))
                peaks_gff_file.write('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n' %\
                                     ('NC_000913','.',exp_name+'_peaks', vals[1], vals[2], float(vals[6]), '+', '.','.'))
            
            peak_region = session.get_or_create(base.GenomeRegion, leftpos=vals[1], rightpos=vals[2], strand='+')

            
            
            peak_data = session.get_or_create(data.ChIPPeakData, data_set_id=exp_id, genome_region_id=peak_region.id,\
                                              value=vals[6], eventpos=position, pval=vals[8])


def load_arraydata(file_path):
    array_data_file = open(file_path)
    
    header = array_data_file.readline().rstrip('\n').split('\t')

    exp_id_map = {}
    for i,name in enumerate(header[2:]):
        try: exp_id_map[i] = session.query(data.ArrayExperiment).filter(func.lower(data.ArrayExperiment.name) == name[:-4].lower()).one().id
        except: 
            print name
            continue
             
    for line in array_data_file.readlines():
        vals = line.split('\t')
        
        ##This code sucks right now and depend on components or load_cuffdiff/load_cuffnorm 
        ##above being run first
        try: gene = session.query(components.Gene).filter_by(locus_id=vals[0]).one()
        except: 
            try: gene = session.query(base.GenomeRegion).filter_by(name=vals[0]).one()
            except: 
                print 'fuck  '+str(vals[1])
                continue
        
        for i,val in enumerate(vals[2:]):
        
            try: value = float(val)
            except: continue
            try:
                session.get_or_create(data.GenomeData, data_set_id=exp_id_map[i],\
                                                       genome_region_id = gene.id,\
                                                       value=value) 
            except: None
        
    


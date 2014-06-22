####A portion of this code is derived from sequtil/make_gff.py written by aebrahim####
#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
from om import base, data, components, settings, timing
from os.path import split
from math import log
from itertools import combinations

import os

from numpy import zeros, roll
from sqlalchemy import func, or_

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
        #if i > 1e3: break
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
                                "data_set_id": data_set_id})

                if i%50000 == 0:
                    genome_data.insert(entries)
                    entries = []
            genome_data.insert(entries)
    if not bulk_file_load:
        genome_data.create_index([("data_set_id", ASCENDING), ("leftpos", ASCENDING)])

    samfile.close()

@timing
def name_based_experiment_loading(exp_name, lab='palsson', institution='UCSD',
								  bulk_file_load=False, norm_factor=1.):
    vals = exp_name.split('_')
    if len(vals) < 5: return

    exp_type = vals[0].split('-')
    try:
        vals[7] = vals[7].split('.')[0]
        name = '_'.join(vals[0:8])
        supplements = vals[7]
    except:
        try: vals[6] = vals[6].split('.')[0]
        except: vals[5] = vals[5].split('.')[0]
        name = '_'.join(vals[0:7])
        supplements = ''

    if norm_factor == 1.:
        norm_method = None
    else: norm_method = 'mean of total mapped reads'

    session = base.Session()

    strain = session.get_or_create(data.Strain, name=vals[1])
    data_source = session.get_or_create(data.DataSource, name=vals[0], lab=lab, institution=institution)
    environment = session.get_or_create(data.InVivoEnvironment, name='_'.join(vals[2:5]+[supplements]), carbon_source=vals[2],\
                                        nitrogen_source=vals[3], electron_acceptor=vals[4], temperature=37,\
                                        supplements=supplements)


    if exp_type[0][0:4] == 'chip':

        experiment = session.get_or_create(data.ChIPExperiment, name=name+'_'+str(norm_factor), replicate=vals[5],\
                                           strain_id=strain.id, data_source_id=data_source.id, environment_id=environment.id,\
                                           protocol_type=exp_type[0], antibody=vals[6],\
                                           target=exp_type[1], file_name=exp_name, normalization_method=norm_method,\
                                                                                   normalization_factor=norm_factor)

        load_samfile_to_db(settings.dropbox_directory+'/crp/data/ChIP/bam/'+exp_name, experiment.id, loading_cutoff=5,\
                           bulk_file_load=bulk_file_load, five_prime=True, norm_factor=norm_factor)



    elif exp_type[0][0:6] == 'RNAseq':
        experiment = session.get_or_create(data.RNASeqExperiment, name=name+'_'+str(norm_factor), replicate=vals[5],\
                                           strain_id=strain.id, data_source_id=data_source.id, environment_id=environment.id,\
                                           machine_id='miseq', sequencing_type='unpaired',\
                                           file_name=exp_name, normalization_method=norm_method,\
                                                               normalization_factor=norm_factor)

        load_samfile_to_db(settings.dropbox_directory+'/crp/data/RNAseq/bam/'+exp_name, experiment.id, loading_cutoff=10,\
                          bulk_file_load=bulk_file_load, five_prime=False, flip=True, norm_factor=norm_factor)

    elif exp_type[0][0:7] == 'affyexp':
        experiment = session.get_or_create(data.ArrayExperiment, name=name, replicate=vals[5],\
                                           strain_id=strain.id, data_source_id=data_source.id, environment_id=environment.id,\
                                           platform=vals[6], file_name=exp_name)
    session.flush()
    session.commit()
    session.close()


def run_cuffquant(exp):

    os.chdir(settings.dropbox_directory+'/crp/data/RNAseq/cxb')
    gtf_file = settings.dropbox_directory+'/crp/data/annotation/e_coli_notRNA_rRNA.gtf'

    out_path = settings.dropbox_directory+'/crp/data/RNAseq/cxb/'+exp.name
    os.system('rm -r '+out_path)
    os.mkdir(out_path)
    os.chdir(out_path)
    exp_file = settings.dropbox_directory+'/crp/data/RNAseq/bam/'+exp.file_name

    os.system('%s -p %d %s %s' % ('cuffquant', 8, gtf_file, exp_file))

@timing
def run_cuffnorm(exp_sets):
    gtf_file = settings.dropbox_directory+'/crp/data/annotation/e_coli_notRNA_rRNA.gtf'
    cxb_dir = settings.dropbox_directory+'/crp/data/RNAseq/cxb/'
    out_path = settings.dropbox_directory+'/crp/data/RNAseq/cuffnorm'

    os.system('rm -r '+out_path)
    os.mkdir(out_path)
    os.chdir(out_path)

    os.system('cuffnorm -p %d --library-type fr-firststrand -L %s %s %s' % (24,\
             ','.join([x[0][0][:-2] for x in exp_sets]), gtf_file,\
             ' '.join([','.join([cxb_dir+x.split('.')[0]+'/abundances.cxb' for x in exp[0]]) for exp in exp_sets])))


def find_single_factor_pairwise_contrasts(data_sets):
    data_set_contrasts = []
    data_set_conditions = {}
    for data_set in data_sets:
        data_set_conditions[(data_set.strain, data_set.environment)] = data_set
    for c in combinations(data_set_conditions.keys(), 2):
        e1, e2 = sorted(c, key=str)
        s1, c1 = e1
        s2, c2 = e2
        if s1 == s2:  # if the strains are the same
           # single shift only
            differences = (c1.carbon_source != c2.carbon_source) + \
                          (c1.nitrogen_source != c2.nitrogen_source) + \
                          (c1.electron_acceptor != c2.electron_acceptor)
            if differences != 1:
                continue
        else:  # if the strains are different
            if not (s1.name == "wt" or s2.name == "wt"):
                continue
            if c1 != c2:  # make sure the conditions are the same
                continue
        exp_1 = data_set_conditions[(s1, c1)]
        exp_2 = data_set_conditions[(s2, c2)]
        data_set_contrasts.append([exp_1,exp_2])
    return data_set_contrasts


def generate_cuffdiff_contrasts(normalized_expression_objects):
    contrasts = find_single_factor_pairwise_contrasts(normalized_expression_objects)
    with open(settings.dropbox_directory+'/crp/data/RNAseq/cuffdiff/contrasts.txt', 'wb') as contrast_file:
        contrast_file.write('condition_A\tcondition_B\n')
        for contrast in contrasts:
            contrast_file.write(contrast[0].name+'\t'+contrast[1].name+'\n')


@timing
def run_cuffdiff(normalized_expression_objects, debug=False):
    gtf_file = settings.dropbox_directory+'/crp/data/annotation/e_coli_notRNA_rRNA.gtf'
    cxb_dir = settings.dropbox_directory+'/crp/data/RNAseq/cxb/'
    out_path = settings.dropbox_directory+'/crp/data/RNAseq/cuffdiff'

    os.system('rm -r '+out_path)
    os.mkdir(out_path)
    os.chdir(out_path)

    generate_cuffdiff_contrasts(normalized_expression_objects)

    if debug:
        print settings.cufflinks+'/cuffdiff -v -p %d --library-type fr-firststrand --FDR 0.05 -C %s -L %s %s %s' % (24, out_path+'/contrasts.txt',
                      ','.join([x.name for x in normalized_expression_objects]), gtf_file,
                      ' '.join([','.join([cxb_dir+x.file_name.split('.')[0]+'/abundances.cxb' for x in exp.children]) for exp in normalized_expression_objects]))
    else:
        os.system(settings.cufflinks+'/cuffdiff -v -p %d --library-type fr-firststrand --FDR 0.05 -C %s -L %s %s %s' % (24, out_path+'/contrasts.txt',
                      ','.join([x.name for x in normalized_expression_objects]), gtf_file,
                      ' '.join([','.join([cxb_dir+x.file_name.split('.')[0]+'/abundances.cxb' for x in exp.children]) for exp in normalized_expression_objects])))


@timing
def run_gem(chip_peak_analysis, control_chip_peak_analysis = None, debug=False):
    default_parameters = {'mrc':20, 'smooth':3, 'nrf':'', 'outNP':''}
    gem_path = settings.home_directory+'/libraries/gem'
    bam_dir = settings.dropbox_directory+'/crp/data/ChIP/bam/'
    """ This could easily be parallelized """


    outdir = chip_peak_analysis.name
    out_path = settings.dropbox_directory+'/crp/data/ChIP_peaks/gem/'+outdir
    os.system('rm -r '+out_path)
    os.mkdir(out_path)
    os.chdir(out_path)

    input_files = ' '.join(['--expt'+x.name+' '+bam_dir+x.file_name for x in chip_peak_analysis.children])

    if control_chip_peak_analysis:
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


    if debug:
        print "java -Xmx5G -jar %s/gem.jar --d %s/Read_Distribution_ChIP-exo.txt --g %s --genome %s %s %s --f SAM %s" %\
                     (gem_path, gem_path, settings.dropbox_directory+'/crp/data/annotation/ec_mg1655.sizes', settings.dropbox_directory+'/crp/data/annotation',\
                     input_files, control_input_files, parameter_string)

    else:
        os.system("java -Xmx5G -jar %s/gem.jar --d %s/Read_Distribution_ChIP-exo.txt --g %s --genome %s %s %s --f SAM %s" %\
                      (gem_path, gem_path, settings.dropbox_directory+'/crp/data/annotation/ec_mg1655.sizes', settings.dropbox_directory+'/crp/data/annotation',\
                       input_files, control_input_files, parameter_string))


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


@timing
def load_cuffnorm():
    cuffnorm_genes = open(settings.dropbox_directory+'/crp/data/RNAseq/cuffnorm/isoforms.attr_table','r')
    cuffnorm_genes.readline()
    cuffnorm_gene_dict = {}
    session = base.Session()
    for line in cuffnorm_genes.readlines():
        vals = line.split('\t')
        cuffnorm_gene_dict[vals[0]] = {'name':vals[4],'leftpos':vals[6].split(':')[1].split('-')[0],\
                                                  'rightpos':vals[6].split(':')[1].split('-')[1]}

    cuffnorm_output = open(settings.dropbox_directory+'/crp/data/RNAseq/cuffnorm/isoforms.fpkm_table','r')
    header = cuffnorm_output.readline().rstrip('\n').split('\t')
    #print [name[:-5]+'_'+str(int(name[-1])+1)+'_1.0' for i,name in enumerate(header[1:])]
    exp_id_map = {i:session.query(data.RNASeqExperiment).filter_by(name=name[:-2]+'_'+str(int(name[-1])+1)+'_1.0').one().id\
                  for i,name in enumerate(header[1:])}

    for line in cuffnorm_output.readlines():
        vals = line.split('\t')

        try: gene = session.query(components.Gene).filter_by(locus_id=vals[0]).one()
        except: continue
            #x = cuffnorm_gene_dict[vals[0]]
            #gene = session.get_or_create(base.GenomeRegion, name=x['name'], leftpos=x['leftpos'],\
            #                                                rightpos=x['rightpos'], strand='+')

        for i,val in enumerate(vals[1:]):

            try: value = float(val)
            except: continue

            session.get_or_create(data.GenomeData, data_set_id=exp_id_map[i-1],\
                                                   genome_region_id = gene.id,\
                                                   value=value)
    session.close()

@timing
def load_cuffdiff():
    cuffdiff_output = open(settings.dropbox_directory+'/crp/data/RNAseq/cuffdiff/gene_exp.diff','r')
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

        if str(vals[4:6]) not in diff_exps.keys():
            x = vals[4].split('_')
            y = vals[5].split('_')
            exp_name = ''
            for i in range(len(x)):
                if x[i] == y[i]:
                    exp_name += x[i]+'_'
                else: exp_name += x[i]+'/'+y[i]+'_'
            exp_name = exp_name.rstrip('_')

            diff_exp = session.get_or_create(data.DifferentialExpression, name=exp_name, norm_method='classic-fpkm',fdr=.05)

            exp1 = session.query(data.Analysis).filter_by(name=vals[4]).one()
            exp2 = session.query(data.Analysis).filter_by(name=vals[5]).one()
            session.get_or_create(data.AnalysisComposition, analysis_id = diff_exp.id, data_set_id = exp1.id)
            session.get_or_create(data.AnalysisComposition, analysis_id = diff_exp.id, data_set_id = exp2.id)
            diff_exps[str(vals[4:6])] = diff_exp.id

        try: gene = session.query(components.Gene).filter_by(locus_id=vals[0]).one()
        except: continue
            #x = {'name':vals[2], 'leftpos':vals[3].split(':')[1].split('-')[0],\
            #                     'rightpos':vals[3].split(':')[1].split('-')[1]}

            #gene = session.get_or_create(base.GenomeRegion, name=x['name'], leftpos=x['leftpos'],\
            #                                                rightpos=x['rightpos'], strand='+')




        session.get_or_create(data.DiffExpData, data_set_id=diff_exps[str(vals[4:6])],\
                                               genome_region_id = gene.id,\
                                               value=value, pval=pvalue)

	session.close()

@timing
def load_gem(chip_peak_analyses):
    gem_path = settings.dropbox_directory+'/crp/data/ChIP_peaks/gem/'
    session = base.Session()
    for chip_peak_analysis in chip_peak_analyses:
        gem_peak_file = open(gem_path+chip_peak_analysis.name+'/out_GPS_events.narrowPeak','r')

        for line in gem_peak_file.readlines():
            vals = line.split('\t')

            position = int(vals[3].split(':')[1])

            peak_region = session.get_or_create(base.GenomeRegion, leftpos=vals[1], rightpos=vals[2], strand='+')

            peak_data = session.get_or_create(data.ChIPPeakData, data_set_id=chip_peak_analysis.id, genome_region_id=peak_region.id,\
                                                value=vals[6], eventpos=position, pval=vals[8])

    session.close()

@timing
def load_arraydata(file_path, type='ec2'):
    array_data_file = open(file_path)

    header = array_data_file.readline().rstrip('\n').split('\t')

    session = base.Session()

    exp_id_map = {}
    for i,name in enumerate(header[2:]):
        try:
            exp_id_map[i] = session.query(data.ArrayExperiment).\
                                    filter(func.lower(data.ArrayExperiment.name) == str(name[:-4]+'_'+type).lower()).one().id
        except: continue

    for line in array_data_file.readlines():
        vals = line.split('\t')


        try: gene = session.query(components.Gene).filter(or_(components.Gene.name == vals[0],\
                                                         components.Gene.locus_id == vals[0])).one()
        except: continue


        for i,val in enumerate(vals[2:]):

            try: value = float(val)
            except: continue
            try:
                session.get_or_create(data.GenomeData, data_set_id=exp_id_map[i],\
                                                       genome_region_id = gene.id,\
                                                       value=value)
            except: None

    session.close()

@timing
def make_genome_region_map():
    session = base.Session()

    data.GenomeRegionMap.__table__.drop()
    data.GenomeRegionMap.__table__.create()

    genome_regions = session.query(data.GenomeRegion).all()

    for genome_region_1 in genome_regions:
        #print genome_region_1
        for genome_region_2 in genome_regions:
            if genome_region_1.id == genome_region_2.id: continue

            midpoint_1 = (genome_region_1.leftpos + genome_region_1.rightpos)/2
            midpoint_2 = (genome_region_2.leftpos + genome_region_2.rightpos)/2
            midpoint_distance = midpoint_1 - midpoint_2
            left_right_distance = genome_region_1.leftpos - genome_region_2.rightpos
            right_left_distance = genome_region_1.rightpos - genome_region_2.leftpos
            if abs(midpoint_distance) < 1000 or abs(left_right_distance) < 1000 or abs(right_left_distance) < 1000:
                session.add(data.GenomeRegionMap(genome_region_1.id, genome_region_2.id, midpoint_distance))

    session.commit()
    session.close()

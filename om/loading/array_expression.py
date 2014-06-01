from os.path import join, isfile
from os import unlink
from subprocess import call
from re import compile

from numpy import array, mean, genfromtxt, zeros_like, arange
from sqlalchemy import distinct
from scipy.stats import ttest_ind

from trnlib import settings, timing
from trnlib.omeORM import *
from itertools import combinations


normalize_script = join(settings.trn_directory, "database", "loading", "normalize.R")

def create_translator(session, annotation_file=settings.data_directory + "annotation/GPL3154.tab"):
    """create dict to translate from affy id to database gene wid"""
    affyid2wid = {}
    bnum2wid = dict(session.query(Gene.bnum, Gene.wid))
    bnum_search = compile(r"b\d{4}")
    with open(annotation_file) as infile:
        for line in infile.readlines():
            if line[0] == '#':
                continue
            vals = line.split('\t')
            if bnum_search.match(vals[1]):
                try:
                    affyid2wid[vals[0]] = bnum2wid[vals[1]]
                except KeyError:
                    None
                    # should not actually look for synonyms
                    #synonym = session.search_by_synonym(Gene,
                    #    vals[1] + " (obsolete)").first()
                    #if synonym is not None:
                    #    affyid2wid[vals[0]] = synonym.wid
    return affyid2wid, bnum2wid

def create_experiments(session, header_line, experiment_set, platform):
    """create the appropriate experiments
    
    input: header_line from normalize script output with CEL filenames
    output: list of created experiment wids"""
    experiments = []
    cols = header_line.rstrip('\n').split('\t')
    for col in cols:
        if col == '': continue
        
        vals = col.rstrip('.CEL').split('_')

        exp_type = vals[0]
        strain = vals[1]
        csource = vals[2]
        nsource = vals[3]
        eacceptor = vals[4]
        replicate = int(vals[5])
        
        condition = get_or_create_condition(session, carbon_source=csource,
            nitrogen_source=nsource, eacceptor=eacceptor)
        strain_object = get_or_create(session, Strain, name=strain)
        # link knockout correctly
        knockout_gene = session.query(Gene).filter_by(
            name=strain.replace("delta-", "")).first()
        if knockout_gene is not None:
            if knockout_gene not in strain_object.knockouts:
                strain_object.knockouts.append(knockout_gene)
                session.add(strain_object)
                session.commit()
    
        tmp = strain.split('-')
        try: target = tmp[1]
        except: target = tmp[0]
        
        # TODO make sure the experiment doesn't already exist
        # if session.query(Dataset).filter_by(name=col).first() is not None:
            # print "already in", session.query(Dataset).filter_by(name=col).first()#raise Exception("already in the database")
        dataset = Dataset(name=col)
        session.add(dataset)
        session.commit()
        
        experiment = ArrayExperiment()
        experiment.platform = platform
        experiment.condition = condition
        experiment.strain = strain_object
        experiment.replicate = replicate
        experiment.name = col  # todo - make link to dataset, and use experiment set
        experiment.experiment_set = experiment_set
        session.add(experiment)
        session.commit()
        assert experiment.wid is not None
        experiments.append(experiment)
    return experiments

def calculate_differential_expression(session, condition1, condition2, strain1, strain2, platform, experiment_set):
    """calculate differential expression (fold change and q)"""
    query_args = {"platform": platform, "experiment_set": experiment_set}
    bnums = array([i[0] for i in session.query(distinct(ExpressionData.bnum)).\
        filter_by(**query_args)])
    n_genes = len(bnums)
    # query data
    d1_query = session.query(ExpressionData.value).filter_by(condition=condition1, strain=strain1, **query_args).order_by(ExpressionData.bnum)
    d2_query = session.query(ExpressionData.value).filter_by(condition=condition2, strain=strain2, **query_args).order_by(ExpressionData.bnum)
    try:
        d1 = array(d1_query.all()).reshape(n_genes, -1)
        d2 = array(d2_query.all()).reshape(n_genes, -1)
    except:
        from IPython import embed; embed()
    # filter the data, using the average of all values in the platform file as a cutoff
    cutoff = mean(genfromtxt("%s/expression/CEL/IG_formatted_%s.tab" % (settings.data_directory, platform))[:, 1:])
    selection = (d1.max(axis=1) > cutoff) + (d2.max(axis=1) > cutoff)  # at least one of the samples must have one value over the cutoff
    d1 = d1[selection, :]
    d2 = d2[selection, :]
    bnums = bnums[selection]
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
    return bnums, fold_change, q

def process_differential_expression(session, experiment_set, platform, experiments=None):
    """determine which experiments to compare with differential expression, then perform them

    if no experiments are given, use all of the experiments in the set"""
    if experiments is None:
        experiments = experiment_set.array_experiments

    experiment_conditions = {}
    for experiment in experiments:
        experiment_conditions[(experiment.strain, experiment.condition)] = experiment.wid
    for c in combinations(experiment_conditions.keys(), 2):
        e1, e2 = sorted(c, key=str)
        s1, c1 = e1
        s2, c2 = e2
        if s1 == s2:  # if the strains are the same
            # single shift only
            differences = (c1.carbon_source != c2.carbon_source) + \
                (c1.nitrogen_source != c2.nitrogen_source) + \
                (c1.eacceptor != c2.eacceptor) + \
                (c1.other != c2.other) + (c1.temperature != c2.temperature)
            if differences != 1:
                continue
        else:  # if the strains are different
            if not (s1.name == "wt" or s2.name == "wt"):
                continue
            if c1 != c2:  # make sure the conditions are the same
                continue
        exp_wid_1 = experiment_conditions[(s1, c1)]
        exp_wid_2 = experiment_conditions[(s2, c2)]
        r = calculate_differential_expression(session, c1, c2, s1, s2, platform, experiment_set)
        upload_differential_expression(session, exp_wid_1, exp_wid_2, *r)

def upload_differential_expression(session, exp_wid_1, exp_wid_2, bnums, fold_change, fdr):
    experiments = {"experiment_wid_1": exp_wid_1, "experiment_wid_2": exp_wid_2}
    for i in range(len(bnums)):
        fc_data = ArrayAnalysis(bnum=bnums[i], type="fold_change", value=fold_change[i], **experiments)
        fdr_data = ArrayAnalysis(bnum=bnums[i], type="ttest", value=fdr[i], **experiments)
        session.add(fc_data)
        session.add(fdr_data)
    session.commit()

@timing
def process_array_expression_data(session, folder, platform, experiment_set=None, differential_expression=True):
    # first need to normalize the data
    normalized_file_path = join(folder, "normalized_data.txt")
    if isfile(normalized_file_path):
        unlink(normalized_file_path)
    call([settings.Rscript, "--slave", normalize_script], cwd=folder)
    # make sure the R script worked
    if not isfile(normalized_file_path):
        raise Exception("R script failed to run")
    # create dict to map from 
    affyid2wid, bnum2wid = create_translator(session)
    datafile = open(normalized_file_path)
    header = datafile.readline()
    experiments = create_experiments(session, header, experiment_set, platform)
    # process the datafiles and upload the data
    re_search = compile(r"b\d{4}")
    for line in datafile.readlines():
        gene_wid = None
        vals = line.strip().split('\t')
        if platform == "asv2":
            found_bnums = re_search.findall(vals[0])
            if len(found_bnums) == 1:
                try:
                    gene_wid = bnum2wid[found_bnums[0]]
                except:
                    continue
        else:
            try:
                gene_wid = affyid2wid[vals[0]]
            except KeyError:
                continue  # not a gene
        if gene_wid is None:
            continue
        query_str = "INSERT INTO array_data(value, experiment_WID) VALUES"
        for i, experiment in enumerate(experiments):
            query_str += " (%s, %i)," % (vals[i + 1], experiment.wid)
        query_str = query_str.strip(",") + " RETURNING wid;"
        array_wids = [i[0] for i in session.execute(query_str).fetchall()]
        query_str = "INSERT INTO array_mapping(gene_WID, array_WID) VALUES"
        for array_wid in array_wids:
            query_str += " (%i, %i)," % (gene_wid, array_wid)
        query_str = query_str.strip(",") + ";"
        session.execute(query_str)
    datafile.close()
    session.commit()
    unlink(normalized_file_path)

    # set up differential expression tests
    if differential_expression:
        process_differential_expression(session, experiment_set, platform, experiments=experiments)

if __name__ == "__main__":
    import os
    print "starting reset"
    os.system("%s < ../schemas/kb_schema_ome_affy.sql > psql.log 2>&1" % (settings.psql_full))
    os.system("%s < ../schemas/kb_schema_ome_views.sql > psql.log 2>&1" % (settings.psql_full))
    print "reset done"
    from trnlib.omeORM import *
    session = Session()
    session.execute("delete from datasets where name like 'affyexp%';")
    session.commit()

    ec2_dataset = session.get_or_create(Dataset,
        name="ec2 affy data", lab="sbrg", institution="ucsd")
    ec2_exp_set = session.get_or_create(ExperimentSet,
        dataset=ec2_dataset, type="affy")
    process_array_expression_data(session, folder, "ec2", experiment_set=ec2_exp_set)

from om import base, settings
from om.data import *
from om.components import *

from IPython.display import HTML
from matplotlib import pylab as plt
import pandas as pd
import cobra
import cPickle as pickle

model = pickle.load(open(settings.data_directory+'/models/iJO1366.pickle', "rb"))
#model = cobra.io.load_matlab_model(settings.data_directory+'/models/iJO1366')



def get_chip_peak_gene_expression(gene_name, factor, condition):

    #if gene_name = 'nan': return None

    chip_regulation = ome.query(ChipPeakGene).filter_by(target=factor, condition=condition, gene_name=str(gene_name)).all()

    if not chip_regulation: return None
    aa = AllAnalysis
    target = 'delta-'+factor[0].lower()+factor[1:]
    array_regulation = ome.query(aa).filter(or_(and_(aa.target1 == target, aa.target2 == 'wt'),\
                                                and_(aa.target2 == target, aa.target1 == 'wt')),\
                                            and_(aa.gene_name == gene_name, aa.condition1 == condition, aa.fdr < .05)).all()

    if len(array_regulation) > 1: return 'error'
    elif len(array_regulation) == 0: return None

    if array_regulation[0].target1 == 'wt': return {'fold_change':float(array_regulation[0].fold_change)}
    else: return {'fold_change': -1*float(array_regulation[0].fold_change)}




def get_regulation_data(cobra_rxn_id, expanded_dataset=False):
    """This function takes a cobra_rxn_id and returns a json regulation structure, by default it works for chip peak gene expression
       and the numeric value always corresponds to positive(activation) and negative(repression)
    {'cobra_rxn_id': [{'value':1., 'type':'transcriptional', 'regulator':'component_1', 'gene':'gene_id_1', 'dataset':dataset_id_1},
                      {'value':-5, 'type':'transcriptional', 'regulator':'component_2', 'gene':'gene_id_2', 'dataset':dataset_id_1},
                      {'value':3., 'type':'translational', 'regulator':'component_id_3', 'mRNA':'rna_id_36', 'dataset':dataset_id_2},
                       ...
                      {'value':n, 'type':'allosteric', 'regulator':'component_id_n', 'protein':'protein_id_n', 'dataset':dataset_id_n}]}
    """
    try:
        cobra_rxn_id.id
        rxn = cobra_rxn_id
    except:
        rxn = model.reactions.get_by_id(cobra_rxn_id)

    cpge = ChIPPeakGeneExpression

    gpr_genes = [ome.query(Gene.name).filter_by(locus_id=gene.id).scalar() for gene in rxn.genes]

    gpr_regulation = ome.query(cpge.expression_value, cpge.target, cpge.gene_name,cpge.diff_exp_id,
                               cpge.carbon_source, cpge.nitrogen_source, cpge.electron_acceptor).\
                                 filter(cpge.gene_name.in_(gpr_genes)).\
                                 group_by(cpge.expression_value, cpge.target, cpge.gene_name,
                                          cpge.diff_exp_id, cpge.carbon_source, cpge.nitrogen_source,
                                          cpge.electron_acceptor).all()


    if not expanded_dataset:
        return [{'value':x.expression_value, 'type':'transcriptional', 'regulator':x.target, 'gene':x.gene_name,
                 'dataset':x.diff_exp_id} for x in gpr_regulation]
    else:
        return [{'value':x.expression_value, 'type':'transcriptional', 'regulator':x.target, 'gene':x.gene_name,
                 'dataset':x.diff_exp_id, 'carbon_source':x.carbon_source, 'nitrogen_source':x.nitrogen_source,
                 'electron_acceptor':x.electron_acceptor} for x in gpr_regulation]



def get_indirect_regulation_frame(rxn, factors=['ArcA','Fnr'], condition=None):
    multi_index = [[],[]]
    for factor in factors:
        for gene in rxn._genes.keys():
            if gene == 's0001': continue
            multi_index[0].append(factor)
            multi_index[1].append(ome.query(Gene.name).filter_by(bnum=gene.id).scalar())

    df = DataFrame(index=multi_index, columns=['value'])
    for factor,gene_name in df.index:
        aa = AllAnalysis
        target = 'delta-'+factor[0].lower()+factor[1:]
        array_regulation = ome.query(aa).filter(or_(and_(aa.target1 == target, aa.target2 == 'wt'),\
                                                    and_(aa.target2 == target, aa.target1 == 'wt')),\
                                                    and_(aa.gene_name == gene_name, aa.condition1 == condition, aa.fdr < .05)).all()

        if len(array_regulation) > 1: return 'error'
        elif len(array_regulation) == 0: continue

        if array_regulation[0].target1 == 'wt':
            df['value'].ix[factor,gene_name] = float(array_regulation[0].fold_change)
        else:
            df['value'].ix[factor,gene_name] = -1*float(array_regulation[0].fold_change)

    return df


def add_gene_group(name, genes):
    session = Session()

    if ome.query(GeneGroup).filter(GeneGroup.name == name).all(): return

    gene_group = GeneGroup(name)
    session.add(gene_group)

    for gene in genes:
        if isinstance(gene, basestring):
            gene = session.query(Gene).filter(or_(Gene.name == gene, Gene.locus_id == gene)).first()

        if gene: session.add(GeneGrouping(gene_group.id, gene.id))

    session.flush()
    session.commit()
    session.close()

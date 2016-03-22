"""Module to implement ORM for the experimental portion of the OME database"""

from ome.base import *

from sqlalchemy.orm import relationship, backref, column_property
from sqlalchemy import Table, MetaData, create_engine, Column, Integer, \
    String, Float, ForeignKey, ForeignKeyConstraint, PrimaryKeyConstraint, \
    select, and_, or_
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.schema import UniqueConstraint
from sqlalchemy.sql.expression import join
from sqlalchemy import func

from math import ceil


class ExpandedEnvironments:
    """
    Make carbon_source, nitrogen_source etc. queryable if the object
    is related to a condition
    """
    @hybrid_property
    def carbon_source(self):
        return self.environment.carbon_source
    @carbon_source.expression
    def carbon_source(cls):
        return select([Environment.carbon_source]).\
            where(Environment.id == cls.environment_id).label("carbon_source")
    @hybrid_property
    def nitrogen_source(self):
        return self.environment.nitrogen_source
    @nitrogen_source.expression
    def nitrogen_source(cls):
        return select([Environment.nitrogen_source]).\
            where(Environment.id == cls.environment_id).label("nitrogen_source")
    @hybrid_property
    def eacceptor(self):
        return self.environment.eacceptor
    @eacceptor.expression
    def eacceptor(cls):
        return select([Environment.eacceptor]).\
            where(Environment.id == cls.environment_id).label("eacceptor")


class Strain(Base):
    __tablename__ = 'strain'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    name = Column(String(100))

    __table_args__ = (UniqueConstraint('name'),{})

    def __repr__(self):
        return "Strain: %s" % (self.name)

    def __init__(self, name):
        self.name = name


class Environment(Base):
    __tablename__ = 'environment'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    name = Column(String(100), unique=True)
    type = Column(String(20))

    __table_args__ = (UniqueConstraint('name'),{})

    __mapper_args__ = {'polymorphic_identity':'environment',
                       'polymorphic_on': type
                       }

    def __repr__(self):
        return "Environment: type:%s, name:%s" % \
            (self.type, self.name)

    def __init__(self, name):
        self.name = name


class InSilicoEnvironment(Environment):
    __tablename__ = 'in_silico_environment'

    id = Column(Integer, ForeignKey('environment.id'), primary_key=True)
    exchanges = Column(String(100))

    __mapper_args__ = { 'polymorphic_identity': 'in_silico' }


    def __init__(self, name, exchanges):
        super(InSilicoEnvironment, self).__init__(name)
        self.exchanges = exchanges


class InVivoEnvironment(Environment):
    __tablename__ = 'in_vivo_environment'

    id = Column(Integer, ForeignKey('environment.id'), primary_key=True)

    carbon_source = Column(String(100))
    nitrogen_source = Column(String(100))
    electron_acceptor = Column(String(100))
    temperature = Column(Float)
    supplements = Column(String(100))

    __table_args__ = (UniqueConstraint('carbon_source','nitrogen_source','electron_acceptor','temperature','supplements'),{})

    __mapper_args__ = { 'polymorphic_identity': 'in_vivo' }

    def __repr__(self):
        return "Environment: C:%s, N:%s, e:%s %s" % \
            (self.carbon_source, self.nitrogen_source, self.electron_acceptor, self.supplements)


    def __init__(self, name, carbon_source, nitrogen_source, electron_acceptor, temperature, supplements):
        super(InVivoEnvironment, self).__init__(name)
        self.carbon_source = carbon_source
        self.nitrogen_source = nitrogen_source
        self.electron_acceptor = electron_acceptor
        self.temperature = temperature
        self.supplements = supplements


class Dataset(Base):
    __tablename__ = 'dataset'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    name = Column(String(100))
    type = Column(String(40))
    replicate = Column(Integer)

    strain_id = Column(Integer, ForeignKey('strain.id', ondelete='CASCADE'))
    strain = relationship("Strain")

    environment_id = Column(Integer, ForeignKey('environment.id', ondelete='CASCADE'))
    environment = relationship("Environment")

    data_source_id = Column(Integer, ForeignKey('data_source.id', ondelete='CASCADE'))
    data_source = relationship("DataSource")

    group_name = Column(String(100))


    __mapper_args__ = {'polymorphic_identity': 'dataset',
                       'polymorphic_on': type}

    __table_args__ = (UniqueConstraint('name','replicate', 'group_name'),{})

    def __repr__(self):
        return "Data Set (#%d):  %s" % \
            (self.id, self.name)


    def __init__(self, name, replicate=1, strain_id=None, environment_id=None, data_source_id=None, group_name=None):

        session = Session()
        if strain_id is None:
            strain_id = session.get_or_create(Strain, name='generic').id

        if environment_id is None:
            environment_id = session.get_or_create(Environment, name='generic').id

        if data_source_id is None:
            data_source_id = session.get_or_create(DataSource, name='generic', lab='generic', institution='generic').id
        session.close()


        self.name = name
        self.replicate = replicate
        self.strain_id = strain_id
        self.environment_id = environment_id
        self.data_source_id = data_source_id
        self.group_name = group_name


class ArrayExperiment(Dataset):
    __tablename__ = 'array_experiment'

    id = Column(Integer, ForeignKey('dataset.id', ondelete='CASCADE'), primary_key=True)
    platform = Column(String(10))

    #terrrible hack right here
    #file_name = Column(String(100))

    __mapper_args__ = { 'polymorphic_identity': 'array_experiment' }

    def __init__(self, name, replicate, strain_id, environment_id, data_source_id,\
                       platform, group_name=None):
        super(ArrayExperiment, self).__init__(name, replicate, strain_id, environment_id, data_source_id, group_name)
        self.platform = platform

    def __repr__(self):
        return "Array Experiment (#%d, %s):  %s  %d" % \
            (self.id, self.name, self.platform, self.replicate)


class RNASeqExperiment(Dataset):
    __tablename__ = 'rnaseq_experiment'

    id = Column(Integer, ForeignKey('dataset.id', ondelete='CASCADE'), primary_key=True)

    sequencing_type = Column(String(20))
    machine_id = Column(String(20))

    normalization_method = Column(String(100))
    normalization_factor = Column(Float)


    #terrrible hack right here
    #file_name = Column(String(100))

    __mapper_args__ = { 'polymorphic_identity': 'rnaseq_experiment' }

    def __repr__(self):
        return "RNASeqExperiment (#%d, %s):  %s  %s" % \
            (self.id, self.name, self.replicate, self.group_name)


    def __init__(self, name, replicate, strain_id, environment_id, data_source_id,\
                       sequencing_type, machine_id, group_name=None, normalization_method=None,\
                       normalization_factor=None):
        super(RNASeqExperiment, self).__init__(name, replicate, strain_id, environment_id, data_source_id, group_name)
        self.sequencing_type = sequencing_type
        self.machine_id = machine_id

        self.normalization_method = normalization_method
        self.normalization_factor = normalization_factor



class ChIPExperiment(Dataset):
    __tablename__ = 'chip_experiment'

    id = Column(Integer, ForeignKey('dataset.id'), primary_key=True)
    antibody = Column(String(20))
    protocol_type = Column(String(20))
    target = Column(String(20))

    normalization_method = Column(String(100))
    normalization_factor = Column(Float)


    #terrrible hacks right here
    #file_name = Column(String(100))
    #directory_path = Column(String(200))

    __mapper_args__ = { 'polymorphic_identity': 'chip_experiment' }

    def __repr__(self):
        return "ChIPExperiment (#%d, %s): %s %s %s" % \
            (self.id, self.name, self.protocol_type, self.target, self.replicate)


    def __init__(self, name, replicate, strain_id, environment_id, data_source_id,\
                       antibody, protocol_type, target, group_name=None, normalization_method=None,\
                       normalization_factor=None):

        super(ChIPExperiment, self).__init__(name, replicate, strain_id, environment_id, data_source_id, group_name)
        self.antibody = antibody
        self.protocol_type = protocol_type
        self.target = target
        self.normalization_method = normalization_method
        self.normalization_factor = normalization_factor



class AnalysisComposition(Base):
    __tablename__ = 'analysis_composition'

    analysis_id = Column(Integer, ForeignKey('analysis.id'), primary_key=True)
    dataset_id = Column(Integer, ForeignKey('dataset.id'), primary_key=True)

    __table_args__ = (UniqueConstraint('analysis_id','dataset_id'),{})

    def __init__(self, analysis_id, dataset_id):
        self.analysis_id = analysis_id
        self.dataset_id = dataset_id


class Analysis(Dataset):
    __tablename__ = 'analysis'

    id = Column(Integer, ForeignKey('dataset.id', ondelete='CASCADE'), primary_key=True)
    analysis_type = Column(String(40))
    children = relationship("Dataset", secondary="analysis_composition",\
                            primaryjoin = id == AnalysisComposition.analysis_id,\
                            backref="parent")

    __mapper_args__ = {'polymorphic_identity': 'analysis',
                       'polymorphic_on': 'analysis_type'}

    def __init__(self, name, replicate=1, strain_id=None, environment_id=None, group_name=None):
        super(Analysis, self).__init__(name, replicate, strain_id, environment_id, group_name=group_name)

    def __repr__(self):
        return "Analysis (#%d):  %s" % \
            (self.id, self.name)


class ChIPPeakAnalysis(Analysis):
    __tablename__ = 'chip_peak_analysis'

    id = Column(Integer, ForeignKey('analysis.id', ondelete='CASCADE'), primary_key=True)
    method = Column(String(30))
    parameters = Column(String(200))


    __mapper_args__ = { 'polymorphic_identity': 'chip_peak_analysis' }

    def __repr__(self):
        return "ChIP Peak Analysis (#%d, %s): %s" % \
                (self.id, self.name, self.environment)

    def __init__(self, name, replicate=1, strain_id=None, environment_id=None, method=None, parameters=None, group_name=None):
        super(ChIPPeakAnalysis, self).__init__(name, replicate, strain_id, environment_id, group_name=group_name)
        self.method = method
        self.parameters = parameters


class NormalizedExpression(Analysis):
    __tablename__ = 'normalized_expression'

    id = Column(Integer, ForeignKey('analysis.id', ondelete='CASCADE'), primary_key=True)
    norm_method = Column(String(40))
    dispersion_method = Column(String(20))
    expression_type = Column(String(20))

    __mapper_args__ = {'polymorphic_identity': 'normalized_expression'}


    def __init__(self, name, replicate=1, strain_id=None, environment_id=None, group_name=None,
                                                                               norm_method=None,
                                                                               dispersion_method=None,
                                                                               expression_type=None):
        super(NormalizedExpression, self).__init__(name, replicate, strain_id, environment_id, group_name)
        self.norm_method = norm_method
        self.dispersion_method = dispersion_method
        self.expression_type = expression_type


    def __repr__(self):
        return "Expression Data (#%d):  %s  %s" % \
            (self.id, self.name, self.expression_type)


class DifferentialExpression(Analysis):
    __tablename__ = 'differential_expression'

    id = Column(Integer, ForeignKey('analysis.id', ondelete='CASCADE'), primary_key=True)
    norm_method = Column(String(20))
    fdr = Column(Float)

    __mapper_args__ = {'polymorphic_identity': 'differential_expression'}


    def __init__(self, name, replicate=1, group_name=None, norm_method=None, fdr=None):
        super(DifferentialExpression, self).__init__(name, replicate, group_name=group_name)
        self.norm_method = norm_method
        self.fdr = fdr


    def __repr__(self):
        return "Differential Expression (#%d): %s" % \
            (self.id, self.name)


class GenomeData(Base):
    __tablename__ = 'genome_data'

    dataset_id = Column(Integer, ForeignKey('dataset.id', ondelete="CASCADE"), primary_key=True)
    dataset = relationship('Dataset')
    genome_region_id = Column(Integer, ForeignKey('genome_region.id'), primary_key=True)
    genome_region = relationship('GenomeRegion', backref='data')
    value = Column(Float)
    type = Column(String(20))

    __table_args__ = (UniqueConstraint('dataset_id','genome_region_id'),{})

    @hybrid_property
    def all_data(self):
        return [x['value'] for x in query_genome_data([self.dataset_id], self.genome_region.leftpos, self.genome_region.rightpos)]


    __mapper_args__ = {'polymorphic_identity': 'genome_data',
                       'polymorphic_on': type}

    def __repr__(self):
        return "%s: %5.2f -- %s" % \
            (self.genome_region, self.value, self.dataset.name, )


    def __init__(self, dataset_id, genome_region_id, value):
        self.dataset_id = dataset_id
        self.genome_region_id = genome_region_id
        self.value = value


class DiffExpData(GenomeData):
    __tablename__ = 'diff_exp_data'

    dataset_id = Column(Integer, primary_key=True)
    genome_region_id = Column(Integer, primary_key=True)

    diff_exp_analysis = relationship('DifferentialExpression')

    pval = Column(Float)

    __table_args__ = (ForeignKeyConstraint(['dataset_id','genome_region_id'],\
                                           ['genome_data.dataset_id', 'genome_data.genome_region_id'], ondelete='CASCADE'),\
                      UniqueConstraint('dataset_id','genome_region_id'),{})

    __mapper_args__ = { 'polymorphic_identity': 'diff_exp_data' }

    def __repr__(self):
        return "Diff Exp Data: %s %5.2f %5.2f %s" % \
            (self.genome_region, self.value, self.pval, self.diff_exp_analysis.name)

    def __init__(self, dataset_id, genome_region_id, value, pval):
        super(DiffExpData, self).__init__(dataset_id, genome_region_id, value)
        self.pval = pval


class ChIPPeakData(GenomeData):
    __tablename__ = 'chip_peak_data'

    dataset_id = Column(Integer, primary_key=True)
    genome_region_id = Column(Integer, primary_key=True)

    peak_analysis = relationship('Analysis')
    eventpos = Column(Integer)
    pval = Column(Float)

    @hybrid_property
    def grouped_eventpos(self):
        return ceil(self.eventpos/400) * 400

    @grouped_eventpos.expression
    def grouped_eventpos(cls):
        return func.ceil(ChIPPeakData.eventpos/400) * 400

    __table_args__ = (ForeignKeyConstraint(['dataset_id','genome_region_id'],\
                                           ['genome_data.dataset_id', 'genome_data.genome_region_id'], ondelete='CASCADE'),\
                      UniqueConstraint('dataset_id','genome_region_id'),{})

    __mapper_args__ = { 'polymorphic_identity': 'chip_peak_data' }

    def __repr__(self):
        return "ChIP Peak: %d-%d %5.2f %s" % \
            (self.genome_region.leftpos, self.genome_region.rightpos, self.value, self.peak_analysis.name)

    def __init__(self, dataset_id, genome_region_id, value, eventpos, pval):
        super(ChIPPeakData, self).__init__(dataset_id, genome_region_id, value)
        self.pval = pval
        self.eventpos = eventpos



class RegulatoryNetwork(Base):
    __tablename__ = 'regulatory_network'

    reg_gene_id = Column(Integer, ForeignKey('genome_region.id', ondelete="CASCADE"), primary_key=True)
    regd_gene_id = Column(Integer, ForeignKey('genome_region.id', ondelete="CASCADE"), primary_key=True)
    direction = Column(String(3))
    evidence = Column(String(40))

    reg_gene = relationship('GenomeRegion',
                             primaryjoin = reg_gene_id == GenomeRegion.id, backref='regulates')

    regd_gene = relationship('GenomeRegion',
                              primaryjoin = regd_gene_id == GenomeRegion.id, backref='regulation')

    __table_args__ = (UniqueConstraint('reg_gene_id','regd_gene_id'),{})


    def __repr__(self):
        return "Regulation: %s --> %s --> %s, %s" % \
            (self.reg_gene.name, self.direction, self.regd_gene.name, self.evidence)


    def __init__(self, reg_gene_id, regd_gene_id, direction, evidence):
        self.reg_gene_id = reg_gene_id
        self.regd_gene_id = regd_gene_id
        self.direction = direction
        self.evidence = evidence



ome = Session()
from ome.components import Gene,TU,TUGenes


grouped = ome.query(func.distinct(ChIPPeakData.dataset_id).label('dataset_id'),\
                    func.max(ChIPPeakData.value).label('value'),\
                    ChIPPeakData.grouped_eventpos.label('eventpos'),\
                    GenomeRegion.strand.label('strand')).\
                    join(GenomeRegion).\
                    group_by(ChIPPeakData.dataset_id, ChIPPeakData.grouped_eventpos, GenomeRegion.strand).subquery()




gene_expression_data = ome.query(
          Gene.id.label('gene_id'),
          Gene.name.label('gene_name'),
          Gene.bigg_id.label('bigg_id'),
          Strain.id.label('strain_id'),
          Strain.name.label('strain'),
          InVivoEnvironment.id.label('environment_id'),
          InVivoEnvironment.carbon_source.label('carbon_source'),
          InVivoEnvironment.nitrogen_source.label('nitrogen_source'),
          InVivoEnvironment.electron_acceptor.label('electron_acceptor'),
          Dataset.type.label('dataset_type'),
          Dataset.group_name.label('group_name'),
          func.max(Dataset.id).label('max_dataset_id'),
          func.avg(GenomeData.value).label('value'),
          func.stddev_pop(GenomeData.value).label('stddev')).\
    join(GenomeData, Dataset, Strain, InVivoEnvironment).\
    group_by(Gene.id, Dataset.group_name, Strain.id, Dataset.type, InVivoEnvironment.id,
             Gene.bigg_id, Gene.name, Strain.name, InVivoEnvironment.carbon_source,
                                                    InVivoEnvironment.nitrogen_source,
                                                    InVivoEnvironment.electron_acceptor).order_by(Dataset.type, InVivoEnvironment.id).subquery()



class GeneExpressionData(Base):
    __table__ = gene_expression_data

    __mapper_args__ = {
        'primary_key':[gene_expression_data.c.gene_id, gene_expression_data.c.max_dataset_id]
    }


    def __repr__(self):
        return "Gene: (%s, %s), Value: %5.2f, std:%5.2f, Condition: %s, %s, %s, Strain: %s, %s %s" % \
            (self.locus_id, self.gene_name, self.value, self.stddev, self.carbon_source,
             self.nitrogen_source, self.electron_acceptor, self.strain, self.dataset_type, self.group_name)



AnalysisComposition2 = aliased(AnalysisComposition)
AnalysisComposition3 = aliased(AnalysisComposition)
NormalizedExpression2 = aliased(NormalizedExpression)


chip_peak = ome.query(ChIPPeakAnalysis.id.label('id'),
                      ChIPPeakAnalysis.environment_id.label('environment_id'),
                      ChIPPeakAnalysis.group_name.label('chip_peak_group'),
                      ChIPExperiment.target.label('target'),
                      ChIPExperiment.antibody.label('antibody'),
                      Strain.id.label('strain_id'),
                      Strain.name.label('strain')).\
                      join(AnalysisComposition, ChIPPeakAnalysis.id == AnalysisComposition.analysis_id).\
                      join(ChIPExperiment, ChIPExperiment.id == AnalysisComposition.dataset_id).\
                      join(Strain, ChIPExperiment.strain_id == Strain.id).\
                      group_by(ChIPPeakAnalysis.id, ChIPPeakAnalysis.environment_id, ChIPPeakAnalysis.group_name,
                               ChIPExperiment.target, ChIPExperiment.antibody, Strain.id).subquery()


Strain2 = aliased(Strain)



diff_exp = ome.query(DifferentialExpression.id.label('id'),
                     Strain.name.label('strain1'),
                     Strain2.name.label('strain2'),
                     InVivoEnvironment.id.label('environment_id'),
                     InVivoEnvironment.carbon_source.label('carbon_source'),
                     InVivoEnvironment.nitrogen_source.label('nitrogen_source'),
                     InVivoEnvironment.electron_acceptor.label('electron_acceptor'),
                     InVivoEnvironment.supplements.label('supplements'),
                     NormalizedExpression.expression_type.label('expression_type')).\
                    join(AnalysisComposition2, DifferentialExpression.id == AnalysisComposition2.analysis_id).\
                    join(AnalysisComposition3, DifferentialExpression.id == AnalysisComposition3.analysis_id).\
                    join(NormalizedExpression, AnalysisComposition2.dataset_id == NormalizedExpression.id).\
                    join(NormalizedExpression2, AnalysisComposition3.dataset_id == NormalizedExpression2.id).\
                    join(Strain, Strain.id == NormalizedExpression.strain_id).\
                    join(Strain2, Strain2.id == NormalizedExpression2.strain_id).\
                    join(InVivoEnvironment, InVivoEnvironment.id == NormalizedExpression.environment_id).\
                    filter(and_(NormalizedExpression.environment_id == NormalizedExpression2.environment_id,
                                NormalizedExpression.strain_id != NormalizedExpression2.strain_id)).subquery()


diff_exp_chip_peak = ome.query(diff_exp.c.id.label('diff_exp_id'),
                               chip_peak.c.id.label('chip_peak_id'),
                               chip_peak.c.chip_peak_group,
                               chip_peak.c.target,
                               chip_peak.c.antibody,
                               diff_exp.c.strain1,
                               diff_exp.c.strain2,
                               diff_exp.c.environment_id,
                               diff_exp.c.carbon_source,
                               diff_exp.c.nitrogen_source,
                               diff_exp.c.electron_acceptor,
                               diff_exp.c.supplements,
                               diff_exp.c.expression_type).\
                               join(chip_peak, and_(chip_peak.c.environment_id == diff_exp.c.environment_id,
                                                    func.lower(chip_peak.c.target) == func.substr(diff_exp.c.strain1, 7, func.length(diff_exp.c.strain1)),
                                                    diff_exp.c.strain2 == 'wt', chip_peak.c.strain == 'wt')).subquery()





GenomeRegion2 = aliased(GenomeRegion)


chip_peak_gene = ome.query(ChIPPeakData.value.label('peak_value'),
                           ChIPPeakData.genome_region_id.label('peak_genome_region_id'),
                           GenomeRegion.leftpos.label('leftpos'),
                           GenomeRegion.rightpos.label('rightpos'),
                           GenomeRegion.strand.label('strand'),
                           chip_peak.c.target.label('target'),
                           chip_peak.c.strain_id.label('strain_id'),
                           chip_peak.c.strain.label('strain'),
                           chip_peak.c.antibody.label('antibody'),
                           chip_peak.c.environment_id.label('environment_id'),
                           InVivoEnvironment.carbon_source.label('carbon_source'),
                           InVivoEnvironment.nitrogen_source.label('nitrogen_source'),
                           InVivoEnvironment.electron_acceptor.label('electron_acceptor'),
                           InVivoEnvironment.supplements.label('supplements'),
                           TU.name.label('tu_name'),
                           Gene.id.label('gene_id'),
                           Gene.name.label('name'),
                           GenomeRegion2.bigg_id.label('gene_bigg_id'),
                           ChIPPeakData.dataset_id.label('chip_peak_id')).\
                    join(GenomeRegionMap, GenomeRegionMap.genome_region_id_1 == ChIPPeakData.genome_region_id).\
                    join(TU, GenomeRegionMap.genome_region_id_2 == TU.genome_region_id).\
                    join(TUGenes, TUGenes.tu_id == TU.id).\
                    join(Gene, Gene.id == TUGenes.gene_id).\
                    join(GenomeRegion, GenomeRegion.id == ChIPPeakData.genome_region_id).\
                    join(chip_peak, chip_peak.c.id == ChIPPeakData.dataset_id).\
                    join(InVivoEnvironment, InVivoEnvironment.id == chip_peak.c.environment_id).\
                    join(GenomeRegion2, GenomeRegion2.id == TUGenes.gene_id).subquery()




chip_peak_gene_expression = ome.query(ChIPPeakData.value.label('peak_value'),
                 ChIPPeakData.genome_region_id.label('peak_genome_region_id'),
                 GenomeRegion.leftpos.label('leftpos'),
                 GenomeRegion.rightpos.label('rightpos'),
                 GenomeRegion.strand.label('strand'),
                 diff_exp_chip_peak.c.target.label('target'),
                 diff_exp_chip_peak.c.strain1.label('strain1'),
                 diff_exp_chip_peak.c.strain2.label('strain2'),
                 diff_exp_chip_peak.c.antibody.label('antibody'),
                 diff_exp_chip_peak.c.environment_id.label('environment_id'),
                 diff_exp_chip_peak.c.carbon_source.label('carbon_source'),
                 diff_exp_chip_peak.c.nitrogen_source.label('nitrogen_source'),
                 diff_exp_chip_peak.c.electron_acceptor.label('electron_acceptor'),
                 TU.name.label('tu_name'),
                 Gene.id.label('gene_id'),
                 Gene.name.label('name'),
                 GenomeRegion2.bigg_id.label('gene_bigg_id'),
                 DiffExpData.type.label('dataset_type'),
                 DiffExpData.genome_region_id.label('gene_genome_region_id'),
                 DiffExpData.value.label('value'),
                 DiffExpData.pval.label('pval'),
                 ChIPPeakData.dataset_id.label('chip_peak_id'),
                 ChIPPeakAnalysis.group_name.label('chip_peak_group'),
                 DiffExpData.dataset_id.label('diff_exp_id')).\
                    join(GenomeRegionMap, GenomeRegionMap.genome_region_id_1 == ChIPPeakData.genome_region_id).\
                    join(TU, GenomeRegionMap.genome_region_id_2 == TU.genome_region_id).\
                    join(TUGenes, TUGenes.tu_id == TU.id).\
                    join(DiffExpData, DiffExpData.genome_region_id == TUGenes.gene_id).\
                    join(diff_exp_chip_peak, and_(DiffExpData.dataset_id == diff_exp_chip_peak.c.diff_exp_id,
                                                  ChIPPeakData.dataset_id == diff_exp_chip_peak.c.chip_peak_id)).\
                    join(Gene, Gene.id == TUGenes.gene_id).\
                    join(GenomeRegion, GenomeRegion.id == ChIPPeakData.genome_region_id).\
                    join(GenomeRegion2, GenomeRegion2.id == TUGenes.gene_id).\
                    join(ChIPPeakAnalysis, ChIPPeakAnalysis.id == ChIPPeakData.dataset_id).subquery()





class ChIPPeakGene(Base):
    __table__ = chip_peak_gene

    def __repr__(self):
        return "TF: %s, Gene: (%s, %s), Condition: %s, %s, %s Peak: %d-%d value:%5.2f" % \
            (self.target, self.gene_name, self.locus_id,
             self.carbon_source, self.nitrogen_source, self.electron_acceptor, self.leftpos, self.rightpos, self.peak_value)



class ChIPPeakGeneExpression(Base):
    __table__ = chip_peak_gene_expression

    def __repr__(self):
        return "TF: %s, Gene: (%s, %s), %5.2f, %5.2f %s-->%s Condition: %s, %s, %s Peak: %d-%d value:%5.2f" % \
            (self.target, self.gene_name, self.locus_id, self.value, self.pval, self.strain1, self.strain2,
             self.carbon_source, self.nitrogen_source, self.electron_acceptor, self.leftpos, self.rightpos, self.peak_value)

InVivoEnvironment2 = aliased(InVivoEnvironment)

differential_gene_expression_data = ome.query(DiffExpData.value.label('value'),
                                           DiffExpData.pval.label('pval'),
                                           DiffExpData.dataset_id.label('diff_exp_id'),
                                           NormalizedExpression.id.label('expression_id'),
                                           NormalizedExpression.expression_type.label('dataset_type'),
                                           Gene.id.label('gene_id'),
                                           Gene.bigg_id.label('bigg_id'),
                                           Gene.name.label('gene_name'),
                                           Strain.name.label('strain1'),
                                           Strain2.name.label('strain2'),
                                           InVivoEnvironment.id.label('environment_id_1'),
                                           InVivoEnvironment.carbon_source.label('carbon_source1'),
                                           InVivoEnvironment.nitrogen_source.label('nitrogen_source1'),
                                           InVivoEnvironment.electron_acceptor.label('electron_acceptor1'),
                                           InVivoEnvironment.supplements.label('supplements1'),
                                           InVivoEnvironment2.id.label('environment_id_2'),
                                           InVivoEnvironment2.carbon_source.label('carbon_source2'),
                                           InVivoEnvironment2.nitrogen_source.label('nitrogen_source2'),
                                           InVivoEnvironment2.electron_acceptor.label('electron_acceptor2'),
                                           InVivoEnvironment2.supplements.label('supplements2'),).\
                                        join(Gene).\
                                        join(DifferentialExpression).\
                                        join(AnalysisComposition, DifferentialExpression.id == AnalysisComposition.analysis_id).\
                       					join(AnalysisComposition2, DifferentialExpression.id == AnalysisComposition2.analysis_id).\
                       					join(NormalizedExpression, AnalysisComposition.dataset_id == NormalizedExpression.id).\
                       					join(NormalizedExpression2, AnalysisComposition2.dataset_id == NormalizedExpression2.id).\
                       					join(Strain, NormalizedExpression.strain_id == Strain.id).\
                       					join(Strain2, NormalizedExpression2.strain_id == Strain2.id).\
                       					join(InVivoEnvironment, NormalizedExpression.environment_id == InVivoEnvironment.id).\
                       					join(InVivoEnvironment2, NormalizedExpression2.environment_id == InVivoEnvironment2.id).\
                        				 filter(or_(InVivoEnvironment.id != InVivoEnvironment2.id, Strain.id != Strain2.id)).subquery()


class DifferentialGeneExpressionData(Base):
    __table__ = differential_gene_expression_data

    def diff_name(self):
        args = dict.fromkeys(['strain','carbon_source','nitrogen_source','electron_acceptor'], '')

        if self.strain1 == self.strain2: args['strain'] = self.strain1
        else: args['strain'] = self.strain1+'/'+self.strain2

        if self.carbon_source1 == self.carbon_source2: args['carbon_source'] = self.carbon_source1
        else: args['carbon_source'] = self.carbon_source1+'/'+self.carbon_source2

        if self.nitrogen_source1 == self.nitrogen_source2: args['nitrogen_source'] = self.nitrogen_source1
        else: args['nitrogen_source'] = self.nitrogen_source1+'/'+self.nitrogen_source2

        if self.electron_acceptor1 == self.electron_acceptor2: args['electron_acceptor'] = self.electron_acceptor1
        else: args['electron_acceptor'] = self.electron_acceptor1+'/'+self.electron_acceptor2

        if self.supplements1 == self.supplements2: args['supplements'] = self.supplements1
        else: args['supplements'] = self.supplements1+'/'+self.supplements2

        if args['supplements'] == '':
            return '_'.join([self.dataset_type, args['strain'], args['carbon_source'], args['nitrogen_source'], args['electron_acceptor']])
        else:
            return '_'.join([self.dataset_type, args['strain'], args['carbon_source'], args['nitrogen_source'], args['electron_acceptor'], args['supplements']])

    def __repr__(self):
        args = dict.fromkeys(['strain','carbon_source','nitrogen_source','electron_acceptor'], '')
        if self.strain1 == self.strain2: args['strain'] = self.strain1
        else: args['strain'] = self.strain1+'/'+self.strain2

        if self.carbon_source1 == self.carbon_source2: args['carbon_source'] = self.carbon_source1
        else: args['carbon_source'] = self.carbon_source1+'/'+self.carbon_source2

        if self.nitrogen_source1 == self.nitrogen_source2: args['nitrogen_source'] = self.nitrogen_source1
        else: args['nitrogen_source'] = self.nitrogen_source1+'/'+self.nitrogen_source2

        if self.electron_acceptor1 == self.electron_acceptor2: args['electron_acceptor'] = self.electron_acceptor1
        else: args['electron_acceptor'] = self.electron_acceptor1+'/'+self.electron_acceptor2

        return "Gene: (%s, %s), %s, %s, %s, %s, Fold Change: %5.2f, FDR: %5.2f" % \
				    	  (self.locus_id, self.gene_name, args['strain'], args['carbon_source'],
					  								  args['nitrogen_source'], args['electron_acceptor'],
					  								  self.value, self.pval)


def query_genome_data(dataset_ids,leftpos=0,rightpos=1000,strand=['+','-'],value=0):
    genome_data = omics_database.genome_data
    return genome_data.find({"$and":
                            [{"leftpos" : {"$gte": leftpos}},
                             {"leftpos": {"$lte": rightpos}},
                             {"strand": {"$in" : strand }},
                             {"dataset_id": {"$in" : dataset_ids}}
                         ]})


def query_reaction_data():
    return


def query_metabolite_data():
    return

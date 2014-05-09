"""Module to implement ORM for the experimental portion of the OME database"""

from PrototypeDB.orm.base import *

from sqlalchemy.orm import relationship, backref, column_property
from sqlalchemy import Table, MetaData, create_engine, Column, Integer, \
    String, Float, ForeignKey, ForeignKeyConstraint, select
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.schema import UniqueConstraint
from sqlalchemy.sql.expression import join
from pymongo import ASCENDING, DESCENDING

import simplejson as json


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
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(100))
    
    __table_args__ = (UniqueConstraint('name'),{})
    
    def __repr__(self):
        return "Strain: %s" % (self.name)
    
    def __repr__dict__(self):
        return {"name":self.name,"wid":self.id,"values":{}}
    
    def __repr__json__(self):
        return json.dumps(self.__repr__dict__())
    
    def __init__(self, name):
        self.name = name
        

class Environment(Base):
    __tablename__ = 'environment'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(100), unique=True)
    type = Column(String(20))
    
    __table_args__ = (UniqueConstraint('name'),{})

    __mapper_args__ = {'polymorphic_identity':'environment',
                       'polymorphic_on': type
                       }
    
    def __repr__(self):
        return "Environment: type:%s, name:%s" % \
            (self.type, self.name)
    
    def __repr__dict__(self):
        return {"name":self.name,"wid":self.id, "type":self.type, "values":{}}
    
    def __repr__json__(self):
        return json.dumps(self.__repr__dict__())    
    
    def __init__(self, name):
        self.name = name
    

class InSilicoEnvironment(Environment):
    __tablename__ = 'in_silico_environment'
    
    id = Column(Integer, ForeignKey('environment.id'), primary_key=True)
    exchanges = Column(String(100))
    
    __mapper_args__ = { 'polymorphic_identity': 'in_silico' }
    
    def __repr__dict(self):
        environment = Environment.__repr__dict__(self)
        environment['values'] = {}
        return environment
    
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
    
    __table_args__ = (UniqueConstraint('carbon_source','nitrogen_source','electron_acceptor','temperature'),{})
    
    __mapper_args__ = { 'polymorphic_identity': 'in_vivo' }    
    
    def __repr__(self):
        return "Environment: C:%s, N:%s, e:%s" % \
            (self.carbon_source, self.nitrogen_source, self.electron_acceptor)
                
    def __repr__dict__(self):
        environment = Environment.__repr__dict__(self)
        environment['values'] = {"carbon_source":self.carbon_source,\
                                 "nitrogen_source":self.nitrogen_source,\
                                 "electron_acceptor":self.electron_acceptor,
                                 "temperature":self.temperature}
        return environment
    
    def __init__(self, name, carbon_source, nitrogen_source, electron_acceptor, temperature):
        super(InVivoEnvironment, self).__init__(name)
        self.carbon_source = carbon_source
        self.nitrogen_source = nitrogen_source
        self.electron_acceptor = electron_acceptor
        self.temperature = temperature


class Protocol(Base):
    __tablename__ = 'protocol'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(100))
    location = Column(String(100))
    
    __table_args__ = (UniqueConstraint('name', name='name'),{})
    
    def __repr__(self):
        return "Protocol (#%d, %s):  %s" % \
            (self.id, self.name, self.location)
    

class DataSet(Base):
    __tablename__ = 'data_set'
    
    id = Column(Integer, primary_key=True, autoincrement=True)   
    name = Column(String(100))
    type = Column(String(40))
    replicate = Column(Integer)
    
    strain_id = Column(Integer, ForeignKey('strain.id'))
    strain = relationship("Strain")
    
    environment_id = Column(Integer, ForeignKey('environment.id'))
    environment = relationship("Environment")
    
    data_source_id = Column(Integer, ForeignKey('data_source.id'))
    data_source = relationship("DataSource", cascade='all,delete')
    
    __mapper_args__ = {'polymorphic_identity': 'data_set',
                       'polymorphic_on': type}
    
    __table_args__ = (UniqueConstraint('name'),{})
    
    def __repr__(self):
        return "Data Set (#%d):  %s" % \
            (self.id, self.name)
    
    def __repr__dict__(self):
        experiment = {"name":self.name,"id":self.id, "type":self.type, "values":{}}
        
        dataset['values']['strain'] = self.strain.__repr__dict__()
        dataset['values']['data_source'] = self.data_source.__repr__dict__()
        dataset['values']['environment'] = self.environment.__repr__dict__()
        
        return dataset
    
    def __repr__json__(self):
        return json.dumps(self.__repr__dict__())

    def __init__(self, name, replicate=1, strain=None, environment=None, data_source=None):
        
        if strain is None:
            strain = get_or_create(Session(), Strain, name='generic')
        
        if environment is None:
            environment = get_or_create(Session(), Environment, name='generic')
        
        if data_source is None:
            data_source = get_or_create(Session(), DataSource, name='generic', lab='generic', institution='generic')
        
        self.name = name
        self.replicate = replicate
        self.strain_id = strain.id
        self.environment_id = environment.id
        self.data_source_id = data_source.id
        

class ArrayExperiment(DataSet):
    __tablename__ = 'array_experiment'
    
    id = Column(Integer, ForeignKey('data_set.id'), primary_key=True)
    platform = Column(String(10))
    
    #terrrible hack right here
    file_name = Column(String(100))
    
    __mapper_args__ = { 'polymorphic_identity': 'array_experiment' }
    
    def __init__(self, name, replicate, strain, environment, data_source,\
                       platform, file_name):
        super(ArrayExperiment, self).__init__(name, replicate, strain, environment, data_source)
        self.platform = platform
        self.file_name = file_name
        
    def __repr__(self):
        return "Array Experiment (#%d, %s):  %s  %d" % \
            (self.id, self.name, self.platform, self.replicate)
    
    def __repr__dict__(self):
        data_set = DataSet.__repr__dict__(self)
                                         
        data_set['values']['platform'] = {"name":self.platform,"values":{}}
        data_set['values']['replicate'] = {"name":self.replicate,"values":{}}
                                                            
        return data_set
    

class RNASeqExperiment(DataSet):
    __tablename__ = 'rna_seq_experiment'
    
    id = Column(Integer, ForeignKey('data_set.id'), primary_key=True)

    sequencing_type = Column(String(20))
    machine_ID = Column(String(20))
    
    #terrrible hack right here
    file_name = Column(String(100))
    
    __mapper_args__ = { 'polymorphic_identity': 'rna_seq_experiment' }          
    
    def __repr__(self):
        return "RNASeqExperiment (#%d, %s):  %s" % \
            (self.id, self.name, self.replicate)

    def __repr__dict__(self):
        data_set = DataSet.__repr__dict__(self)
                                         
        data_set['values']['sequencing_type'] = {"name":self.sequencing_type,"values":{}}
        data_set['values']['machine_id'] = {"name":self.machine_id,"values":{}}
        data_set['values']['replicate'] = {"name":self.replicate,"values":{}}               
        
        return data_set
    
    def __init__(self, name, replicate, strain, environment, data_source,\
                       sequencing_type, machine_id, file_name):
        super(RNASeqExperiment, self).__init__(name, replicate, strain, environment, data_source)
        self.sequencing_type = sequencing_type
        self.machine_id = machine_id
        self.file_name = file_name
    
    
class ChIPExperiment(DataSet):
    __tablename__ = 'chip_experiment'
    
    id = Column(Integer, ForeignKey('data_set.id'), primary_key=True)

    antibody = Column(String(20))
    protocol_type = Column(String(20))
    target = Column(String(20))
    
    #terrrible hack right here
    file_name = Column(String(100))
    
    __mapper_args__ = { 'polymorphic_identity': 'chip_experiment' }
    
    def __repr__(self):
        return "ChIPExperiment (#%d, %s): %s %s %s" % \
            (self.id, self.name, self.protocol_type, self.target, self.replicate)

    def __repr__dict__(self):
        data_set = Dataset.__repr__dict__(self)
                                         
        data_set['values']['antibody'] = {"name":self.antibody,"values":{}}
        data_set['values']['protocol_type'] = {"name":self.protocol_type,"values":{}}
        data_set['values']['target'] = {"name":self.target,"values":{}} 
        data_set['values']['replicate'] = {"name":self.replicate,"values":{}}               
        
        return dataset
    
    def __init__(self, name, replicate, strain, environment, data_source,\
                       antibody, protocol_type, target, file_name):
        super(ChIPExperiment, self).__init__(name, replicate, strain, environment, data_source)
        self.antibody = antibody
        self.protocol_type = protocol_type
        self.target = target
        self.file_name = file_name


class AnalysisComposition(Base):
    __tablename__ = 'analysis_composition'
    
    analysis_id = Column(Integer, ForeignKey('analysis.id'), primary_key=True)
    data_set_id = Column(Integer, ForeignKey('data_set.id'), primary_key=True)

    __table_args__ = (UniqueConstraint('analysis_id','data_set_id'),{})

    def __init__(self, analysis_id, data_set_id):
        self.analysis_id = analysis_id
        self.data_set_id = data_set_id
        

class Analysis(DataSet):
    __tablename__ = 'analysis'
    
    id = Column(Integer, ForeignKey('data_set.id'), primary_key=True)
    type = Column(String(40))
    children = relationship("DataSet", secondary="analysis_composition",\
                            primaryjoin = id == AnalysisComposition.analysis_id,\
                            backref="parent")
    
    __mapper_args__ = {'polymorphic_identity': 'analysis',
                       'polymorphic_on': 'type'}
    
    def __init__(self, name, replicate=1, strain=None, environment=None):
        super(Analysis, self).__init__(name, replicate, strain, environment)
    
    def __repr__(self):
        return "Analysis (#%d):  %s" % \
            (self.id, self.name)


class ChIPPeak(Analysis):
    __tablename__ = 'chip_peak_analysis'
    
    id = Column(Integer, ForeignKey('analysis.id'), primary_key=True)
    method = Column(String(30))
    parameters = Column(String(100))

    __mapper_args__ = { 'polymorphic_identity': 'chip_peak_analysis' }

    def __repr__(self):
        return "ChIP Peak Analysis (#%d, %s): %s %s" % \
                (self.id, self.name, self.method, self.parameters)

    def __init__(self, name, replicate=1, strain=None, environment=None, method=None, parameters=None):
        super(ChIPPeak, self).__init__(name, replicate, strain, environment)
        self.method = method
        self.parameters = parameters


class NormalizedExpression(Analysis):
    __tablename__ = 'normalized_expression'
    
    id = Column(Integer, ForeignKey('analysis.id'), primary_key=True)
    norm_method = Column(String(40))
    dispersion_method = Column(String(20))
    
    __mapper_args__ = {'polymorphic_identity': 'normalized_expression'}

    
    def __init__(self, name, replicate=1, strain=None, environment=None, norm_method=None, dispersion_method=None):
        super(NormalizedExpression, self).__init__(name, replicate, strain, environment)
        self.norm_method = norm_method
        self.dispersion_method = dispersion_method
        
    def __repr__(self):
        return "Expression Data (#%d):  %s" % \
            (self.id, self.name)
    

class DifferentialExpression(Analysis):
    __tablename__ = 'differential_expression'
    
    id = Column(Integer, ForeignKey('analysis.id'), primary_key=True)
    norm_method = Column(String(20))
    fdr = Column(Float)       
    
    __mapper_args__ = {'polymorphic_identity': 'differential_expression'}


    def __init__(self, name, norm_method, fdr):
        super(DifferentialExpression, self).__init__(name)
        self.norm_method = norm_method
        self.fdr = fdr
        
        
    def __repr__(self):
        return "Differential Expression (#%d): %s" % \
            (self.id, self.name)


class GenomeData(Base):
    __tablename__ = 'genome_data'
    
    data_set_id = Column(Integer, ForeignKey('data_set.id'), primary_key=True)
    data_set = relationship('DataSet')
    genome_region_id = Column(Integer, ForeignKey('genome_region.id'), primary_key=True)
    genome_region = relationship('GenomeRegion', backref='data')
    value = Column(Float)
    type = Column(String(20))
    
    __table_args__ = (UniqueConstraint('data_set_id','genome_region_id'),{})
    
    @hybrid_property
    def all_data(self):
        return [x['value'] for x in query_genome_data([self.data_set_id], self.genome_region.leftpos, self.genome_region.rightpos)]
    

    __mapper_args__ = {'polymorphic_identity': 'genome_data',
                       'polymorphic_on': type}

    def __repr__(self):
        return "%s: %5.2f -- %s" % \
            (self.genome_region, self.value, self.data_set.name, )
        

    def __init__(self, data_set_id, genome_region_id, value):
        self.data_set_id = data_set_id
        self.genome_region_id = genome_region_id
        self.value = value


class DiffExpData(GenomeData):
    __tablename__ = 'diff_exp_data'
    
    data_set_id = Column(Integer, primary_key=True)
    diff_exp_analysis = relationship('DifferentialExpression')
    genome_region_id = Column(Integer, primary_key=True)
    
    pval = Column(Float)
    
    __table_args__ = (ForeignKeyConstraint(['data_set_id','genome_region_id'],\
                                           ['genome_data.data_set_id', 'genome_data.genome_region_id']),\
                      UniqueConstraint('data_set_id','genome_region_id'),{})
    
    __mapper_args__ = { 'polymorphic_identity': 'diff_exp_data' }

    def __repr__(self):
        return "Diff Exp Data: %s %5.2f %5.2f %s" % \
            (self.genome_region, self.value, self.pval, self.diff_exp_analysis.name)
    
    def __init__(self, data_set_id, genome_region_id, value, pval):
        super(DiffExpData, self).__init__(data_set_id, genome_region_id, value)
        self.pval = pval


class ChIPPeakData(GenomeData):
    __tablename__ = 'chip_peak_data'
    
    data_set_id = Column(Integer, primary_key=True)
    genome_region_id = Column(Integer, primary_key=True)
    peak_analysis = relationship('Analysis')
    eventpos = Column(Integer)
    pval = Column(Float)
    
    __table_args__ = (ForeignKeyConstraint(['data_set_id','genome_region_id'],\
                                           ['genome_data.data_set_id', 'genome_data.genome_region_id']),\
                      UniqueConstraint('data_set_id','genome_region_id'),{})
    
    __mapper_args__ = { 'polymorphic_identity': 'chip_peak_data' }

    def __repr__(self):
        return "ChIP Peak: %d-%d %5.2f %s" % \
            (self.genome_region.leftpos, self.genome_region.rightpos, self.value, self.peak_analysis.name)
    
    def __init__(self, data_set_id, genome_region_id, value, eventpos, pval):
        super(ChIPPeakData, self).__init__(data_set_id, genome_region_id, value)
        self.pval = pval
        self.eventpos = eventpos
    

def _load_data(collection,data):
    collection.insert(data)
    
    
def load_genome_data(file_path, data_set_id, bulk_file_load=False, loading_cutoff=0):
    genome_data = omics_database.genome_data
    if file_path[-3:] == 'gff':
        entries = []
        for cntr,line in enumerate(open(file_path,'r').readlines()):
            
            if line[0] == '#': continue
            data = line.rstrip('\n').split('\t')
            
            #filter out low read count data
            if float(data[5]) < loading_cutoff: continue
            
            entries.append({
                "leftpos": int(data[3]),
                "rightpos": int(data[4]),
                "value": float(data[5]),
                "strand": data[6],
                "data_set_id": data_set_id})
            
            if cntr%10000 == 0:
                genome_data.insert(entries)
                entries = []
                
    if not bulk_file_load: 
        genome_data.create_index([("data_set_id", ASCENDING), ("leftpos", ASCENDING)])
   
   
def query_genome_data(data_set_ids,leftpos=0,rightpos=1000,strand=['+','-'],value=0):
    genome_data = omics_database.genome_data
    return genome_data.find({"$and": 
                            [{"leftpos" : {"$gte": leftpos}}, 
                             {"leftpos": {"$lte": rightpos}}, 
                             {"strand": {"$in" : strand }},
                             {"data_set_id": {"$in" : data_set_ids}}
                         ]})
    
    
def query_reaction_data():
    return
    
    
def query_metabolite_data():
    return
      


    
    
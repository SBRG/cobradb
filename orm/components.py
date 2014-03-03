from PrototypeDB.orm.base import *

from sqlalchemy.orm import relationship
from sqlalchemy import Table, MetaData, create_engine, Column, Integer, \
    String, Float, ForeignKey, select
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.schema import UniqueConstraint,PrimaryKeyConstraint

import simplejson as json



class GenomeRegion(Base):
    __tablename__ = 'genome_region'
    __table_args__ = (UniqueConstraint('leftpos','rightpos','strand'),)

    
    id = Column(Integer, primary_key=True)
    leftpos = Column(Integer, nullable=False)
    rightpos = Column(Integer, nullable=False)
    strand = Column(String(1), nullable=False)
    
    
    def __repr__(self):
        return "GenomeRegion (#%d): %d-%d (%s strand)" % \
                (self.id, self.leftpos, self.rightpos, self.strand)
                
    def __init__(self, leftpos, rightpos, strand):
        self.leftpos = leftpos
        self.rightpos = rightpos
        self.strand = strand
        
        
class Component(Base):
    __tablename__ = 'component'

    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    type = Column(String(10), nullable=False)

    __mapper_args__ = {'polymorphic_identity': 'component', 
                       'polymorphic_on': type
                       }

    def __init__(self, name):
        self.name = name
        
    def __repr__(self):
        return "Component (#%d):  %s" % \
            (self.id, self.name)


class DNA(Component):
    __tablename__ = 'dna'
    #genome_region = relationship(GenomeRegion)
    
    __mapper_args__ = { 'polymorphic_identity': 'DNA' }
    
    component_id = Column(Integer, ForeignKey('component.id'), primary_key=True)
    dna_type = Column(String(20))
    genome_region_id = Column(Integer, ForeignKey('genome_region.id'))
    
    
    def __init__(self, name, dna_type, leftpos, rightpos, strand):
        genome_region = get_or_create(Session(), GenomeRegion, leftpos=leftpos, rightpos=rightpos, strand=strand)
        super(DNA, self).__init__(name)
        self.dna_type = dna_type
    
    def __repr__(self):
        return "DNA (#%d, %s)" % \
            (self.component_id, self.name)   
    
    
class Gene(DNA):
    __tablename__ = 'gene'

    __mapper_args__ = { 'polymorphic_identity': 'gene' }
    
    component_id = Column(Integer, ForeignKey('dna.component_id'), primary_key=True)
    locus_id = Column(String(10))
    name = Column(String(10))
    long_name = Column(String(100))
    
    def __init__(self, name, dna_type, locus_id):
        super(Gene, self).__init__(name, dna_type)
        self.locus_id = locus_id


class DnaBindingSite(DNA):
    __tablename__ = 'dna_binding_site'
    __table_args__ = (PrimaryKeyConstraint('component_id','bound_component_id'),)


    __mapper_args__ = { 'polymorphic_identity': 'binding_site' }
    
    component_id = Column(Integer, ForeignKey('dna.component_id'))
    bound_component_id = Column(Integer, ForeignKey('component.id'))
    
    def __init__(self, name, dna_type, bound_component_id):
        super(DnaBindingSite, self).__init__(name, dna_type)
        self.bound_component_id = bound_component_id

    
class RNA(Component):
    __tablename__ = 'rna'
     
    __mapper_args__ = { 'polymorphic_identity': 'rna' }
     
    component_id = Column(Integer, ForeignKey('component.id'), primary_key=True)
    rna_type = Column(String(20))
    genome_region_id = Column(Integer, ForeignKey('genome_region.id'))
    
    
    def __init__(self, name, rna_type, leftpos, rightpos, strand):
        genome_region = get_or_create(Session(), GenomeRegion, leftpos=leftpos, rightpos=rightpos, strand=strand)
        super(RNA, self).__init__(name)
        self.rna_type = rna_type
     
    def __repr__(self):
        return "RNA (#%d, %s)" % \
            (self.id, self.name)      


class TU(RNA):
    __tablename__ = 'tu'

    __mapper_args__ = { 'polymorphic_identity': 'tu' }
    
    component_id = Column(Integer, ForeignKey('rna.component_id'), primary_key=True)
    tss = Column(Integer)
    genes = Column(String(10))
    name = Column(String(10))
    long_name = Column(String(100))
    
    def __init__(self, name):
        super(TU, self).__init__(name)
        self.name = name
   
     
class Protein(Component):
    __tablename__ = 'protein'
     
    __mapper_args__ = { 'polymorphic_identity': 'protein' }
     
    component_id = Column(Integer, ForeignKey('component.id'), primary_key=True)
    protein_type = Column(String(20))
    name = Column(String(10))
    long_name = Column(String(100))
    gene = Column(String(10))
    
    def __init__(self, name, protein_type, long_name):
        super(Protein, self).__init__(name)
        self.protein_type = dna_type
    
    def __repr__(self):
        return "Protein (#%d, %s)" % \
            (self.id, self.name)             
     
     
class Metabolite(Component):
    __tablename__ = 'metabolite'
     
    __mapper_args__ = { 'polymorphic_identity': 'metabolite' }
     
    component_id = Column(Integer, ForeignKey('component.id'), primary_key=True)
    metabolite_type = Column(String(20))
    name = Column(String(10))
    long_name = Column(String(100))
    formula = Column(String(100))
    
    def __init__(self, name, metabolite_type):
        super(Metabolite, self).__init__(name)
        self.metabolite_type = metabolite_type
     
     
    def __repr__(self):
        return "Metabolite (#%d, %s)" % \
            (self.id, self.name)     
    

Base.metadata.create_all()




    
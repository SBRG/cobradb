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
    type = Column(String(10))
    
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

    id = Column(Integer, ForeignKey('component.id'), primary_key=True)
    dna_type = Column(String(20))
    leftpos = Column(Integer)
    rightpos = Column(Integer)
    strand = Column(String(1))
    """The way to have this as genome_region would require altering the get_or_create()
    function to take in multiple class types and automatically join them when doing
    the filtering. This seems like overkill for now and its also not clear if genome_region
    should exist by itself"""
    #genome_region_id = Column(Integer, ForeignKey('genome_region.id'))
    #genome_region = relationship("GenomeRegion")
    
    __mapper_args__ = { 'polymorphic_identity': 'dna',
                        'polymorphic_on': dna_type
                      }
    
    
    def __init__(self, dna_type='dna', name=None, leftpos=None, rightpos=None, strand=None):
        super(DNA, self).__init__(name)
        self.dna_type = dna_type
        self.leftpos = leftpos
        self.rightpos = rightpos
        self.strand = strand
    
    def __repr__(self):
        return "DNA (#%d, %s) %d-%d %s"% \
            (self.id, self.name, self.leftpos, self.rightpos, self.strand)   
    
    
class Gene(DNA):
    __tablename__ = 'gene'

    __mapper_args__ = { 'polymorphic_identity': 'gene' }
    
    id = Column(Integer, ForeignKey('dna.id'), primary_key=True)
    locus_id = Column(String(10))
    name = Column(String(10))
    info = Column(String(200))
    long_name = Column(String(100))
    
    def __init__(self, name, leftpos, rightpos, strand, locus_id, info=None, long_name=None):
        super(Gene, self).__init__('gene', name, leftpos, rightpos, strand)
        self.locus_id = locus_id
        self.info = info
        self.long_name = long_name


class DnaBindingSite(DNA):
    __tablename__ = 'dna_binding_site'

    __mapper_args__ = { 'polymorphic_identity': 'binding_site' }
    
    id = Column(Integer, ForeignKey('dna.id'), primary_key=True)
    #bound_components = relationship("Component", secondary=dna_binding_bound_component_association,\
    #                                backref="dna_binding_site")
    
    
    def __init__(self, name, leftpos, rightpos, strand):
        super(DnaBindingSite, self).__init__('binding_site', name, leftpos, rightpos, strand)
        

   
class ComplexComposition(Base):
    __tablename__ = 'complex_composition'
    
    complex_id = Column(Integer, ForeignKey('complex.id'), primary_key=True)
    component_id = Column(Integer, ForeignKey('component.id'), primary_key=True)
    stoichiometry = Column(Integer)

    def __init__(self, complex_id, component_id, stoichiometry):
        self.complex_id = complex_id
        self.component_id = component_id
        self.stoichiometry = stoichiometry
        

class Complex(Component):
    __tablename__ = 'complex'
    
    __mapper_args__ = {'polymorphic_identity': 'complex'}
    
    id = Column(Integer, ForeignKey('component.id'), primary_key=True)

    children = relationship("Component", secondary="complex_composition",\
                            primaryjoin = id == ComplexComposition.complex_id,\
                            backref="parent")
    
    @hybrid_property
    def all_children(self):
        session = Session()
        included_components = session.query(
                                    ComplexComposition.complex_id,
                                    ComplexComposition.component_id).\
                            filter(ComplexComposition.complex_id == self.id).\
                            cte(name="included_components", recursive=True)

        
        incl_alias = aliased(included_components, name="incl_cplx")
        complex_alias = aliased(ComplexComposition, name="cplx")
        included_components = included_components.union_all(
                                                session.query(
                                                    complex_alias.complex_id,
                                                    complex_alias.component_id).\
                                                    filter(complex_alias.complex_id==incl_alias.c.component_id)
                                                )
        
        return session.query(Component).join(included_components, Component.id == included_components.c.component_id).all()
        
    #all_children = get_all_children(self)
        
    def __init__(self, name):
        super(Complex, self).__init__(name)

    
   
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
    name = Column(String(10))
    long_name = Column(String(100))
    #gene = Column(String(10))
    
    def __init__(self, name, long_name=None):
        super(Protein, self).__init__(name)
        self.long_name = long_name
    
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




    
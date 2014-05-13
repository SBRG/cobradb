from om.orm.base import *

from sqlalchemy.orm import relationship, backref
from sqlalchemy import Table, MetaData, create_engine, Column, Integer, \
    String, Float, ForeignKey, select
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.schema import UniqueConstraint,PrimaryKeyConstraint

import simplejson as json


class Gene(GenomeRegion):
    __tablename__ = 'gene'

    id = Column(Integer, ForeignKey('genome_region.id'), primary_key=True)
    locus_id = Column(String(10))
    info = Column(String(200))
    long_name = Column(String(100))
    
    __mapper_args__ = { 'polymorphic_identity': 'gene' }
    
    def __repr__(self):
        return "Gene: (%s, %s) %d-%d (%s)"% \
            (self.locus_id, self.name, self.leftpos, self.rightpos,\
                                 self.strand)   
    
    
    def __init__(self, name, leftpos, rightpos, strand, locus_id, info=None, long_name=None):
        super(Gene, self).__init__(leftpos, rightpos, strand, name)
        self.locus_id = locus_id
        self.info = info
        self.long_name = long_name


class Motif(GenomeRegion):
    __tablename__ = 'motif'
    
    id = Column(Integer, ForeignKey('genome_region.id'), primary_key=True)
    pval = Column(Float)
    
    bound_component_id = Column(Integer, ForeignKey('component.id'))
    bound_component = relationship("Component")
    
    def __repr__(self):
        return "Motif (%s) %d-%d %s %5.2f"% \
            (self.bound_component.name, self.leftpos, self.rightpos,\
                                 self.strand, self.pval)   
    
    def __init__(self, leftpos, rightpos, strand, pval, info=None):
        super(Motif, self).__init__(leftpos, rightpos, strand)
        self.pval = pval
        self.info = info


class Component(Base):
    __tablename__ = 'component'

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(100))
    type = Column(String(20))
   
    __table_args__ = (UniqueConstraint('name'),{})
    
    __mapper_args__ = {'polymorphic_identity': 'component', 
                       'polymorphic_on': type
                      }

    def __init__(self, name):
        self.name = name
        
    def __repr__(self):
        return "Component (#%d):  %s" % \
            (self.id, self.name)

   
class ComplexComposition(Base):
    __tablename__ = 'complex_composition'
    
    complex_id = Column(Integer, ForeignKey('complex.id'), primary_key=True)
    component_id = Column(Integer, ForeignKey('component.id'), primary_key=True)
    stoichiometry = Column(Integer)

    __table_args__ = (UniqueConstraint('complex_id','component_id'),{})

    def __init__(self, complex_id, component_id, stoichiometry):
        self.complex_id = complex_id
        self.component_id = component_id
        self.stoichiometry = stoichiometry
        

class Complex(Component):
    __tablename__ = 'complex'
    
    __mapper_args__ = {'polymorphic_identity': 'complex'}
    
    id = Column(Integer, ForeignKey('component.id'), primary_key=True)
    
    long_name = Column(String(200)) 
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
        
        
    def __repr__(self):
        return "Complex (#%d):  %s" % \
            (self.id, self.long_name)
                
        
    def __init__(self, name, long_name=None):
        super(Complex, self).__init__(name)
        self.long_name = long_name


class DNA(Component):
    __tablename__ = 'dna'

    id = Column(Integer, ForeignKey('component.id'), primary_key=True)
    type = Column(String(20))

    genome_region_id = Column(Integer, ForeignKey('genome_region.id'))
    genome_region = relationship('GenomeRegion', backref=backref('dna', lazy='dynamic'))
    
    __mapper_args__ = { 'polymorphic_identity': 'dna',
                        'polymorphic_on': type
                      }
    
    
    def __init__(self, name=None, leftpos=None, rightpos=None, strand=None):
        super(DNA, self).__init__(name)
        session = Session()
        self.genome_region_id = session.get_or_create(GenomeRegion, leftpos=leftpos,\
                                              rightpos=rightpos, strand=strand).id
        session.close()
        
    def __repr__(self):
        return "DNA (#%d, %s) %d-%d %s"% \
            (self.id, self.name, self.genome_region.leftpos, self.genome_region.rightpos,\
                                 self.genome_region.strand)   
    

class DnaBindingSite(DNA):
    __tablename__ = 'dna_binding_site'

    __mapper_args__ = { 'polymorphic_identity': 'binding_site' }
    
    id = Column(Integer, ForeignKey('dna.id'), primary_key=True)
    centerpos = Column(Integer)
    width = Column(Integer)
    
    
    def __init__(self, name, leftpos, rightpos, strand, centerpos, width):
        super(DnaBindingSite, self).__init__(name, leftpos, rightpos, strand)
        self.centerpos = centerpos
        self.width = width
        
    
class RNA(Component):
    __tablename__ = 'rna'
     
    __mapper_args__ = { 'polymorphic_identity': 'rna' }
     
    id = Column(Integer, ForeignKey('component.id'), primary_key=True)
    rna_type = Column(String(20))
    genome_region_id = Column(Integer, ForeignKey('genome_region.id'))
    
    
    def __init__(self, name=None, leftpos=None, rightpos=None, strand=None):
        super(RNA, self).__init__(name)
        session = Session()
        self.genome_region_id = session.get_or_create(GenomeRegion, leftpos=leftpos,\
                                              rightpos=rightpos, strand=strand).id
        session.close()
        
    def __repr__(self):
        return "RNA (#%d, %s)" % \
            (self.id, self.name)      


class TUGenes(Base):
    __tablename__ = 'tu_genes'
    
    tu_id = Column(Integer, ForeignKey('tu.id'), primary_key=True)
    gene_id = Column(Integer, ForeignKey('gene.id'), primary_key=True)
    
    __table_args__ = (UniqueConstraint('tu_id','gene_id'),{})
    
    def __init__(self, tu_id, gene_id):
        self.tu_id = tu_id
        self.gene_id = gene_id


class TU(RNA):
    __tablename__ = 'tu'

    __mapper_args__ = { 'polymorphic_identity': 'tu' }
    
    id = Column(Integer, ForeignKey('rna.id'), primary_key=True)
    genome_region = relationship("GenomeRegion")
    genes = relationship("Gene", secondary="tu_genes",\
                                 primaryjoin = id == TUGenes.tu_id,\
                                 backref="tu")
    
    name = Column(String(100))
    
    """
    @hybrid_property
    def tss(self):
        if self.genome_region.strand == '+':
            return self.genome_region.leftpos
        else:
            return self.genome_region.rightpos
    """
    
    def __init__(self, name, leftpos, rightpos, strand):
        super(TU, self).__init__(name, leftpos, rightpos, strand)
        self.name = name
   
     
class Protein(Component):
    __tablename__ = 'protein'
     
    __mapper_args__ = { 'polymorphic_identity': 'protein' }
     
    id = Column(Integer, ForeignKey('component.id'), primary_key=True)
    long_name = Column(String(200))
    gene_id = Column(Integer, ForeignKey('gene.id'))
    gene = relationship('Gene', backref='protein')
    
    def __init__(self, name, gene_id=None, long_name=None):
        super(Protein, self).__init__(name)
        self.gene_id = gene_id
        self.long_name = long_name
    
    def __repr__(self):
        return "Protein (#%d, %s)" % \
            (self.id, self.long_name)             
     
     
class SmallMolecule(Component):
    __tablename__ = 'small_molecule'
     
    __mapper_args__ = { 'polymorphic_identity': 'small_molecule' }
     
    id = Column(Integer, ForeignKey('component.id'), primary_key=True)

    long_name = Column(String(100))
    formula = Column(String(100))
    smiles = Column(String(200))
    def __init__(self, name, long_name, formula="", smiles=""):
        super(SmallMolecule, self).__init__(name)
        self.long_name = long_name
        self.formula = formula
        self.smiles = smiles
     
     
    def __repr__(self):
        return "Small Molecule (#%d, %s)" % \
            (self.id, self.long_name)     
    






    
from ome.orm.base import *

from sqlalchemy.orm import relationship
from sqlalchemy import Table, MetaData, create_engine, Column, Integer, \
    String, Float, ForeignKey, select
from sqlalchemy.ext.hybrid import hybrid_property

import simplejson as json


class Component(Base):
    __table__ = make_table("component")

    __mapper_args__ = {
        'polymorphic_identity': 'component',
        'polymorphic_on': __table__.c.component_type_id}

    def __repr__(self):
        return "Component (#%d):  %s" % \
            (self.id, self.name)


class DNA(Component):
    __table__ = make_table("dna")
    
    __mapper_args__ = { 
        'polymorphic_identity': 1,
        'polymorphic_on': __table__.c.dna_type_id }
    
    def __repr__(self):
        return "DNA (#%d, %s)" % \
            (self.id, self.name)   


class DnaBindingSite(DNA):
    __table__ = make_table("dna_binding_site")
    
    __mapper_args__ = { 'polymorphic_identity': 1 }
    

class RNA(Component):
    __table__ = make_table("rna")
    
    __mapper_args__ = { 'polymorphic_identity': 2 }
    
    def __repr__(self):
        return "RNA (#%d, %s)" % \
            (self.id, Component.name)      
    
    
class Protein(Component):
    __table__ = make_table("protein")
    
    __mapper_args__ = { 'polymorphic_identity': 3 }
    
    def __repr__(self):
        return "Protein (#%d, %s)" % \
            (self.id, Component.name)   
    
    
    
class Metabolite(Component):
    __table__ = make_table("metabolite")
    
    __mapper_args__ = { 'polymorphic_identity': 4 }
    
    def __repr__(self):
        return "Metabolite (#%d, %s)" % \
            (self.id, Component.name)     
    

class GenomeRegion(Base):
    __table__ = make_table("genome_region")
    
    


class Gene(Base):
    __table__ = make_table("gene")



    
from om.orm.base import *

from sqlalchemy.orm import relationship, backref
from sqlalchemy import Table, MetaData, create_engine, Column, Integer, \
    String, Float, ForeignKey, select
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.schema import UniqueConstraint,PrimaryKeyConstraint

import simplejson as json


class Compartment(Base):
    __tablename__ = 'compartment'

    id = Column(Integer, primary_key=True)
    name = Column(String(20))
    abbreviation = Column(String(5))
    
    __table_args__ = (UniqueConstraint('name'),{})

    def __repr__(self):
        return "Compartment (#%d): %s, %s"% \
            (self.id, self.name, self.abbreviation)
    
    
    def __init__(self, name, abbreviation):
        self.name = name
        self.abbreviation = abbreviation

        
class ModelComponent(Base):
    __tablename__ = 'model_component'

    id = Column(Integer, primary_key=True)
    component_id = Column(Integer, ForeignKey('component.id'))
    component = relationship('Component')
    compartment_id = Column(Integer, ForeignKey('compartment.id'))
    compartment = relationship('Compartment')
    model_id = Column(Integer, ForeignKey('model.id'))
    model = relationship('Model')
    
    __table_args__ = (UniqueConstraint('component_id','compartment_id','model_id'),{})
    
    def __repr__(self):
        return "Model Component (#%d): %s, %s, %s"% \
            (self.id, self.model, self.component, self.compartment)   
    
    
    def __init__(self, component_id, compartment_id, model_id):
        self.component_id = component_id
        self.compartment_id = compartment_id
        self.model_id = model_id


class ReactionMatrix(Base):
    __tablename__ = 'reaction_matrix'

    model_reaction_id = Column(Integer, ForeignKey('model_reaction.id'), primary_key=True)
    model_component_id = Column(Integer, ForeignKey('model_component.id'), primary_key=True)
    stoichiometry = Column(Float)
        
    __table_args__ = (UniqueConstraint('model_reaction_id','model_component_id'),{})

    def __repr__(self):
        return "Reaction Matrix (#%d): %d, %d"% \
            (self.model_reaction_id, self.model_component_id, self.stoichiometry)   
    
    
    def __init__(self, model_reaction_id, model_component_id, stoichiometry):
        self.model_reaction_id = model_reaction.id
        self.model_component_id = model_component.id
        self.stoichiometry = stoichiometry


class ModelReaction(Base):
    __tablename__ = 'model_reaction'

    id = Column(Integer, primary_key=True, autoincrement=True)
    reaction_id = Column(Integer, ForeignKey('reaction.id'))
    reaction = relationship('Reaction')
    model_id = Column(Integer, ForeignKey('model.id'))
    model = relationship('Model')
    
    __table_args__ = (UniqueConstraint('reaction_id','model_id'),{})
    
    def __repr__(self):
        return "Model Reaction (#%d): %s, %s"% \
            (self.id, self.reaction, self.model)   
    
    
    def __init__(self, reaction_id, model_id):
        self.reaction_id = reaction_id
        self.model_id = model_id
        
        
class GPR(Base):
    __tablename__ = 'gpr'
    
    model_reaction_id = Column(Integer, ForeignKey('model_reaction.id'), primary_key=True)
    gene_id = Column(Integer, ForeignKey('genome_region.id'), primary_key=True)

    __table_args__ = (UniqueConstraint('model_reaction_id','gene_id'),{})

    def __repr__(self):
        return "GPR: %d, %d"% \
            (self.model_reaction_id, self.gene_id)   
    
    def __init__(self, model_reaction_id, gene_id):
        self.model_reaction_id
        self.gene_id
        

class Model(Base):
    __tablename__ = 'model'

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(100))
    bigg_id = Column(String(100))
    version = Column(Float)
    
    __table_args__ = (UniqueConstraint('name'),{})

    def __repr__(self):
        return "Model (#%d): %s, %s, %5.2f"% \
            (self.id, self.name, self.bigg_id, self.version)   
    
    def __init__(self, name, bigg_id, version):
        self.name = name
        self.bigg_id = bigg_id
        self.version = version


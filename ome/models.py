from sqlalchemy import create_engine, ForeignKey, Column, Integer, String, Numeric, Table, MetaData, DateTime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.schema import UniqueConstraint
from sqlalchemy.schema import Sequence
from ome.base import *
from ome.components import *

#session.execute("CREATE EXTENSION pg_trgm;")

"""
CREATE EXTENSION pg_trgm;
CREATE INDEX gene_locus_trigram_idx ON gene USING gin (to_tsvector('english',locus_id));
CREATE INDEX reaction_name_trigram_idx ON reaction USING gin (to_tsvector('english',name));
CREATE INDEX genome_name_trigram_idx ON genome_region USING gin (to_tsvector('english',name));
CREATE INDEX model_biggid_trigram_idx ON model USING gin (to_tsvector('english',bigg_id));
CREATE INDEX genome_organism_trigram_idx ON genome USING gin (to_tsvector('english',organism));
CREATE INDEX metabolite_name_trigram_idx ON metabolite USING gin (to_tsvector('english',name));
"""
class Model(Base):
    __tablename__='model'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    bigg_id = Column(String)
    first_created = Column(DateTime)
    genome_id = Column(Integer, ForeignKey('genome.id'))
    genome = relationship('Genome', backref='model')
    notes = Column(String)

    __table_args__ = (UniqueConstraint('bigg_id'),{})

    def __repr__(self):
        return "Model (#%d) %s %s" % (self.id, self.bigg_id, self.first_created)



class ModelGene(Base):
    __tablename__='model_gene'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    model_id = Column(Integer, ForeignKey('model.id'), nullable=False)
    gene_id = Column(Integer, ForeignKey('gene.id'), nullable=False)



class ModelReaction(Base):
    __tablename__='model_reaction'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    reaction_id = Column(Integer, ForeignKey('reaction.id'), nullable=False)
    model_id = Column(Integer, ForeignKey('model.id'), nullable=False)
    name = Column(String)
    upperbound = Column(Numeric)
    lowerbound = Column(Numeric)
    gpr = Column(String)
    __table_args__ = (UniqueConstraint('reaction_id', 'model_id'),{})
    #UniqueConstraint('reaction_id', 'model_id')



class GPRMatrix(Base):
    __tablename__='gpr_matrix'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    model_gene_id = Column(Integer, ForeignKey('model_gene.id'), nullable=False)
    model_reaction_id = Column(Integer, ForeignKey('model_reaction.id'), nullable=False)



class CompartmentalizedComponent(Base):
    __tablename__='compartmentalized_component'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    component_id = Column(Integer, ForeignKey('component.id'), nullable=False)
    compartment_id = Column(Integer, ForeignKey('compartment.id'), nullable=False)
    #UniqueConstraint('compartment_id', 'component_id')
    __table_args__ = (UniqueConstraint('compartment_id', 'component_id'),{})


class ModelCompartmentalizedComponent(Base):
    __tablename__='model_compartmentalized_component'
    id = Column(Integer, Sequence('wids'), primary_key=True)
    model_id = Column(Integer, ForeignKey('model.id'), nullable=False)
    compartmentalized_component_id = Column(Integer, ForeignKey('compartmentalized_component.id'), nullable=False)
    compartment_id = Column(Integer, ForeignKey('compartment.id'), nullable=False)



class Compartment(Base):
    __tablename__='compartment'
    id = Column(Integer, Sequence('wids'), primary_key=True)
    name = Column(String, unique = True)



class ReactionMatrix(Base):
    __tablename__='reaction_matrix'
    id = Column(Integer, Sequence('wids'), primary_key=True)
    reaction_id = Column(Integer, ForeignKey('reaction.id'), nullable=False)
    compartmentalized_component_id = Column(Integer, ForeignKey('compartmentalized_component.id'), nullable=False)
    stoichiometry = Column(Numeric)
    #UniqueConstraint('reaction_id', 'compartmentalized_component')
    __table_args__ = (UniqueConstraint('reaction_id', 'compartmentalized_component_id'),{})


class EscherMap(Base):
    __tablename__='escher_map'
    id = Column(Integer, primary_key=True)
    biggid = Column(String)
    category = Column(String)
    model_name = Column(String)



class Comments(Base):
    __tablename__ = 'comments'
    id= Column(Integer, primary_key=True)
    kegg_id=Column(String)
    cas_number=Column(String)
    name = Column(String)
    formula = Column(String)
    text = Column(String)



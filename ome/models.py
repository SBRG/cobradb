from sqlalchemy import (create_engine, ForeignKey, Column, Integer, String,
                        Numeric, Table, MetaData, DateTime, LargeBinary)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.schema import UniqueConstraint
from sqlalchemy.schema import Sequence
from ome.base import *
from ome.components import *


class Model(Base):
    __tablename__='model'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    bigg_id = Column(String, nullable=False)
    genome_id = Column(Integer, ForeignKey('genome.id', onupdate="CASCADE", ondelete="CASCADE"))
    genome = relationship('Genome', backref='model')
    description = Column(String, nullable=False)
    published_filename = Column(String, nullable=True)

    __table_args__ = (
        UniqueConstraint('bigg_id', 'genome_id'),
    )

    def __repr__(self):
        return '<ome Model(id={self.id}, bigg_id={self.bigg_id})>'.format(self=self)


class ModelGene(Base):
    __tablename__='model_gene'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    model_id = Column(Integer,
                      ForeignKey('model.id', onupdate="CASCADE", ondelete="CASCADE"),
                      nullable=False)
    gene_id = Column(Integer,
                     ForeignKey('gene.id', onupdate="CASCADE", ondelete="CASCADE"),
                     nullable=False)

    __table_args__ = (
        UniqueConstraint('model_id', 'gene_id'),
    )


class ModelReaction(Base):
    __tablename__='model_reaction'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    reaction_id = Column(Integer,
                         ForeignKey('reaction.id', onupdate="CASCADE", ondelete="CASCADE"),
                         nullable=False)
    model_id = Column(Integer,
                      ForeignKey('model.id', onupdate="CASCADE", ondelete="CASCADE"),
                      nullable=False)
    copy_number = Column(Integer, nullable=False)

    objective_coefficient = Column(Numeric, nullable=False)
    lower_bound = Column(Numeric, nullable=False)
    upper_bound = Column(Numeric, nullable=False)
    gene_reaction_rule = Column(String, nullable=False)
    original_gene_reaction_rule = Column(String, nullable=True)

    __table_args__ = (
        UniqueConstraint('reaction_id', 'model_id', 'copy_number'),
    )

    def __repr__(self):
        return ('<ome ModelReaction(id={self.id}, reaction_id={self.reaction_id}, model_id={self.model_id}, copy_number={self.copy_number})>'
                .format(self=self))


class GeneReactionMatrix(Base):
    __tablename__='gene_reaction_matrix'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    model_gene_id = Column(Integer,
                           ForeignKey('model_gene.id', onupdate="CASCADE", ondelete="CASCADE"),
                           nullable=False)
    model_reaction_id = Column(Integer,
                               ForeignKey('model_reaction.id', onupdate="CASCADE", ondelete="CASCADE"),
                               nullable=False)

    __table_args__ = (UniqueConstraint('model_gene_id', 'model_reaction_id'), {})

    def __repr__(self):
        return ('<ome GeneReactionMatrix(id=%d, model_gene_id=%d, model_reaction_id=%d)>' %
                (self.id, self.model_gene_id, self.model_reaction_id))


class CompartmentalizedComponent(Base):
    __tablename__='compartmentalized_component'
    id = Column(Integer, Sequence('wids'), primary_key=True)
    component_id = Column(Integer,
                          ForeignKey('component.id', onupdate="CASCADE", ondelete="CASCADE"),
                          nullable=False)
    compartment_id = Column(Integer,
                            ForeignKey('compartment.id', onupdate="CASCADE", ondelete="CASCADE"),
                            nullable=False)

    __table_args__ = (
        UniqueConstraint('compartment_id', 'component_id'),
    )


class ModelCompartmentalizedComponent(Base):
    __tablename__='model_compartmentalized_component'
    id = Column(Integer, Sequence('wids'), primary_key=True)
    model_id = Column(Integer,
                      ForeignKey('model.id', onupdate="CASCADE", ondelete="CASCADE"),
                      nullable=False)
    compartmentalized_component_id = Column(Integer,
                                            ForeignKey('compartmentalized_component.id'),
                                            nullable=False)

    __table_args__ = (
        UniqueConstraint('compartmentalized_component_id', 'model_id'),
    )


class Compartment(Base):
    __tablename__ = 'compartment'
    id = Column(Integer, Sequence('wids'), primary_key=True)
    bigg_id  = Column(String, unique = True)
    name = Column(String)


class ReactionMatrix(Base):
    __tablename__ = 'reaction_matrix'
    id = Column(Integer, Sequence('wids'), primary_key=True)
    reaction_id = Column(Integer, ForeignKey('reaction.id'), nullable=False)
    compartmentalized_component_id = Column(Integer,
                                            ForeignKey('compartmentalized_component.id',
                                                       onupdate="CASCADE", ondelete="CASCADE"),
                                            nullable=False)
    stoichiometry = Column(Numeric)

    __table_args__ = (
        UniqueConstraint('reaction_id', 'compartmentalized_component_id'),
    )


class EscherMap(Base):
    __tablename__ = 'escher_map'
    id = Column(Integer, Sequence('wids'), primary_key=True)
    map_name = Column(String, nullable=False)
    map_data = Column(LargeBinary, nullable=False)
    model_id = Column(Integer, ForeignKey(Model.id), nullable=False)
    priority = Column(Integer, nullable=False)

    __table_args__ = (
        UniqueConstraint('map_name'),
    )


class EscherMapMatrix(Base):
    __tablename__ = 'escher_map_matrix'
    id = Column(Integer, Sequence('wids'), primary_key=True)
    ome_id = Column(Integer, nullable=False)
    escher_map_id = Column(Integer, ForeignKey(EscherMap.id), nullable=False)
    # the reaction id or node id
    escher_map_element_id = Column(String(50))
    type = Column(String, nullable=False)

    __table_args__ = (
        UniqueConstraint('ome_id', 'escher_map_id'),
    )


class Comments(Base):
    __tablename__ = 'comments'
    id= Column(Integer, primary_key=True)
    email = Column(String)
    text = Column(String)
    date_created = Column(DateTime)
    type = Column(String)
    url = Column(String)

class ModelCount(Base):
    __tablename__='model_count'
    id = Column(Integer, primary_key=True)
    model_id = Column(Integer,
                      ForeignKey('model.id', onupdate="CASCADE", ondelete="CASCADE"),
                      nullable=False)
    reaction_count = Column(Integer)
    gene_count = Column(Integer)
    metabolite_count = Column(Integer)

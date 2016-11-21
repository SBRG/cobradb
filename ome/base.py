"""Module to implement ORM to the ome database"""

from ome import settings

from sqlalchemy.orm import sessionmaker, relationship, aliased
from sqlalchemy.orm.session import Session as _SA_Session
from sqlalchemy import (Table, MetaData, create_engine, Column, Integer, String,
                        Float, Numeric, ForeignKey, Boolean, Enum, DateTime)
from sqlalchemy.schema import UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
from types import MethodType
from os import system
from contextlib import contextmanager
import logging


# connect to postgres
engine = create_engine(settings.db_connection_string)
Base = declarative_base(bind=engine)
metadata = MetaData(bind=engine)
Session = sessionmaker(bind=engine, class_=_SA_Session)


# make the enums
_enum_l = [
    Enum('component', 'reaction', 'gene', 'compartmentalized_component',
         name='synonym_type'),
    Enum('pmid', 'doi',
         name='reference_type'),
    Enum('model_reaction', 'model_compartmentalized_component', 'model_gene',
         name='old_id_synonym_type'),
    Enum('is_version', name='is_version'),
    Enum('component', 'reaction', name='deprecated_id_types'),
    Enum('model_compartmentalized_component', 'model_reaction',
         name='escher_map_matrix_type')
]
custom_enums = { x.name: x for x in _enum_l }


# exceptions
class NotFoundError(Exception):
    pass


class DatabaseVersion(Base):
    __tablename__ = 'database_version'

    is_version = Column(custom_enums['is_version'], primary_key=True)
    date_time = Column(DateTime, nullable=False)

    __table_args__ = (
        UniqueConstraint('is_version'),
    )

    def __init__(self, date_time):
        self.is_version = 'is_version'
        self.date_time = date_time


class Genome(Base):
    __tablename__ = 'genome'

    id = Column(Integer, primary_key=True)
    accession_type = Column(String(200), nullable=False)
    accession_value = Column(String(200), nullable=False)
    organism = Column(String(200), nullable=True)
    taxon_id = Column(String(200), nullable=True)
    ncbi_assembly_id = Column(String(200), nullable=True)

    __table_args__ = (
        UniqueConstraint('accession_type', 'accession_value'),
    )

    def __repr__(self):
        return ('<ome Genome(id={self.id}, accession_type={self.accession_type}, '
                'accession_value={self.accession_value})>'.format(self=self))


class Chromosome(Base):
    __tablename__ = 'chromosome'

    id = Column(Integer, primary_key=True)
    ncbi_accession = Column(String(200))
    genome_id = Column(Integer, ForeignKey('genome.id'))

    __table_args__ = (
        UniqueConstraint('ncbi_accession', 'genome_id'),
    )

    def __repr__(self):
        return ('<ome Chromosome(id={self.id}, ncbi_accession={self.ncbi_accession}, genome_id={self.genome_id})>'
                .format(self=self))


class GenomeRegion(Base):
    __tablename__ = 'genome_region'
    id = Column(Integer, primary_key=True)
    chromosome_id = Column(Integer, ForeignKey('chromosome.id'))
    bigg_id = Column(String, nullable=False)
    leftpos = Column(Integer, nullable=True)
    rightpos = Column(Integer, nullable=True)
    strand = Column(String(1), nullable=True)
    type = Column(String(20))

    __table_args__ = (
        UniqueConstraint('bigg_id', 'chromosome_id'),
    )

    __mapper_args__ = {
        'polymorphic_identity': 'genome_region',
        'polymorphic_on': type
    }

    def __repr__(self):
        return ('<ome GenomeRegion(id={self.id}, leftpos={self.leftpos}, rightpos={self.rightpos})>'
                .format(self=self))


class Component(Base):
    __tablename__ = 'component'

    id = Column(Integer, primary_key=True)
    bigg_id = Column(String)
    name = Column(String, nullable=True)
    type = Column(String(20))

    __table_args__ = (UniqueConstraint('bigg_id'), {})

    __mapper_args__ = {
        'polymorphic_identity': 'component',
        'polymorphic_on': type
    }

    def __repr__(self):
        return "Component (#%d):  %s" % \
            (self.id, self.name)


class Reaction(Base):
    __tablename__ = 'reaction'

    id = Column(Integer, primary_key=True)
    type = Column(String(20))
    bigg_id = Column(String, nullable=False)
    name = Column(String, nullable=True)
    reaction_hash = Column(String, nullable=False)
    pseudoreaction = Column(Boolean, default=False)

    __table_args__ = (
        UniqueConstraint('bigg_id'),
    )

    __mapper_args__ = {
        'polymorphic_identity': 'reaction',
        'polymorphic_on': type
    }

    def __repr__(self):
        return ('<ome Reaction(id=%d, bigg_id=%s%s)>' %
                (self.id, self.bigg_id, ', pseudoreaction' if self.pseudoreaction else ''))


class DataSource(Base):
    __tablename__ = 'data_source'

    id = Column(Integer, primary_key=True)
    bigg_id = Column(String, nullable=False)
    name = Column(String(100))
    url_prefix = Column(String)

    __table_args__ = (
        UniqueConstraint('bigg_id'),
    )

    def __repr__(self):
        return (
            '<ome DataSource(id={self.id}, bigg_id={self.bigg_id}, '
            'name={self.name}, url_prefix={self.url_prefix})>'
        ).format(self=self)


class Synonym(Base):
    __tablename__ = 'synonym'
    id = Column(Integer, primary_key=True)
    ome_id = Column(Integer)
    synonym = Column(String)
    type = Column(custom_enums['synonym_type'])
    data_source_id = Column(Integer, ForeignKey('data_source.id', ondelete='CASCADE'))

    __table_args__ = (
        UniqueConstraint('ome_id', 'synonym', 'type', 'data_source_id'),
    )

    def __repr__(self):
        return ('<ome Synonym(id=%d, synonym="%s", type="%s", ome_id=%d, data_source_id=%d)>' %
                (self.id, self.synonym, self.type, self.ome_id, self.data_source_id))


class Publication(Base):
    __tablename__ = "publication"
    id = Column(Integer, primary_key=True)
    reference_type = Column(custom_enums['reference_type'])
    reference_id = Column(String)

    __table_args__=(
        UniqueConstraint('reference_type', 'reference_id'),
    )


class PublicationModel(Base):
    __tablename__ = "publication_model"
    model_id = Column(Integer,
                      ForeignKey('model.id', ondelete='CASCADE'),
                      primary_key=True)
    publication_id = Column(Integer,
                            ForeignKey('publication.id', ondelete='CASCADE'),
                            primary_key=True)

    __table_args__ = (
        UniqueConstraint('model_id', 'publication_id'),
    )


class OldIDSynonym(Base):
    __tablename__ = "old_id_model_synonym"
    id = Column(Integer, primary_key=True)
    type = Column(custom_enums['old_id_synonym_type'])
    synonym_id = Column(Integer,
                        ForeignKey('synonym.id', ondelete='CASCADE'),
                        nullable=False)
    ome_id = Column(Integer, nullable=False)

    __table_args__ = (
        UniqueConstraint('synonym_id', 'ome_id'),
    )

    def __repr__(self):
        return ('<ome OldIDSynonym(id=%d, type="%s", ome_id=%d, synonym_id=%d)>' %
                (self.id, self.type, self.ome_id, self.synonym_id))


class GenomeRegionMap(Base):
        __tablename__ = 'genome_region_map'

        genome_region_id_1 = Column(Integer, ForeignKey('genome_region.id'), primary_key=True)
        genome_region_id_2 = Column(Integer, ForeignKey('genome_region.id'), primary_key=True)
        distance = Column(Integer)

        __table_args__ = (
            UniqueConstraint('genome_region_id_1','genome_region_id_2'),
        )

        def __repr__(self):
            return "GenomeRegionMap (%d <--> %d) distance:%d" % (self.genome_region_id_1, self.genome_region_id_2, self.distance)

class DeprecatedID(Base):
    __tablename__ = 'deprecated_id'

    id = Column(Integer, primary_key=True)
    type = Column(custom_enums['deprecated_id_types'])
    deprecated_id = Column(String)
    ome_id = Column(Integer)

    __table_args__ = (
        UniqueConstraint('type', 'deprecated_id'),
    )

    def __repr__(self):
        return ('<ome DeprecatedID(type="%s", deprecated_id="%s", component_id=%d)>' %
                (self.type, self.deprecated_id, self.component_id))

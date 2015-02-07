"""Module to implement ORM to the ome database"""

from ome import settings

from types import MethodType
from os import system
from sqlalchemy.orm import sessionmaker, relationship, aliased
from sqlalchemy.orm.session import Session as _SA_Session
from sqlalchemy import (Table, MetaData, create_engine, Column, Integer, String,
                        Float, Numeric, ForeignKey, Boolean)
from sqlalchemy.schema import UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
from contextlib import contextmanager
from sqlalchemy.schema import Sequence
from warnings import warn

try:
    import pymongo
    MONGO_INSTALLED = True
except ImportError:
    warn('pymongo not installed')
    MONGO_INSTALLED = False

engine = create_engine("postgresql://%s:%s@%s/%s" %
    (settings.postgres_user, settings.postgres_password, settings.postgres_host, settings.postgres_database))
Base = declarative_base(bind=engine)
metadata = MetaData(bind=engine)

if MONGO_INSTALLED:
    try:
        connection = pymongo.Connection()
        omics_database = connection.omics_database
    except Exception as e:
        warn("Failed to connect to mongo with error: " + e.message)
        omics_database = None
else:
    omics_database = None

class Genome(Base):
    __tablename__ = 'genome'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    bioproject_id = Column(String(200))
    organism = Column(String(200))

    __table_args__ = (UniqueConstraint('bioproject_id'),{})


    def __repr__(self):
        return "Genome (#%d) %s %s" % (self.id, self.bioproject_id, self.organism)

    def __init__(self, bioproject_id, organism):
        self.bioproject_id = bioproject_id
        self.organism = organism


class Chromosome(Base):
    __tablename__ = 'chromosome'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    genome_id = Column(Integer, ForeignKey('genome.id'))
    genome = relationship('Genome', backref='chromosomes')
    genbank_id = Column(String(100))
    ncbi_id = Column(String(100))

    __table_args__ = (UniqueConstraint('genome_id', 'genbank_id'),{})


    def __repr__(self):
        return "Chromosome %s -- %s" % (self.ncbi_id, self.genome)


    def __init__(self, genome_id, genbank_id, ncbi_id):
        self.genome_id = genome_id
        self.genbank_id = genbank_id
        self.ncbi_id = ncbi_id


class GenomeRegion(Base):
    __tablename__ = 'genome_region'
    id = Column(Integer, Sequence('wids'), primary_key=True)
    chromosome_id = Column(Integer, ForeignKey('chromosome.id'))
    bigg_id = Column(String(20))
    leftpos = Column(Integer)
    rightpos = Column(Integer)
    strand = Column(String(1))
    type = Column(String(20))
    __table_args__ = (UniqueConstraint('bigg_id','leftpos','rightpos','strand','chromosome_id'),{})

    __mapper_args__ = {'polymorphic_identity': 'genome_region',
                       'polymorphic_on': type
                      }

    def __repr__(self):
        return "GenomeRegion: %d-%d (%s)" % \
                (self.leftpos, self.rightpos, self.strand)


    def __init__(self, leftpos, rightpos, strand, chromosome_id, bigg_id):
        self.leftpos = leftpos
        self.rightpos = rightpos
        self.strand = strand
        self.chromosome_id = chromosome_id
        self.bigg_id = bigg_id

        
class Component(Base):
    __tablename__ = 'component'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    bigg_id = Column(String)
    name = Column(String)
    type = Column(String(20))

    __table_args__ = (UniqueConstraint('bigg_id'), {})

    __mapper_args__ = {
        'polymorphic_identity': 'component',
        'polymorphic_on': type
    }

    def __init__(self, bigg_id, name):
        self.bigg_id = bigg_id
        self.name = name

    def __repr__(self):
        return "Component (#%d):  %s" % \
            (self.id, self.name)

        
class Reaction(Base):
    __tablename__ = 'reaction'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    bigg_id = Column(String)
    name = Column(String)
    type = Column(String(20))
    notes = Column(String)
    reaction_hash = Column(String)
    objective_coefficient = Column(Numeric)
    lower_bound = Column(Numeric)
    upper_bound = Column(Numeric)
    __table_args__ = (UniqueConstraint('bigg_id'),{})

    __mapper_args__ = {
        'polymorphic_identity': 'reaction',
        'polymorphic_on': type
    }

    def __init__(self, name, notes, reaction_hash, bigg_id):
        self.name = name
        self.bigg_id = bigg_id
        self.notes = notes
        self.reaction_hash = reaction_hash

    def __repr__(self):
        return "Reaction (#%d):  %s" % \
            (self.id, self.name)


class DataSource(Base):
    __tablename__ = 'data_source'

    id = Column(Integer, Sequence('wids'), primary_key=True)
    name = Column(String(100))
    lab = Column(String(100))
    institution = Column(String(100))
    #data_sets = relationship("DataSet")

    __table_args__ = (UniqueConstraint('name'),{})

    def __repr__(self):
        return "Data Source %s (#%d)" % (self.name, self.id)

    def __repr__dict__(self):
        return {"name":self.name,"wid":self.id,"values":{"lab":self.lab,"institution":self.institution}}

    def __repr__json__(self):
        return json.dumps(self.__repr__dict__())

    def __init__(self, name, lab=None, institution=None):
        self.name = name
        self.lab = lab
        self.institution = institution


class Synonym(Base):
    __tablename__ = "synonym"
    id = Column(Integer, Sequence('wids'), primary_key=True)
    ome_id = Column(Integer)
    synonym = Column(String)
    type = Column(String)
    synonym_data_source_id = Column(Integer, ForeignKey('data_source.id', ondelete='CASCADE'))
    synonym_data_source = relationship("DataSource")


    __table_args__ = (UniqueConstraint('ome_id','synonym','type'),{})

    def __repr__(self):
        return "%s in (%s)" % (self.synonym, self.synonym_data_source)

    def __init__(self, ome_id, synonym, type, synonym_data_source_id):
        self.ome_id = ome_id
        self.synonym = synonym
        self.type = type
        self.synonym_data_source_id = synonym_data_source_id

class Publication(Base):
    __tablename__="publication"
    id = Column(Integer, Sequence('wids'), primary_key=True)
    pmid = Column(Integer)
    
    __table_arges__=(UniqueConstraint('pmid'),{})
    
    def __init_(self, pmid):
        self.pmid = pmid

class PublicationModel(Base):
    __tablename__="publication_model"
    model_id = Column(Integer, ForeignKey('model.id', ondelete='CASCADE'), primary_key=True)
    publication_id = Column(Integer, ForeignKey('publication.id', ondelete='CASCADE'), primary_key=True)
    __table_args__ = (UniqueConstraint('model_id', 'publication_id'),{})
    
    def __init__(self, model_id, publication_id):
        self.model_id = model_id
        self.publication_id = publication_id
        
class LinkOut(Base):
    __tablename__="link_out"
    id = Column(Integer, Sequence('wids'), primary_key=True) 
    external_id = Column(String)
    external_source = Column(String)
    ome_id = Column(Integer)
    type = Column(String)
        
class GenomeRegionMap(Base):
        __tablename__ = 'genome_region_map'

        genome_region_id_1 = Column(Integer, ForeignKey('genome_region.id'), primary_key=True)
        genome_region_id_2 = Column(Integer, ForeignKey('genome_region.id'), primary_key=True)
        distance = Column(Integer)

        __table_args__ = (UniqueConstraint('genome_region_id_1','genome_region_id_2'),{})


        def __repr__(self):
            return "GenomeRegionMap (%d <--> %d) distance:%d" % (self.genome_region_id_1, self.genome_region_id_2, self.distance)


        def __init__(self, genome_region_id_1, genome_region_id_2, distance):
            self.genome_region_id_1 = genome_region_id_1
            self.genome_region_id_2 = genome_region_id_2
            self.distance = distance


class _Session(_SA_Session):
    """an sqlalchemy session object to interact with the OME database

    This object can used to make queries against the ome database. For
    example, a query without using any ORM looks like this
    >>> session = Session()
    >>> session.execute("SELECT name from genes where bnum='b0001'").fetchone()
    (u'thrL',)
    Using the sqlalchemy ORM gives more descriptive objects. For example:
    >>> b0001 = session.query(Gene).filter(Gene.bnum=="b0001").first()
    >>> b0001.name
    u'thrL'
    Raw queries which return ORM objects are also possible:
    >>> sql_statement = "SELECT * from genes where bnum='b0001'"
    >>> b0001 = session.query(Gene).from_statement(sql_statement).first()
    >>> b0001.name
    u'thrL'

    The Session will automatically set the search_path to settings.schema
    """


    def __init__(self, *args, **kwargs):
        super(_Session, self).__init__(*args, **kwargs)
        #self.execute("set search_path to %s;" % (settings.schema))
        self.commit()
        self.get_or_create = MethodType(get_or_create, self)
        #self.search_by_synonym = MethodType(search_by_synonym, self)


    def __repr__(self):
        return "OME session %d" % (self.__hash__())


def get_or_create(session, class_type, **kwargs):
    """gets an object using filter_by on the unique kwargs. If no such object
    is found in the database, a new one will be created which satisfies
    these constraints. This is why every class that wants to use this
    method to be instantiated needs to have a UniqueConstraint defined.
    """

    for constraint in list(class_type.__table_args__):
        if constraint.__class__.__name__ == 'UniqueConstraint':
            unique_cols = constraint.columns.keys()

    inherited_result = True
    if '__mapper_args__' in class_type.__dict__ and 'inherits' in class_type.__mapper_args__:
        inherited_class_type = class_type.__mapper_args__['inherits']
        for constraint in list(inherited_class_type.__table_args__):
            if constraint.__class__.__name__ == 'UniqueConstraint':
				        inherited_unique_cols = constraint.columns.keys()

        inherited_result = session.query(inherited_class_type).filter_by(**{k: kwargs[k] for k in inherited_unique_cols}).first()


    result = session.query(class_type).filter_by(**{k: kwargs[k] for k in unique_cols}).first()

    if not result or not inherited_result:
        result = class_type(**kwargs)
        session.add(result)
        session.commit()

    return result


def update(session, object, **kwargs):
    """Ideally this would only search on the primary key columns so
    that an update could be made in one call. However, its not currently
    clear how to do that so necessary to pass in the actual object and
    update following a call to get_or_create() There is probably some
    way to do this with class_mapper but its hard right now
    """
    #result = session.query(class_type).filter_by(**kwargs).first()
    #result = session.query(class_type).filter_by(name=kwargs['name']).first()
    #if result is None: return

    for key,value in kwargs.iteritems():
        setattr(object,key,value)
    session.add(object)
    session.commit()

    return object


@contextmanager
def create_Session():
    session = Session()
    try:
        yield session
        session.commit()
    except:
        session.rollback()
        raise
    finally:
        print "close"
        session.close()



Session = sessionmaker(bind=engine, class_=_Session)


if __name__ == "__main__":
    session = Session()


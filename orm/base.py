"""Module to implement ORM to the ome database"""

from types import MethodType
from os import system

from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.session import Session as _SA_Session
from sqlalchemy import Table, MetaData, create_engine,Column, Integer, \
    String, Float, ForeignKey, and_, or_, not_, distinct, select
from sqlalchemy.ext.declarative import declarative_base
import trnlib.settings as settings
import pymongo

Base = declarative_base()

engine = create_engine("postgresql://%s:%s@%s/%s" %
    (settings.user, settings.password, settings.host, settings.dev_database))
metadata = MetaData(bind=engine, schema=settings.schema)

connection = pymongo.Connection()
bigg_database = connection.bigg_database


def make_table(table_name):
    """function to create a table with the default parameters"""
    return Table(table_name, metadata, autoload=True)


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
        self.execute("set search_path to %s;" % (settings.schema))
        self.commit()
        self.get_or_create = MethodType(get_or_create, self)
        #self.search_by_synonym = MethodType(search_by_synonym, self)

    def __repr__(self):
        return "OME session %d" % (self.__hash__())


def get_or_create(session, class_type, **kwargs):
    """gets an object using filter_by on the kwargs. If no such object
    is found in the database, a new one will be created which satisfies
    these constraints"""
    result = session.query(class_type).filter_by(**kwargs).first()
    if result is None:
        result = class_type()
        for key, value in kwargs.iteritems():
            setattr(result, key, value)
        session.add(result)
        session.commit()
    return result


def load_tables(module_name=None):
    if module_name == 'data':
        system("%s < %sdatabase/schemas/ome_data.sql > psql.log 2>&1" % (settings.psql_full, settings.trn_directory))
        connection = pymongo.Connection()
        genome_data = connection.bigg_database.genome_data
        genome_data.drop()
    elif module_name == 'models':
        system("%s < %sdatabase/schemas/ome_models.sql > psql.log 2>&1" % (settings.psql_full, settings.trn_directory))
        
        
Session = sessionmaker(bind=engine, class_=_Session)


if __name__ == "__main__":
    session = Session()
    
    
"""Module to implement ORM to the ome database"""

from types import MethodType
from os import system

from sqlalchemy.orm import sessionmaker, relationship, aliased
from sqlalchemy.orm.session import Session as _SA_Session
from sqlalchemy import Table, MetaData, create_engine,Column, Integer, \
    String, Float, ForeignKey, and_, or_, not_, distinct, select, func
from sqlalchemy.schema import UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property


engine = create_engine("postgresql://%s:%s@%s/%s" %
    ('sfederow', 'ome', 'zak.ucsd.edu', 'ome'))
Base = declarative_base(bind=engine)
metadata = MetaData(bind=engine, schema='escherichia')


        
def make_table(table_name):
    """function to create a table with the default parameters"""
    return Table(table_name, metadata, autoload=True)


def get_me_component(chemical_formula):
    components = []
    #remove H
    chemical_formula.pop('H', None)
    for component in Session().query(Component).all():
        formula = component.formula
        try: formula.pop('H', None)
        except: continue
        if formula == chemical_formula:
            components.append(component)
    return components


class Compartment(Base):
    __table__ = make_table('compartment')

    def __repr__(self):
        return "Compartment (#%d):  %s" % \
            (self.id, self.compartment_name)        

class Component(Base):
    __table__ = make_table('component')
    
    elements = relationship('ComponentChemicalFormula')
    
    @hybrid_property
    def formula(self):
        if self.elements:
            return {str(x.element):float(x.stoichiometry) for x in self.elements}
            #return ''.join([str(x.element)+str(x.stoichiometry) for x in self.elements])
        else:
            return None
        
    def __repr__(self):
        return "Component (#%d):  %s" % \
            (self.id, self.descriptive_name)


class ComponentChemicalFormula(Base):
    __table__ = make_table('component_chemical_formula')
   
    def __repr__(self):
        return "%d: %s %d" % \
            (self.component, self.element, self.stoichiometry) 

    
class ChemicalElement(Base):
    __table__ = make_table('chemical_element')
    
    def __repr__(self):
        return "%s: %s %5.2f" % \
            (self.id, self.element_name, self.atomic_mass)


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
        self.execute("set search_path to %s;" % ('escherichia'))
        self.commit()
        #self.get_or_create = MethodType(get_or_create, self)
        #self.search_by_synonym = MethodType(search_by_synonym, self)


    def __repr__(self):
        return "OME session %d" % (self.__hash__())


Session = sessionmaker(bind=engine, class_=_Session)

    
    
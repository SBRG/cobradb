from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, backref, Query
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.orm.collections import attribute_mapped_collection
from sqlalchemy import create_engine, Column, String, Float, ForeignKey, Boolean

from warnings import warn

engine = create_engine('sqlite:///:memory:', echo=False)
Base = declarative_base()
Session = sessionmaker(bind=engine)


class Reaction(Base):
    __tablename__ = "reactions"
    id = Column(String(200), primary_key=True)
    reversibility = Column(Boolean)
    # the reaction_metabolites are indexed by metabolite
    _reaction_metabolites = relationship("_ReactionMetabolites",
        backref="reaction",
        collection_class=attribute_mapped_collection("metabolite"))
    # the metabolites are a dict. The values will be the stoichiometry of
    # the reaction_metabolite. The values are already indexed by
    # the metabolite
    metabolites = association_proxy("_reaction_metabolites", "stoichiometry",
        creator=lambda metabolite, stoichiometry:  # key is metabolite
            _ReactionMetabolites(stoichiometry=stoichiometry,
                                metabolite=metabolite))

    def __init__(self, *args, **kwargs):
        if len(args) > 1:
            raise TypeError("Too many arguments supplied")
        if len(args) == 1:
            if "id" in kwargs:
                raise ValueError("id specified by both arg and kwarg")
            kwargs["id"] = args[0]
        Base.__init__(self, **kwargs)


    @property
    def reaction(self):
        """Generate a human readable reaction string."""
        reactant_dict = {}
        product_dict = {}
        def coefficient_to_string(number):
            if number == 1:
                return ""
            if number == int(number):
                return str(int(number))
            return "%.2f" % number
        for metabolite, coefficient in self.metabolites.items():
            id = metabolite.id
            if coefficient > 0:
                product_dict[id] = coefficient_to_string(coefficient)
            else:
                reactant_dict[id] = coefficient_to_string(abs(coefficient))
        reactant_string = " + ".join(['%s %s' % (coefficient_str, metabolite) for metabolite, coefficient_str in reactant_dict.items()])
        if not self.reversibility:
            arrow = ' -> '                
        else:
            arrow = ' <=> '
        product_string = " + ".join(['%s %s' % (coefficient_str, metabolite) for metabolite, coefficient_str in product_dict.items()])
        reaction_string = reactant_string + arrow + product_string
        return reaction_string

class Metabolite(Base):
    __tablename__ = "metabolite"
    id = Column(String(200), primary_key=True)
    name = Column(String(200))
    is_non_metabolite = Column(Boolean)
    _reaction_metabolites = relationship("_ReactionMetabolites",
        backref="metabolite")
    #formula = Column(String(400))
    reactions = relationship(Reaction,
        secondary=lambda: _ReactionMetabolites.__table__, viewonly=True)

    def __init__(self, *args, **kwargs):
        if len(args) > 1:
            raise TypeError("Too many arguments supplied")
        if len(args) == 1:
            if "id" in kwargs:
                raise ValueError("id specified by both arg and kwarg")
            kwargs["id"] = args[0]
        Base.__init__(self, **kwargs)
    
    def __str__(self):
        return str(self.id)

    def __repr__(self):
        return str(self)

class Chemical_Element(Base):
    __tablename__ = "chemical_element"
    id = Column(String(2),primary_key=True)
    element_name = Column(String(50))
    atomic_mass = Column(Float)
    
class Metabolite_Chemical_Formulat(Base):
    __tablename__ = "metabolite_chemical_formula"
    metabolite = relationship(Metabolite,primary_key=True)
    chemical_element = relationship(Chemical_Element,primary_key=True)
    stoichiometry = Column(Integer)


class Compartment(Base):
    __tablename__ = "compartment"
    id = Column(String(20), primary_key=True)


class _ReactionMetabolites(Base):
    __tablename__ = "reaction_matrix"
    reaction_id = Column(String(200),
        ForeignKey("reactions.id"), primary_key=True)
    metabolite_id = Column(String(200),
        ForeignKey("metabolites.id"), primary_key=True)
    stoichiometry = Column(Float)
    compartment_id = Column(String(20), ForeignKey("compartments.id"), primary_key=True)
    compartment = relationship(Compartment)
    def __repr__(self):
        return "(%s, %s) %f" % \
            (self.reaction_id, self.metabolite_id, self.stoichiometry)



Base.metadata.create_all(engine)

if __name__ == "__main__":
    c = Compartment()
    c.id = "cytosol"
    g6p = Metabolite()
    g6p.id = "g6p"
    f6p = Metabolite()
    f6p.id = "f6p"
    r = Reaction()
    r.id = "pgi"
    from IPython import embed; embed()

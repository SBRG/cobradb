"""Module to implement ORM to the ome database"""

from types import MethodType

from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy.orm.session import Session as _SA_Session
from sqlalchemy import Table, MetaData, create_engine, Column, Integer, \
    String, Float, ForeignKey, and_, or_, not_, distinct, select
from sqlalchemy.orm.collections import attribute_mapped_collection
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.sql.expression import text
import settings
import simplejson


engine = create_engine("postgresql://%s:%s@%s/%s" %
    (settings.user, settings.password, settings.host, settings.database))
Base = declarative_base()
metadata = MetaData(bind=engine, schema=settings.schema)

def _make_table(table_name):
    """function to create a table with the default parameters"""
    return Table(table_name, metadata, autoload=True)

# tables and classes which are used for linking only
_gene_wid_protein_wid = _make_table("gene_wid_protein_wid")
_citation_list = _make_table("citation_list")

_operon_components = _make_table("operon_components")
_tu_components = _make_table("tu_components")
_knockouts = _make_table("knockouts")
_pathway_genes = _make_table("pathway_genes")
_reactants = _make_table("reactants")
_products = _make_table("products")
class _RuleRegulation(Base):
    __table__ = _make_table("rule_regulation")
class _RuleComponent(Base):
    __table__ = _make_table("rule_components")

# class to make carbon_source, nitrogen_source etc. queryable if the object
# is related to a condition
class ExpandedConditions:
    @hybrid_property
    def carbon_source(self):
        return self.condition.carbon_source
    @carbon_source.expression
    def carbon_source(cls):
        return select([Condition.carbon_source]).\
            where(Condition.wid == cls.condition_wid).label("carbon_source")
    @hybrid_property
    def nitrogen_source(self):
        return self.condition.nitrogen_source
    @nitrogen_source.expression
    def nitrogen_source(cls):
        return select([Condition.nitrogen_source]).\
            where(Condition.wid == cls.condition_wid).label("nitrogen_source")
    @hybrid_property
    def eacceptor(self):
        return self.condition.eacceptor
    @eacceptor.expression
    def eacceptor(cls):
        return select([Condition.eacceptor]).\
            where(Condition.wid == cls.condition_wid).label("eacceptor")

# class to give a sequence to objects with a leftpos, rightpos, and strand
class WithSequence:
    @hybrid_property
    def sequence(self):
        t = text("SELECT base FROM genome WHERE position >= :pos1 AND position <= :pos2")
        result = engine.execute(t, pos1=self.leftpos, pos2=self.rightpos)
        if self.strand == '+': return str(''.join([b.base for b in result]))
        else: return reverse_complement(str(''.join([b.base for b in result])))

def _create_citation_relation(__table__):
    """links an object to its citations"""
    return relationship("Citation", secondary=_citation_list,
        primaryjoin=__table__.c.wid == _citation_list.c.table_wid,
        secondaryjoin=_citation_list.c.citation_wid == Citation.wid,
        foreign_keys=[_citation_list.c.table_wid,
            _citation_list.c.citation_wid])

def _create_regulator_rule_relation(object_name):
    """links a regulated object (i.e. TU) to any regulatory rules"""
    return relationship(Rule,
        secondary=_RuleRegulation.__table__, uselist=True,
        primaryjoin="%s.wid == _RuleRegulation.regulates_wid" % (object_name),
        secondaryjoin=_RuleRegulation.rule_wid == Rule.wid,
        foreign_keys=[_RuleRegulation.regulates_wid,
            _RuleRegulation.rule_wid],
        backref="regulates_%ss" % object_name.lower())

def _create_current_regulator_rule_relation(table_name):
    """links a regulated object (i.e. TU) to its *current* regulatory rule"""
    return relationship(Rule,
        secondary=_RuleRegulation.__table__, uselist=False,
        primaryjoin="%s.wid == _RuleRegulation.regulates_wid" % (table_name),
        secondaryjoin="and_(_RuleRegulation.rule_wid == Rule.wid," + \
                          " Rule.replaced_by_wid == None)",
        foreign_keys=[_RuleRegulation.regulates_wid,
            _RuleRegulation.rule_wid],
        backref="currently_regulates_%s" % table_name, viewonly=True)

def _create_rule_component_relation(__table__):
    """links an object  to a rule it may participate in"""
    # TODO - use this!
    return relationship(Rule, secondary=_rule_components, uselist=True,
        primaryjoin=__table__.c.wid == _rule_components.c.component_wid,
        secondaryjoin=_rule_components.c.rule_wid == Rule.wid,
        foreign_keys=[_rule_components.c.component_wid,
            _rule_components.c.rule_wid],
        backref=backref("component_%s" % __table__.name))

def _create_synonym_list(__table__):
    """creates a list of synonyms by looking in wid2otherid"""
    return relationship(wid2otherid,
        primaryjoin=__table__.c.wid == wid2otherid.wid,
        foreign_keys=[wid2otherid.wid])

def _create_reaction_relation(__table__, secondary, other, backref_name=None):
    """link a reaction to either its reactants or products"""
    if backref_name is None:
        backref_name = "%s" % __table__.name
    return relationship(other, secondary=secondary,
        backref=backref(backref_name),
        primaryjoin=__table__.c.wid == secondary.c.reaction_wid,
        secondaryjoin=secondary.c.other_wid == other.wid,
        foreign_keys=[secondary.c.reaction_wid,
            secondary.c.other_wid])


#These functions should probably get moved somewhere

def complement(seq): 
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    complseq = [complement[base] for base in seq] 
    return complseq

def reverse_complement(seq): 
    seq = list(seq) 
    seq.reverse() 
    return ''.join(complement(seq))
    
class Citation(Base):
    __table__ = _make_table("citations")
    def __repr__(self):
        if self.date is None: year = ""
        else: year = str(self.date.year)
        repr = u"Citation: %s %s (PMID %d)" % \
            (self.first_author, year, self.pmid)
        return repr.encode('ascii', 'replace')

class wid2otherid(Base):
    __table__ = _make_table("wid2otherid")
    dataset = relationship("Dataset")
    def __repr__(self):
        return "%s in (%s)" % (self.otherid, str(self.dataset))

class Protein(Base):
    __table__ = _make_table('proteins')
    genes = relationship("Gene", secondary=_gene_wid_protein_wid,
        backref="proteins")
    _rule_components = relationship("_RuleComponent",
        primaryjoin="Protein.wid == _RuleComponent.component_wid",
        backref="component_protein",
        foreign_keys=[_RuleComponent.__table__.c.component_wid])
    associated_rules = association_proxy("_rule_components",
        "rule",
        creator=lambda new_rule:
            _RuleComponent(type="protein", rule=new_rule))
    # associated_rules = _create_rule_component_relation(__table__)
    synonyms = _create_synonym_list(__table__)
    citations = _create_citation_relation(__table__)
    def __repr__(self):
        return "Protein: %s" % (self.name)

class Rule(Base):
    """A regulatory rule for an entity such as a TU

    All active rules have not yet been replaced, so replaced_by = None
    In order to replace a rule, create a new rule, and set replaced_by to
    the new rule, then set what the new rule should regulate. No rule that
    has been replaced will show up as a current rule."""
    __table__ = _make_table("rules")
    citations = _create_citation_relation(__table__)
    replaced_by = relationship("Rule")
    _rule_components = relationship("_RuleComponent", backref="rule")
    component_proteins = association_proxy("_rule_components",
        "component_protein",
        creator=lambda new_protein:
            _RuleComponent(type="protein", component_protein=new_protein))
    def __repr__(self):
        return "Rule: %s" % (self.rule)

class GeneTag(Base):
    __table__ = _make_table("tags")
    strain = relationship("Strain")
    gene = relationship("Gene")

class Strain(Base):
    __table__ = _make_table("strains")
    knockouts = relationship("Gene", secondary=_knockouts)
    tags = relationship("Gene", secondary=GeneTag.__table__)
    def __repr__(self):
        return "Strain: %s" % (self.name)

class EnzymaticReaction(Base):
    __table__ = _make_table("enzymatic_reactions")
    synonyms = _create_synonym_list(__table__)

class MotifInfo(Base):
    __table__ = _make_table("motif_info")
    @hybrid_property
    def sequence(self):
        return self.binding_site.sequence

class BindingSite(Base, WithSequence):
    __table__ = _make_table("binding_sites")
    motif_info = relationship("MotifInfo", backref="binding_site")
    # super broken way of determining the dataset...
    dataset = relationship("Dataset", secondary=wid2otherid.__table__,
        primaryjoin="BindingSite.wid == wid2otherid.wid",
        secondaryjoin="wid2otherid.dataset_wid == Dataset.wid",
        foreign_keys=[wid2otherid.wid, wid2otherid.dataset_wid], uselist=False)
    def __repr__(self):
        return "BindingSite: (%s)%d-%d" % \
            (self.strand, self.leftpos, self.rightpos)

class TSS(Base):
    __table__ = _make_table("tss")
    def __repr__(self):
        return "TSS: %s(%d)" % (self.name, self.position)

class Gene(Base, WithSequence):
    __table__ = _make_table("genes")
    rules = _create_regulator_rule_relation("Gene")
    current_rule = _create_current_regulator_rule_relation("Gene")
    citations = _create_citation_relation(__table__)
    synonyms = _create_synonym_list(__table__)
    array_data = relationship("ArrayData",
        secondary=_make_table("array_mapping"),
        backref=backref("gene", uselist=False))

    def __repr__(self):
        display_name = self.name
        if display_name.startswith("y"):
            for synonym in self.synonyms:
                if synonym.dataset.name == "ecocyc":
                    syn = synonym.otherid
                    if not syn[:1].isdecimal():
                        display_name += "/" + syn
                        break
        return "Gene: %s (%s)" % (display_name, self.bnum)
    def _repr_html_(self):
        url = "http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=" + self.bnum
        return "<a href='%s' target=\"_blank\">%s</a>" % (url, repr(self)) 
    
class TUComponents(Base):
    __table__ = _make_table("tu_components")

class TU(Base, WithSequence):
    __table__ = _make_table("tu")
    rules = _create_regulator_rule_relation("TU")
    current_rule = _create_current_regulator_rule_relation("TU")
    citations = _create_citation_relation(__table__)
    components = relationship(TUComponents, backref="TU")
    genes = relationship("Gene", secondary=_tu_components, backref="TU",
        primaryjoin=_tu_components.c.tu_wid == __table__.c.wid,
        secondaryjoin=_tu_components.c.other_wid == Gene.wid,
        foreign_keys=[_tu_components.c.tu_wid,
                      _tu_components.c.other_wid])
    tss = relationship("TSS", secondary=_tu_components,
        backref=backref("TU", uselist=False), uselist=False,
        primaryjoin=_tu_components.c.tu_wid == __table__.c.wid,
        secondaryjoin=_tu_components.c.other_wid == TSS.wid,
        foreign_keys=[_tu_components.c.tu_wid,
                      _tu_components.c.other_wid])

    def __repr__(self):
        return "TU: (name:%s, wid:%d)" % (self.name, self.wid)

class ComplexComponents(Base):
    __table__ = _make_table("complex_components")

class Chemical(Base):
    __table__ = _make_table("chemicals")

def _make_component_relationship(obj_name):
    return relationship(obj_name,
        secondary=ComplexComponents.__table__,
        primaryjoin="ComplexComponents.complex_wid == Complex.wid",
        secondaryjoin="ComplexComponents.other_wid == %s.wid" % (obj_name),
        foreign_keys=[ComplexComponents.other_wid,
                      ComplexComponents.complex_wid])

class Complex(Base):
    __table__ = _make_table("complex")
    components = relationship(ComplexComponents, backref="complex")
    protein_components = _make_component_relationship("Protein")
    chemical_components = _make_component_relationship("Chemical")
    binding_site_components = _make_component_relationship("BindingSite")
    complex_components = _make_component_relationship("Complex")
    citations = _create_citation_relation(__table__)

class Operon(Base, WithSequence):
    __table__ = _make_table('operons')
    tus = relationship("TU", secondary=_operon_components, backref="operons",
        primaryjoin=_operon_components.c.operon_wid == __table__.c.wid,
        secondaryjoin=_operon_components.c.other_wid == TU.wid,
        foreign_keys=[_operon_components.c.operon_wid,
                      _operon_components.c.other_wid])
    binding_sites = relationship("BindingSite",
        secondary=_operon_components, backref="operons",
        primaryjoin=_operon_components.c.operon_wid == __table__.c.wid,
        secondaryjoin=_operon_components.c.other_wid == BindingSite.wid,
        foreign_keys=[_operon_components.c.operon_wid,
                      _operon_components.c.other_wid])
    def __repr__(self):
        return "Operon %s (%s:%d-%d)" % (self.name, self.strand,
            self.leftpos, self.rightpos)

class Dataset(Base):
    __table__ = _make_table('datasets')
    citations = _create_citation_relation(__table__)
    def __repr__(self):
        return "Dataset %s (#%d)" % (self.name, self.wid)

class Condition(Base):
    __table__ = _make_table("conditions")
    def __repr__(self):
        return "Condition: C:%s, N:%s, e:%s" % \
            (self.carbon_source, self.nitrogen_source, self.eacceptor)

class ExperimentSet(Base):
    __table__ = _make_table("experiment_set")
    dataset = relationship(Dataset)
    def __repr__(self):
        return "ExperimentSet (#%d, %s):  %s" % \
            (self.wid, self.type, self.dataset)

class GrowthRateExperiment(Base):
    __table__ = _make_table("growth_rate_experiment")
    dataset = relationship(Dataset)
    condition = relationship(Condition)
    strain = relationship(Strain)
    experiment_set = relationship(ExperimentSet)

class GFFExperiment(Base, ExpandedConditions):
    __table__ = _make_table("gff_experiment")
    dataset = relationship(Dataset)
    condition = relationship(Condition)
    strain = relationship(Strain)
    experiment_set = relationship(ExperimentSet, backref="gff_experiments")
    def __repr__(self):
        return "GFFExperiment: %s on %s in %s" % \
            (self.strain, self.condition, self.dataset)
    def _repr_json_(self):
        encoder = simplejson.JSONEncoder()
        return encoder.encode({
            "strain_name": self.strain.name,
            "dataset_name": self.dataset.name,
            "condition_wid": self.condition.wid,
            "carbon_source": self.condition.carbon_source})

class GFFData(Base):
    __table__ = _make_table("gff_data")
    experiment = relationship(GFFExperiment)
    experiment_set = association_proxy("experiment", "experiment_set")
    def __repr__(self):
        return "GFFData: value of %.3f %s%d-%d" % \
            (self.value, self.strand, self.leftpos, self.rightpos)

class TranscriptionalReaction(Base):
    __table__ = _make_table("transcriptional_reactions")
    citations = _create_citation_relation(__table__)
    product_tu = _create_reaction_relation(__table__, _products, TU)

class Pathway(Base):
    __table__ = _make_table("pathways")
    genes = relationship(Gene, secondary=_pathway_genes, backref="pathways")
    citations = _create_citation_relation(__table__)
    def __repr__(self):
        return "Pathway: %s" % (self.name)

class ArrayExperiment(Base, ExpandedConditions):
    __table__ = _make_table("array_experiment")
    condition = relationship(Condition)
    strain = relationship(Strain)
    experiment_set = relationship(ExperimentSet, backref="array_experiments")
    def __repr__(self):
        return "ArrayExperiment: %s on %s in %s" % \
            (self.strain, self.condition, self.dataset)

class ArrayData(Base):
    __table__ = _make_table("array_data")
    experiment = relationship(ArrayExperiment)
    strain = association_proxy("experiment", "strain")
    condition = association_proxy("experiment", "condition")
    def __repr__(self):
        return "ArrayData: value %.2f on %s for %s" % (self.value,
            self.condition, self.strain)

class ArrayAnalysis(Base):
    __table__ = _make_table("array_analysis")
    gene = relationship("Gene", backref="array_analysis")
    experiment1 = relationship("ArrayExperiment",
        primaryjoin="ArrayAnalysis.experiment_wid_1 == ArrayExperiment.wid")
    experiment2 = relationship("ArrayExperiment",
        primaryjoin="ArrayAnalysis.experiment_wid_2 == ArrayExperiment.wid")
    condition1 = association_proxy("experiment1", "condition")
    condition2 = association_proxy("experiment2", "condition")
    strain1 = association_proxy("experiment1", "strain")
    strain2 = association_proxy("experiment2", "strain")
    @property
    def _strain_str(self):
        if self.strain1 == self.strain2:
            return str(self.strain1)
        else:
            strain_str = "%s-->%s" % (str(self.strain1), str(self.strain2))
            strain_str = strain_str.replace("Strain: ", "")
            return "Strain: " + strain_str
    @property
    def _condition_str(self):
        if self.condition1 == self.condition2:
            condition_str = str(self.condition1)
        else:
            condition_str = "Condition: "
            attrs = ["carbon_source", "nitrogen_source", "eacceptor"]
            abbrevs = ["C", "N", "e"]  # abbreviations
            for attr, abbrev in zip(attrs, abbrevs):
                attr1 = getattr(self.condition1, attr)
                attr2 = getattr(self.condition2, attr)
                if attr1 == attr2:
                    condition_str += "%s:%s" % (abbrev, attr1)
                else:
                    condition_str += "%s:%s-->%s" % (abbrev, attr1, attr2)
                condition_str += ", "
            condition_str = condition_str.rstrip(", ")
        return condition_str
    
    def __repr__(self):
        return "ArrayAnalysis (%s): %f for %s %s on %s" % \
            (self.type, self.value, self.gene, self._strain_str, self._condition_str)
    def _repr_html_(self):
        return "ArrayAnalysis (%s): %f for %s %s on %s" % \
            (self.type, self.value, self.gene._repr_html_(), self._strain_str,
             self._condition_str.replace("-->", "&rarr;"))
    
    
# views
class AllAnalysis(Base):
    __table__ = Table("all_analysis", metadata,
        Column('bnum', String, ForeignKey("genes.bnum"), primary_key=True),
        Column('target1', String, primary_key=True),
        Column('target2', String, primary_key=True),
        Column('c1_wid', Integer, ForeignKey("conditions.wid"),
            primary_key=True),
        Column('c2_wid', Integer, ForeignKey("conditions.wid"),
            primary_key=True),
        Column('aa_fdr_wid', Integer, ForeignKey("array_analysis.wid")),
        Column('aa_fc_wid', Integer, ForeignKey("array_analysis.wid")),
        autoload=True)
    gene = relationship(Gene, backref="all_array_analysis", viewonly=True)
    array_analysis_fold_change = relationship("ArrayAnalysis", viewonly=True,
        primaryjoin="ArrayAnalysis.wid == AllAnalysis.aa_fc_wid")
    array_analysis_fdr = relationship("ArrayAnalysis", viewonly=True,
        primaryjoin="ArrayAnalysis.wid == AllAnalysis.aa_fdr_wid")
    condition1 = relationship(Condition, viewonly=True,
        primaryjoin="AllAnalysis.c1_wid == Condition.wid")
    condition2 = relationship(Condition, viewonly=True,
        primaryjoin="AllAnalysis.c2_wid == Condition.wid")
    strain1 = association_proxy("array_analysis_fold_change", "strain1")
    strain2 = association_proxy("array_analysis_fold_change", "strain2")
    def __repr__(self):
        s = repr(self.array_analysis_fold_change)
        s = s[s.find("for"):]
        s = "AllAnalysis (fdr: %f, fold change %f): %s" % \
            (self.fdr, self.fold_change, s)
        return s
    def _repr_html_(self):
        s = self.array_analysis_fold_change._repr_html_()
        s = s[s.find("for"):]
        color = "99FF66" if self.fold_change < 0 else "FF9999"
        s = 'AllAnalysis (fdr: %f, <span style="background-color: #%s">fold change %f</span>): %s' % \
            (self.fdr, color, self.fold_change, s)
        return s

class AveragedExpressionData(Base):
   __table__ = Table("averaged_expression_data", metadata,
       Column('gene_wid', Integer, ForeignKey("genes.wid"), primary_key=True),
       Column('condition_wid', Integer, ForeignKey("conditions.wid"), primary_key=True),
       Column('exp_set_wid', Integer, ForeignKey('experiment_set.wid'), primary_key=True),
       #Column('strain_wid', Integer, ForeignKey("strains.wid"), primary_key=True),
       autoload=True)
   gene = relationship(Gene, viewonly=True)
   condition = relationship(Condition, viewonly=True)
   #strain = relationship(Strain, viewonly=True)
   experiment_set = relationship(ExperimentSet, viewonly=True)


class ExpressionData(Base):
    __table__ = Table("expression_data", metadata,
        Column('gene_wid', Integer, ForeignKey("genes.wid"), primary_key=True),
        Column('bnum', String),
        Column('condition_wid', Integer, ForeignKey("conditions.wid"), primary_key=True),
        Column('replicate', Integer, primary_key=True),
        Column('exp_set_wid', Integer, ForeignKey('experiment_set.wid'), primary_key=True),
        Column('value', Float),
        Column('strain_wid', Integer, ForeignKey("strains.wid"), primary_key=True),
        autoload=True)
    gene = relationship(Gene, viewonly=True)
    condition = relationship(Condition, viewonly=True)
    strain = relationship(Strain, viewonly=True)
    experiment_set = relationship(ExperimentSet, viewonly=True)

       
#class ExpressionComparison(Base):
#    __table__ = Table("expression_comparison", metadata,
#        Column('bnum', String, ForeignKey("genes.bnum"), primary_key=True),
#        Column('strain1', String, primary_key=True),
#        Column('strain2', String, primary_key=True),
#        Column('c1_wid', Integer, ForeignKey("conditions.wid"),
#            primary_key=True),
#        Column('c2_wid', Integer, ForeignKey("conditions.wid"),
#            primary_key=True),
#        autoload=True)

class BindingComplex(Base):
    __table__ = Table("binding_complex", metadata,
        Column('complex_wid', Integer, ForeignKey("complex.wid"), primary_key=True),
        Column('binding_site_wid', Integer, ForeignKey("binding_sites.wid")),
        Column('dataset_wid', Integer, ForeignKey('datasets.wid')),
        Column('protein_wid', Integer, ForeignKey('proteins.wid')),
        autoload=True)
    complex = relationship("Complex")
    binding_site = relationship("BindingSite")
    dataset = relationship("Dataset")
    protein = relationship("Protein")
    def __repr__(self):
        return "Binding complex: %s" % self.complex_name

#class ChipPeak(Base):
    #__table__ = Table("chip_peaks", metadata,
        #Column('binding_site_wid', Integer, ForeignKey("binding_sites.wid"), primary_key=True),
        #Column('experiment_wid', Integer, ForeignKey("gff_experiment.wid")),
        #Column('condition_wid', Integer, ForeignKey("conditions.wid")),
        #autoload=True)
    #condition = relationship(Condition, viewonly=True)
    #binding_site = relationship(BindingSite, viewonly=True)

class ChipPeakData(Base):
    __table__ = Table("chip_peak_data", metadata,
        Column('binding_site_wid', Integer, ForeignKey("binding_sites.wid"), primary_key=True),
        Column('experiment_wid', Integer, ForeignKey("gff_experiment.wid")),
        Column('condition_wid', Integer, ForeignKey("conditions.wid")),
        autoload=True)
    condition = relationship(Condition, viewonly=True)
    binding_site = relationship(BindingSite, viewonly=True)

class ChipPeakGene(Base, ExpandedConditions):
    __table__ = Table("chip_peak_gene", metadata,
        Column('binding_site_wid', Integer, ForeignKey("binding_sites.wid"), primary_key=True),
        Column('gene_wid', Integer, ForeignKey("genes.wid"), primary_key=True),
        Column('condition_wid', Integer, ForeignKey("conditions.wid"), primary_key=True),
        autoload=True)
    gene = relationship(Gene, viewonly=True, backref="chip_peaks")
    condition = relationship(Condition, viewonly=True)
    binding_site = relationship(BindingSite, viewonly=True)
    def __repr__(self):
        return "ChipPeakGene (%s): (%s, %s, %s)" % \
            (self.target, str(self.gene), str(self.condition), str(self.binding_site))

class ChipPeakTUGene(Base):
    __table__ = Table("chip_peak_tu_gene", metadata,
        Column('binding_site_wid', Integer, ForeignKey("binding_sites.wid"), primary_key=True),
        Column('gene_wid', Integer, ForeignKey("genes.wid"), primary_key=True),
        Column('condition_wid', Integer, ForeignKey("conditions.wid"), primary_key=True),
        Column('tu_wid', Integer, ForeignKey("tu.wid"), primary_key=True),
        autoload=True)
    gene = relationship(Gene, viewonly=True)
    tu = relationship(TU, viewonly=True)
    condition = relationship(Condition, viewonly=True)
    binding_site = relationship(BindingSite, viewonly=True)

class ChipPeakGeneExpression(Base):
    __table__ = Table("chip_peak_gene_expression", metadata,
        Column('binding_site_wid', Integer, ForeignKey("binding_sites.wid"), primary_key=True),
        Column('gene_wid', Integer, ForeignKey("genes.wid"), primary_key=True),
        Column('condition_wid', Integer, ForeignKey("conditions.wid"), primary_key=True),
        Column('fdr', Float),
        Column('fold_change', Float),
        autoload=True)
    gene = relationship(Gene, backref="chip_peak_gene_expression", viewonly=True)
    condition = relationship(Condition, viewonly=True)
    binding_site = relationship(BindingSite, viewonly=True)
    pathways = association_proxy("gene", "pathways")

    @hybrid_property
    def effect(self):
        # if the knockout is lower than wt, then it is an activator (without
        # the activator, expression is lower)
        # knockout < wt ==> activator
        # knockout > wt ==> repressor
        if self.fold_change == 0:
            return "error - not differentially expressed"
        if self.target1 == "wt":
            # fold change > 0 ==> wt > knockout ==> knockout < wt
            if self.fold_change > 0:
                return "activator"
            else:
                return "repressor"
        if self.target2 == "wt":
            # fold change < 0 ==> wt > knockout
            if self.fold_change < 0:
                return "activator"
            else:
                return "repressor"
    def __repr__(self):
        if self.fdr < 0.01:
            fdr_str = "%.2e" % self.fdr
        else:
            fdr_str = "%0.2f" % self.fdr
        return "ChipPeakGeneExpression: %s (%s) for %s on %s (fold change %.2f, fdr=%s)" % \
            (self.target, self.effect, str(self.gene), str(self.condition), self.fold_change, fdr_str)

    def _repr_html_(self):
        if self.fdr < 0.01:
            fdr_str = "%.2e" % self.fdr
        else:
            fdr_str = "%0.2f" % self.fdr
        return "ChipPeakGeneExpression: %s (%s) for %s on %s (fold change %.2f, fdr=%s)" % \
            (self.target, self.effect, self.gene._repr_html_(), str(self.condition), self.fold_change, fdr_str)

#class NimbleScanPeaks(Base):
#    __table__ = Table("nimblescan_peaks", metadata,
#        Column("wid", Integer, primary_key=True), autoload=True)
#    def __repr__(self):
#        return "NimbleScan peak: target: '%s', leftpos: %i, rightpos: %i" % \
#            (self.target, self.leftpos, self.rightpos)

class RegulatoryNetwork(Base):
    __table__ = Table("regulatory_network", metadata,
        Column("reg_bnum", String, ForeignKey("genes.bnum"), primary_key=True),
        Column("regd_bnum", String, ForeignKey("genes.bnum"), primary_key=True),
        autoload=True)
    reg = relationship(Gene, primaryjoin="RegulatoryNetwork.reg_bnum == Gene.bnum")
    regd = relationship(Gene, primaryjoin="RegulatoryNetwork.regd_bnum == Gene.bnum")
    def __repr__(self):
        if self.direction == "+":
            arrow = ">"
        elif self.direction == "-":
            arrow = "|"
        elif self.direction == "+-":
            arrow = "?"
        else:
            arrow = self.direction
        return "RegulonDB %s --%s %s (%s)" % \
            (self.reg.name, arrow, self.regd.name, self.quality)

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
        self.search_by_synonym = MethodType(search_by_synonym, self)

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

def get_or_create_condition(session, carbon_source, nitrogen_source,
                            eacceptor, temperature=37, other=None):
    """normalize names in condition creation"""
    if carbon_source.lower() == "glc": carbon_source = "glucose"
    if nitrogen_source == "NH4CL": nitrogen_source = "NH4Cl"
    if nitrogen_source == "Adenine": nitrogen_source = "adenine"
    if eacceptor == "aerobic": eacceptor = "O2"
    return get_or_create(session, Condition, carbon_source=carbon_source,
        nitrogen_source=nitrogen_source, eacceptor=eacceptor,
        temperature=temperature, other=other)

def search_by_synonym(session, class_type, synonym):
    """performs a query to find an object of class_type with the given
    synonym in its synonyms field (which was comes from wid2otherid)"""
    return session.query(class_type).filter(class_type.synonyms.any(
        wid2otherid.otherid.ilike(synonym)))


Session = sessionmaker(bind=engine, class_=_Session)

if __name__ == "__main__":
    session = Session()
    gene = session.query(Gene).filter(Gene.name.ilike("arca")).first()
    print gene
    print "Array Data"
    for ad in gene.array_data: print ad
    print "Array Analysis"
    for aa in gene.array_analysis: print aa

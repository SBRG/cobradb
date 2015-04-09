from cobra.io import sbml
import pytest
import numpy
from os.path import abspath, dirname, join
def test_different_sbml():
    sbml_file_name = "iJO1366"
    database_model = sbml.create_cobra_model_from_sbml_file(join("/Users/dbuser/bigg_stage/bigg2/bigg2/static/model_dumps", sbml_file_name + ".xml"))
    published_model = sbml.create_cobra_model_from_sbml_file(join("/Users/dbuser/bigg_stage/bigg2/bigg2/static/published_models",sbml_file_name + "_published.xml"))
    database_model.reactions.get_by_id("EX_glc_e").lower_bound = -10
    solution1 = database_model.optimize()
    solution2 = published_model.optimize()
    assert solution1.f == solution2.f
    assert len(database_model.reactions) == len(published_model.reactions)
    assert len(database_model.metabolites) == len(published_model.metabolites)
    assert len(database_model.genes) == len(published_model.genes)
"""    
def test_all_models():
    for model in session.query(Model).all():
        mcc = session.query(ModelCompartmentalizedComponent).join(Component).filter(Model.id == model.id).filter(Component.bigg_id =="pyr").first()
        if mcc == None:
            logging.warning(model.id + " has no pyr")
"""
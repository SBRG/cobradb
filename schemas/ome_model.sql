SET search_path TO ecoli;

-- components

DROP TABLE IF EXISTS model CASCADE;
DROP TABLE IF EXISTS reaction CASCADE;
DROP TABLE IF EXISTS model_reaction CASCADE;
DROP TABLE IF EXISTS compartment CASCADE;
DROP TABLE IF EXISTS compartmentalized_component CASCADE;
DROP TABLE IF EXISTS reaction_matrix CASCADE;



CREATE TABLE model (
    ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
    dataset_ID INT, FOREIGN KEY (dataset_ID) REFERENCES dataset(ID) ON DELETE CASCADE,
    name VARCHAR(100));


CREATE TABLE reaction (
    ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
    dataset_ID INT, FOREIGN KEY (dataset_ID) REFERENCES dataset(ID) ON DELETE CASCADE,
    abbreviation VARCHAR(100) NOT NULL UNIQUE,
    name VARCHAR(300));


CREATE TABLE model_reaction (
    model_ID INT, FOREIGN KEY (model_ID) REFERENCES model(ID) ON DELETE CASCADE,
    reaction_ID INT, FOREIGN KEY (reaction_ID) REFERENCES reaction(ID) ON DELETE CASCADE);


CREATE TABLE compartment (
    ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
    name VARCHAR(45) NOT NULL UNIQUE,
    abbreviation VARCHAR(45));


CREATE TABLE compartmentalized_component (
    ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
    component_ID INT, FOREIGN KEY (component_ID) REFERENCES component(ID) ON DELETE CASCADE,
    compartment_ID INT, FOREIGN KEY (compartment_ID) REFERENCES compartment(ID) ON DELETE CASCADE);


CREATE TABLE reaction_matrix (
  reaction_ID INT, FOREIGN KEY (reaction_ID) REFERENCES reaction(ID) ON DELETE CASCADE, 
  component_ID INT, FOREIGN KEY (component_ID) REFERENCES compartmentalized_component(ID) ON DELETE CASCADE,
  stoichiometry NUMERIC not null ,
  PRIMARY KEY (reaction_ID, component_ID));


    
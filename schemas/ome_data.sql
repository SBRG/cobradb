SET search_path TO ecoli;


DROP TABLE IF EXISTS data_source CASCADE;
DROP TABLE IF EXISTS environment_types CASCADE;
DROP TABLE IF EXISTS environments CASCADE;
DROP TABLE IF EXISTS in_vivo_environments CASCADE;
DROP TABLE IF EXISTS in_silico_environments CASCADE;
DROP TABLE IF EXISTS strains CASCADE;
DROP TABLE IF EXISTS protocols CASCADE;
DROP TABLE IF EXISTS datasets CASCADE;
DROP TABLE IF EXISTS dataset_types CASCADE;
DROP TABLE IF EXISTS array_experiment CASCADE;
DROP TABLE IF EXISTS chip_experiment CASCADE;
DROP TABLE IF EXISTS rna_seq_experiment CASCADE;
DROP TABLE IF EXISTS metabolomics_experiment CASCADE;
DROP TABLE IF EXISTS tss_experiment CASCADE;
DROP TABLE IF EXISTS gcrma_normalization CASCADE;
DROP TABLE IF EXISTS fpkm CASCADE;


CREATE TABLE data_source ( --TODO insert more information
	ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
	name VARCHAR(100),
	lab VARCHAR(100),
	institution VARCHAR(100));
	
	
CREATE TABLE environment_types (
	ID INT NOT NULL PRIMARY KEY,
	name VARCHAR(10));	

insert into environment_types (ID, name) values
  (1, 'in vivo'),
  (2, 'in silico');
	
	
CREATE TABLE environments (
	ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
	env_type_ID INT, FOREIGN KEY (env_type_ID) REFERENCES environment_types(ID) ON DELETE CASCADE, 
	name VARCHAR(100) UNIQUE);
	--PRIMARY KEY (ID, env_type_ID));


CREATE TABLE in_vivo_environments (
	environment_ID INT NOT NULL PRIMARY KEY,
	env_type_ID INT CHECK (env_type_ID = 1),
	carbon_source VARCHAR(100),
	nitrogen_source VARCHAR(100),
	electron_acceptor VARCHAR(10),
	temperature DECIMAL,
	FOREIGN KEY (environment_ID) REFERENCES environments (ID));


CREATE TABLE in_silico_environments (
	environment_ID INT NOT NULL PRIMARY KEY,
	env_type_ID INT CHECK (env_type_ID = 2),	
	FOREIGN KEY (environment_ID) REFERENCES environments (ID));
	
	
CREATE TABLE strains ( --TODO link out to graph DB
	ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
	name VARCHAR(20));


CREATE TABLE protocols (
	ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
	name VARCHAR(100),
	location VARCHAR(100));


CREATE TABLE dataset_types (
	ID INT NOT NULL PRIMARY KEY,
	name VARCHAR(100));
	
	
insert into dataset_types (ID, name) values
  (1, 'Array_experiment'),
  (2, 'RNAseq_experiment'),
  (3, 'ChIP_experiment'),
  (4, 'Metabolomics_experiment'),
  (5, 'TSS_experiment'),
  (6, 'gcrma normalization'),
  (7, 'FPKM'),
  (8, 'ME processing'),
  (9, 'metabolomic normalization'),
  (10, 'peak calling');

  
  	
CREATE TABLE datasets (
	ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
	dataset_type_ID INT, FOREIGN KEY (dataset_type_ID) REFERENCES dataset_types(ID) ON DELETE CASCADE,
	name VARCHAR(100),
	data_source_ID INT, FOREIGN KEY (data_source_ID) REFERENCES data_source(ID) ON DELETE CASCADE,
	dataset_IDs INT[]);


CREATE TABLE array_experiments (
	dataset_ID INT NOT NULL PRIMARY KEY,
	dataset_type_ID INT CHECK (dataset_type_ID = 1),
	environment_ID INT, FOREIGN KEY (environment_ID) REFERENCES environments(ID) ON DELETE CASCADE,
	strain_ID INT, FOREIGN KEY (strain_ID) REFERENCES strains(ID) ON DELETE CASCADE,
	platform VARCHAR(8),
	replicate INT NOT NULL,
	FOREIGN KEY (dataset_ID, dataset_type_ID) REFERENCES datasets (ID, dataset_type_ID));
	
	
CREATE TABLE rna_seq_experiments (
	dataset_ID INT NOT NULL PRIMARY KEY,
	dataset_type_ID INT CHECK (dataset_type_ID = 2),
	environment_ID INT, FOREIGN KEY (environment_ID) REFERENCES environments(ID) ON DELETE CASCADE,
	strain_ID INT, FOREIGN KEY (strain_ID) REFERENCES strains(ID) ON DELETE CASCADE,
	sequencing_type VARCHAR(20),
	machine_ID VARCHAR(20),
	replicate INT NOT NULL,
	FOREIGN KEY (experiment_ID, exp_type_ID) REFERENCES experiments (ID, exp_type_ID));


CREATE TABLE chip_experiments (
	dataset_ID INT NOT NULL PRIMARY KEY,
	dataset_type_ID INT CHECK (dataset_type_ID = 3),
	environment_ID INT, FOREIGN KEY (environment_ID) REFERENCES environments(ID) ON DELETE CASCADE,
	strain_ID INT, FOREIGN KEY (strain_ID) REFERENCES strains(ID) ON DELETE CASCADE,
	antibody VARCHAR(20),
	protocol_type VARCHAR(20),
	target VARCHAR(20), --should reference protein table
	replicate INT NOT NULL,
	FOREIGN KEY (experiment_ID, exp_type_ID) REFERENCES experiments (ID, exp_type_ID));


CREATE TABLE metabolomics_experiments (
	dataset_ID INT NOT NULL PRIMARY KEY,
	dataset_type_ID INT CHECK (dataset_type_ID = 4),
	environment_ID INT, FOREIGN KEY (environment_ID) REFERENCES environments(ID) ON DELETE CASCADE,
	strain_ID INT, FOREIGN KEY (strain_ID) REFERENCES strains(ID) ON DELETE CASCADE,
	name VARCHAR(100),
	mass_spec_ID VARCHAR(20),
	extraction_method VARCHAR(20),
	replicate INT NOT NULL,
	FOREIGN KEY (experiment_ID, exp_type_ID) REFERENCES experiments (ID, exp_type_ID));


CREATE TABLE tss_experiment (
	dataset_ID INT NOT NULL PRIMARY KEY,
	dataset_type_ID INT CHECK (dataset_type_ID = 5),
	environment_ID INT, FOREIGN KEY (environment_ID) REFERENCES environments(ID) ON DELETE CASCADE,
	strain_ID INT, FOREIGN KEY (strain_ID) REFERENCES strains(ID) ON DELETE CASCADE,
	name VARCHAR(100),
	protocol_type VARCHAR(20),
	machine_ID VARCHAR(20),
	replicate INT NOT NULL,
	FOREIGN KEY (experiment_ID, exp_type_ID) REFERENCES experiments (ID, exp_type_ID));
	

CREATE TABLE gcrma_normalization (
	dataset_ID INT NOT NULL PRIMARY KEY,
	dataset_type_ID INT CHECK (dataset_type_ID = 6),
	type VARCHAR(20),
	value DECIMAL,
	leftpos INT,
	rightpos INT,
	algorithm VARCHAR(10));


CREATE TABLE fpkm (
	dataset_ID INT NOT NULL PRIMARY KEY,
	dataset_type_ID INT CHECK (dataset_type_ID = 7),
	type VARCHAR(20),
	value DECIMAL,
	leftpos INT,
	rightpos INT,
	percent_cutoff DECIMAL);


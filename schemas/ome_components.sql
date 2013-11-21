SET search_path TO ecoli;

-- components

DROP TABLE IF EXISTS component_types CASCADE;
DROP TABLE IF EXISTS component CASCADE;
DROP TABLE IF EXISTS genome_region CASCADE;
DROP TABLE IF EXISTS genes CASCADE;
DROP TABLE IF EXISTS proteins CASCADE;
DROP TABLE IF EXISTS gene_ID_protein_ID CASCADE;
DROP TABLE IF EXISTS motifs CASCADE;
DROP TABLE IF EXISTS binding_sites CASCADE;
DROP TABLE IF EXISTS binding_site_motif CASCADE;
DROP TABLE IF EXISTS TU CASCADE;
DROP TABLE IF EXISTS TU_components CASCADE;
DROP TABLE IF EXISTS tss CASCADE;
DROP TABLE IF EXISTS complex CASCADE;
DROP TABLE IF EXISTS complex_components CASCADE;
DROP TABLE IF EXISTS operons CASCADE;
DROP TABLE IF EXISTS operon_components CASCADE;
DROP TABLE IF EXISTS pathways CASCADE;
DROP TABLE IF EXISTS regulatory_network CASCADE;
DROP TABLE IF EXISTS genome CASCADE;



CREATE TABLE genome_region (
    ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
    leftpos INT NOT NULL,
    rightpos INT NOT NULL,
    strand varchar(1), CHECK (strand = '+' OR strand = '-' OR strand = ''),
    UNIQUE(leftpos, rightpos, strand));
    

CREATE TABLE component_types (
	ID INT NOT NULL PRIMARY KEY,
	name VARCHAR(100));
	
	
insert into component_types (ID, name) values
  (1, 'DNA'),
  (2, 'RNA'),
  (3, 'Protein'),
  (4, 'Metabolite');


CREATE TABLE component (
	ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
	component_type_ID INT, FOREIGN KEY (component_type_ID) REFERENCES component_types(ID) ON DELETE CASCADE,
	name VARCHAR(100));


CREATE TABLE dna_types (
    ID INT NOT NULL PRIMARY KEY,
    name VARCHAR(100));
    

CREATE TABLE DNA (
    component_ID INT PRIMARY KEY REFERENCES component(ID),
    component_type_ID INT CHECK (component_type_ID = 1),
    dna_type_ID INT, FOREIGN KEY (dna_type_ID) REFERENCES dna_types(ID) ON DELETE CASCADE,
    genome_region_ID INT, FOREIGN KEY (genome_region_ID) REFERENCES genome_region(ID) ON DELETE CASCADE);   
    
    
insert into dna_types (ID, name) values
    (1, 'binding_site'),
    (2, 'gene');    


CREATE TABLE DNA_binding_site (
    dna_ID INT, FOREIGN KEY (dna_ID) REFERENCES DNA(component_ID) ON DELETE CASCADE,
    dna_type_ID INT CHECK (dna_type_ID = 1),
    bound_component_ID INT, FOREIGN KEY (bound_component_ID) REFERENCES component(ID) ON DELETE CASCADE,
    PRIMARY KEY (dna_ID, bound_component_ID));
        
        
CREATE TABLE gene (
    dna_ID INT PRIMARY KEY REFERENCES DNA(component_ID) ON DELETE CASCADE,
    dna_type_ID INT CHECK (dna_type_ID = 2),
    locus_ID VARCHAR(10),
    name VARCHAR(10),
    long_name VARCHAR(100),
    info text);
    

CREATE TABLE binding_site_motif (
	binding_site_ID INT,
	motif_ID INT,
	PRIMARY KEY (binding_site_ID, motif_ID));
    
    
CREATE TABLE rna_types (
    ID INT NOT NULL PRIMARY KEY,
    name VARCHAR(100));
    
    
CREATE TABLE RNA (
    component_ID INT PRIMARY KEY REFERENCES component(ID),
    component_type_ID INT CHECK (component_type_ID = 2),
    rna_type_ID INT, FOREIGN KEY (rna_type_ID) REFERENCES rna_types(ID) ON DELETE CASCADE,
    genome_region_ID INT, FOREIGN KEY (genome_region_ID) REFERENCES genome_region(ID) ON DELETE CASCADE);
    

CREATE TABLE protein_types (
    ID INT NOT NULL PRIMARY KEY,
    name VARCHAR(100));
 
 
CREATE TABLE protein (
    component_ID INT PRIMARY KEY REFERENCES component(ID),
    component_type_ID INT CHECK (component_type_ID = 3),
    gene_ID INT, FOREIGN KEY (gene_ID) REFERENCES gene(dna_ID) ON DELETE CASCADE,	
    name text,
	long_name text,
	sequence text);  
	
	
CREATE TABLE metabolite_types (
    ID INT NOT NULL PRIMARY KEY,
    name VARCHAR(100));
    
    	
CREATE TABLE metabolite (
    component_ID INT PRIMARY KEY REFERENCES component(ID),
    component_type_ID INT CHECK (component_type_ID = 4),
    metabolite_type_ID INT, FOREIGN KEY (metabolite_type_ID) REFERENCES metabolite_types(ID) ON DELETE CASCADE,
    name VARCHAR(100));  
    
    
insert into rna_types (ID, name) values
    (1, '5prime_triphosphate'),
    (2, 'rRNA'),
    (3, 'tRNA'),
    (4, 'stable');    


CREATE TABLE mRNA (
    rna_ID INT PRIMARY KEY REFERENCES RNA(component_ID) ON DELETE CASCADE,
    rna_type_ID INT CHECK (rna_type_ID = 1),
    gene_ID INT, FOREIGN KEY (gene_ID) REFERENCES gene(dna_ID));
    
    
CREATE TABLE rRNA (
    rna_ID INT PRIMARY KEY REFERENCES RNA(component_ID) ON DELETE CASCADE,
    rna_type_ID INT CHECK (rna_type_ID = 2),
    gene_ID INT, FOREIGN KEY (gene_ID) REFERENCES gene(dna_ID));

    
CREATE TABLE tRNA (
    rna_ID INT PRIMARY KEY REFERENCES RNA(component_ID) ON DELETE CASCADE,
    rna_type_ID INT CHECK (rna_type_ID = 3));
    
    
CREATE TABLE stable_RNA (
    rna_ID INT PRIMARY KEY REFERENCES RNA(component_ID) ON DELETE CASCADE,
    rna_type_ID INT CHECK (rna_type_ID = 4));

    
CREATE TABLE TU (
	ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
	name VARCHAR(100),
	strand CHAR(1));


 CREATE TABLE TU_component (
	tu_ID INT, FOREIGN KEY (tu_ID) REFERENCES TU(ID) ON DELETE CASCADE,
	other_ID INT,
	type VARCHAR(15),
	PRIMARY KEY (tu_ID, other_ID));


CREATE TABLE tss (
	ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
	name VARCHAR(15),
	position INT,
	strand CHAR(1));


 CREATE TABLE complex (
	ID INT NOT NULL PRIMARY KEY DEFAULT nextval('ids'),
	name text,
	type text);


CREATE TABLE complex_components (
	complex_ID INT, FOREIGN KEY (complex_ID) REFERENCES complex(ID) ON DELETE CASCADE,
	other_ID INT,
	type VARCHAR(25),
	coefficient INT,
	PRIMARY KEY (complex_ID, other_ID));   


CREATE TABLE gene_protein (
	gene_ID INT, FOREIGN KEY (gene_ID) REFERENCES gene(dna_ID) ON DELETE CASCADE,
	protein_ID INT, FOREIGN KEY (protein_ID) REFERENCES protein(component_ID) ON DELETE CASCADE,
	PRIMARY KEY(gene_ID, protein_ID));
	
	
CREATE TABLE genome (
	position INT NOT NULL PRIMARY KEY,
	base CHAR(1),
	strand CHAR(1));

import sys

import os
import math
import re
import cobra
#from cobra.io.sbml import create_cobra_model_from_sbml_file
# import copy
# from copy import deepcopy
import cPickle as pickle

from sqlalchemy import text
from sqlalchemy import or_

from PrototypeDB.lib import settings, timing


####Auxiliary functions####

def get_nc_RNA():
    with open(settings.data_directory + "/annotation/ncbi/NC_000913.rnt") as infile:
        nc_rna_bnums = {}
        for i in range(3):
            infile.readline()
        for line in infile:
            nc_rna_bnums[line.split('\t')[5]] = None
    return nc_rna_bnums


def getAttributes(file_name):
    file = open(settings.data_directory + '/annotation/15.1/data/'+file_name,'r')
    atts = []
    attflag = 0
    for line in file.readlines():
        if line[0] != '#': return atts
        if line[0:13] == '# Attributes:': 
            attflag = 1
            continue

        if attflag == 1 and line != '#\n': atts.append(line.lstrip('#    ').rstrip('\n'))


def parseEco_dat(file_name):
    file = open(settings.data_directory + '/annotation/17.1/data/'+file_name,'r')
    master_dict = {}
    attributes = {}
    unique_id = ''

    for line in file.readlines():
        if line[0] == '#': continue

        vals = line.rstrip('\n').split(" - ")

        if vals[0] == 'UNIQUE-ID': 
            master_dict[unique_id] = attributes
            unique_id = vals[1]
            attributes = {}    
        try: 
            tmp = attributes[vals[0]] 
            tmp.append(vals[1])
            attributes[vals[0]] = tmp
        except: 
            try: attributes[vals[0]] = [vals[1]]
            except: None

    return master_dict


def parseEco_col(file_name):
    file = open(settings.data_directory + '/annotation/17.1/data/'+file_name,'r')

    master_dict = {}
    attributes = {}
    unique_id = ''
    header = []

    for line in file.readlines():
        if line[0] == '#': continue

        if line[0:9] == 'UNIQUE-ID': 
            header = line.rstrip('\n').split('\t')
            continue

        vals = line.rstrip('\n').split('\t')

        for i in range(len(vals)):
            if vals[i] == '': continue

            try: 
                tmp = attributes[header[i]] 
                tmp.append(vals[i])
                attributes[header[i]] = tmp
            except: 
                try: attributes[header[i]] = [vals[i]]
                except: None

        master_dict[vals[0]] = attributes
        attributes = {}

    return master_dict


def parseEco_fsa(file_name):
    file = open(settings.data_directory + '/annotation/15.1/data/'+file_name,'r')
    id,seq = '',''
    seq_dict = {}
    for line in file.readlines():
        if line[0] == '>': 
            seq_dict[id] = seq
            seq = ''
            id = line.split()[0].lstrip('>')
        else: seq+=line.rstrip('\n') 
    return seq_dict


def make_citations(ome, ecocyc_entry, id):
    if "CITATIONS" not in ecocyc_entry:
        return
    citations = ecocyc_entry['CITATIONS']
    citation_pmid_list = []
    for citation_entry in citations:
        citation_pmid = citation_entry.strip().strip("[").strip("]").split(":")[0]
        try:
            citation_pmid_list.append(int(citation_pmid))
        except:
            pass
    for citation_pmid in set(citation_pmid_list):
        citation = get_or_create(session, Citation, pmid=citation_pmid)
        ome.execute("INSERT INTO citation_list (table_id, citation_id) " + \
                        "VALUES (%d, %d);" % (id, citation.id))


def checkEqual(list): return list[1:] == list[:-1]


def get_start_strand_bnum(ome, components):
    strand = []
    leftpos = []
    rightpos = []
    bnums = []

    for c in components:
        id_and_type = find_id_and_type(ome, c)
        if id_and_type == '': continue
        if id_and_type[1] == 'gene':
            tmp = ome.execute("SELECT leftpos, rightpos, strand, bnum FROM genes WHERE id = "+str(id_and_type[0])+";")
            result = tmp.fetchone()

            leftpos.append(result[0])
            rightpos.append(result[1])
            strand.append(result[2])
            bnums.append(result[3])

    if not checkEqual(strand) or not strand: return ''
    leftpos.sort()
    rightpos.sort()
    if strand[0] == '+': return [leftpos[0], strand[0], bnums[0]]
    else: return [rightpos[-1], strand[0], bnums[-1]]


def find_id_and_type(ome, ecoid):
    result = ome.execute("SELECT id, type FROM id2otherid WHERE otherID = '%s';"%(ecoid))
    try: return result.fetchall()[0]
    except: return ''


def uniquify(seq, idfun=None): 
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result


def find_id(ome, id):
    result = ome.execute("SELECT id FROM id2otherid WHERE otherID = '"+id+"';")
    try: return result.fetchall()[0][0]
    except: return ''


def make_id(ome, genes):
    global id

    ids = []
    for g in genes: ids.append(find_id(g))

    ome.execute("INSERT INTO complex(id, name) VALUES (%i, '%s');"%(id, ''))

    ome.execute("INSERT INTO id2otherID(id, otherID, type, dataset_id) VALUES " + \
                        "(%i, '%s', '%s', %i);"%(id, 0, 'complex', -1))


    ids = uniquify(ids)
    for w in ids:
        if w == '': continue

        ome.execute("INSERT INTO complex_components(complex_id, other_id, coefficient) VALUES " + \
                        "(%i, %i, %i);"%(id, int(w), 1))

    id += 1
    return id-1


def parse_gpr(gpr):
    gprs = []
    elems = gpr.split(' or ')
    for e in elems: gprs.append(e.replace(')','').replace('(','').split(' and '))
    return gprs    


@timing
def load_genome(genome_filepath=None):
    """use the psql \\copy command to load the genome"""
    if genome_filepath is None:
        genome_filepath = settings.data_directory + "/data/genome/Escherichia_coli_MG1655_genome.tab"
    # using \copy is 6.5x faster
    genome_filepath = os.path.abspath(genome_filepath).replace("\\", "\\\\")
    with open("tmp_command.sql", "w") as outfile:
        outfile.write("set search_path to %s;" % (settings.schema))
        outfile.write(r"\copy genome from '%s'" % (genome_filepath))
    os.system("%s < tmp_command.sql > psql.log 2>&1" % (settings.psql_full))
    os.remove("tmp_command.sql")


def scrub_ecocyc_entry(entry, args=['UNIQUE-ID','COMMON-NAME','TYPES'],extra_args=[]):
    try:
        return {arg: entry[arg] for arg in args+extra_args}
    except:
        return None


def get_or_create_ecocyc_gene(session, base, components, gene_entry):
    vals = scrub_ecocyc_entry(gene_entry, args=['UNIQUE-ID','ACCESSION-1','LEFT-END-POSITION',\
                                                'RIGHT-END-POSITION','TRANSCRIPTION-DIRECTION',\
                                                'COMMON-NAME'])
    if vals is None: return None
    
    gene = session.query(components.Gene).filter(or_(components.Gene.locus_id == vals['ACCESSION-1'][0],\
                                                     components.Gene.name == vals['COMMON-NAME'][0])).first()
    if gene is not None:
        return gene
    else:
        return None
        #session.get_or_create(components.Gene, name=vals['UNIQUE-ID'][0], long_name=vals['COMMON-NAME'][0],\
        #                        leftpos=vals['LEFT-END-POSITION'], rightpos=vals['RIGHT-END-POSITION'],\
        #                        strand=vals['TRANSCRIPTION-DIRECTION'])
    """
    gene = session.get_or_create(components.Gene, )
    if name not in synonyms and name != '':
        synonyms.append(name)
    if ID not in synonyms:
        synonyms.append(ID)
    for synonym in synonyms:
        if session.query(base.id2otherid).filter_by(id=gene.id, other_id=synonym).first() is not None:
                continue  # already exists
        id_entry = session.get_or_create(base.id2otherid,id=gene.id, other_id=synonym,\
                                              type='gene', data_source_id=ecocyc_ID)
    """
    return session.get_or_create(components.Gene, name=vals['UNIQUE-ID'][0], long_name=vals['COMMON-NAME'][0])
            
    
def get_or_create_ecocyc_protein(session, base, components, protein_entry):
    vals = scrub_ecocyc_entry(protein_entry,extra_args=['GENE'])
    if vals is None: return None
    try: gene = get_or_create_ecocyc_gene(session, base, components, ecocyc_genes[vals['GENE'][0]])
    except: return None
    if gene is None: return None
    
    return session.get_or_create(components.Protein, name=vals['UNIQUE-ID'][0], long_name=vals['COMMON-NAME'][0], gene_id=gene.id)


def get_or_create_ecocyc_ligand(session, base, components, ligand_entry):
    vals = scrub_ecocyc_entry(ligand_entry,extra_args=['SMILES'])
    if vals is None: return None
    
    name = vals['COMMON-NAME'][0].replace('\'','[prime]')
    markup = re.compile("<.+?>")
    name = markup.sub("", name)
    
    return session.get_or_create(components.SmallMolecule, name=vals['UNIQUE-ID'][0],\
                                 long_name=name, smiles=vals['SMILES'][0])


def get_or_create_ecocyc_protein_complex(session, base, components, protein_complex_entry):
    vals = scrub_ecocyc_entry(protein_complex_entry,extra_args=['COMPONENTS'])
    if vals is None: return None
    
    protein_complex = session.get_or_create(components.Complex, name=vals['UNIQUE-ID'][0], long_name=vals['COMMON-NAME'][0])
    
    for component in vals['COMPONENTS']:
        try: component_vals = scrub_ecocyc_entry(ecocyc_proteins[component])
        except: 
            try: component_vals = scrub_ecocyc_entry(ecocyc_ligands[component])
            except: continue
        if component_vals is None: continue
        
        if 'Protein-Complexes' in component_vals['TYPES']: 
            complex_component = get_or_create_ecocyc_protein_complex(ecocyc_proteins[component])
        
        elif 'Polypeptides' in component_vals['TYPES']:
            complex_component = get_or_create_ecocyc_protein(ecocyc_proteins[component]) 
            
        else:
            try:
                complex_component = get_or_create_ecocyc_ligand(ecocyc_ligands[component])
            except:
                #TODO: Add in the rest of complex type additions here
                continue
        if complex_component is None: continue
        
        session.get_or_create(components.ComplexComposition, complex_id=protein_complex.id,\
                                                                 component_id=complex_component.id,\
                                                                 stoichiometry=1.)
    return protein_complex    
    
    
def load_motifs(base, components):
    motifs = open(settings.data_directory + '/annotation/ec_annotation_from_Ecocyc_2010July19_wpseudo.gff','r')
 
    
@timing
def load_genes(base,components):
    session = base.Session()
    ##First load annotation file containing merger of ecocyc and NCBI
    ecocyc_ID = session.get_or_create(base.DataSource, name="ecocyc").id
    genes = open(settings.data_directory + '/annotation/ec_annotation_from_Ecocyc_2010July19_wpseudo.gff','r')
    noncoding_genes = get_nc_RNA()
    for line in genes:
        if line[0] == '#':
            continue
        vals = line.rstrip('\r\n').split('\t')
        info = vals[8].replace('\'','_prime_').lstrip('"').rstrip('"').split(';')
        bnum = info[0]
        if bnum == '':
            continue

        gene = session.get_or_create(components.Gene, name=info[1], leftpos=int(vals[3]), rightpos=int(vals[4]), strand=vals[6], locus_id=bnum)
        gene.noncoding = bnum in noncoding_genes
    
        id_entry = session.get_or_create(base.id2otherid,id=gene.id, other_id=bnum,\
                                              type='gene', data_source_id=ecocyc_ID)

    """
    for unique_id,entry in ecocyc_genes.iteritems():    
        
        vals = scrub_ecocyc_entry(entry)
        if vals is None: continue
        
        if 'Protein-Complexes' in vals['TYPES']:
            get_or_create_ecocyc_protein_complex(entry)
            
        elif 'Polypeptides' in vals['TYPES']:
            get_or_create_ecocyc_protein(entry)         
    """

ecocyc_genes = parseEco_dat('genes.dat')        
@timing
def load_proteins(base, components):
    session = base.Session()
    ##First load annotation file containing merger of ecocyc and NCBI
    ecocyc_ID = session.get_or_create(base.DataSource, name="ecocyc").id
    ecocyc_proteins = parseEco_dat('proteins.dat')
    ecocyc_ligands = parseEco_dat('compounds.dat')
    ecocyc_protein_cplxs = parseEco_dat('protligandcplxes.dat')
    protein_sequences = parseEco_fsa('protseq.fsa')
            

    
    for unique_id,entry in ecocyc_proteins.iteritems():    
        
        vals = scrub_ecocyc_entry(entry)
        if vals is None: continue
        
        if 'Protein-Complexes' in vals['TYPES']:
            get_or_create_ecocyc_protein_complex(session, base, components, entry)
            
        elif 'Polypeptides' in vals['TYPES']:
            get_or_create_ecocyc_protein(session, base, components, entry)    
            
                        
    for unique_id,entry in ecocyc_protein_cplxs.iteritems():    
        
        vals = scrub_ecocyc_entry(entry)
        if vals is None: continue
        
        if 'Protein-Complexes' in vals['TYPES'] or 'Protein-Small-Molecule-Complexes' in vals['TYPES']:
            get_or_create_ecocyc_protein_complex(session, base, components, entry)
            
        elif 'Polypeptides' in vals['TYPES']:
            get_or_create_ecocyc_protein(session, base, components, entry)    
            

@timing
def load_chemicals(ome, ecocyc_ID=None):
    if ecocyc_ID is None:
        ecocyc_ID = session.get_or_create(Dataset, name="ecocyc").id
    chemicals = parseEco_dat('compounds.dat')
    markup = re.compile("<.+?>")
    for c in chemicals:
        try:
            ID = chemicals[c]['UNIQUE-ID'][0]
            name = chemicals[c]['COMMON-NAME'][0].replace('\'','[prime]')
            smiles = chemicals[c]['SMILES'][0]
        except:
            continue
        name = markup.sub("", name)
        result = ome.execute("INSERT INTO chemicals(name, smiles) VALUES " + \
                        "('%s', '%s') RETURNING id;" % (name, smiles))
        id = result.fetchone()[0]

        ome.execute("INSERT INTO id2otherID(id, otherID, type, dataset_id) VALUES " + \
                        "(%i, '%s', '%s', %i);"%(id, ID, 'chemical', ecocyc_ID))
        make_citations(ome, chemicals[c], id)
        ome.commit()


@timing
def load_promoters(ome, ecocyc_ID=None):
    if ecocyc_ID is None:
        ecocyc_ID = session.get_or_create(Dataset, name="ecocyc").id
    promoters = parseEco_dat('promoters.dat')

    for p in promoters:
        try:
            ID = promoters[p]['UNIQUE-ID'][0]
            position = int(promoters[p]['ABSOLUTE-PLUS-1-POS'][0])
        except: continue
        try: name = promoters[p]['COMMON-NAME'][0]
        except: name = ''

        result = ome.execute("INSERT INTO tss(name, position) VALUES ('%s', %i) RETURNING id;"%(name, position))
        id = result.fetchone()[0]

        ome.execute("INSERT INTO id2otherID(id, otherID, type, dataset_id) VALUES " + \
                        "(%i, '%s', '%s', %i);"%(id, ID, 'tss', ecocyc_ID))
        make_citations(ome, promoters[p], id)
    ome.commit()


@timing
def load_transcription_units(base, components):
    session = base.Session()
    ##First load annotation file containing merger of ecocyc and NCBI
    ecocyc_ID = session.get_or_create(base.DataSource, name="ecocyc").id
    
    transcription_units = parseEco_dat('transunits.dat')

    for t in transcription_units:
        try: 
            ID = transcription_units[t]['UNIQUE-ID'][0]
            components = transcription_units[t]['COMPONENTS']
        except: continue
        try: name = transcription_units[t]['COMMON-NAME'][0]
        except: name = ''

        start_strand_bnum = get_start_strand_bnum(ome, components)
        if start_strand_bnum == '': continue

        result = ome.execute("INSERT INTO TU(name, strand) VALUES ('%s', '%s') RETURNING id;"%(name, start_strand_bnum[1]))
        tu_id = result.fetchone()[0]
        ome.execute("INSERT INTO id2otherID(id, otherID, type, dataset_id) VALUES " + \
                    "(%i, '%s', '%s', %i);"%(tu_id, ID, 'TU', ecocyc_ID))

        has_tss_flag = 0
        genes = []
        for c in components:
            id_and_type = find_id_and_type(ome, c)

            if id_and_type == '': continue
            if id_and_type[1] == 'binding_site': continue
            if id_and_type[1] == 'gene': genes.append(id_and_type[0])
            if id_and_type[1] == 'tss': has_tss_flag = 1

            ome.execute("INSERT INTO TU_components(tu_id, other_id, type) VALUES (%i, %i, '%s');"%(tu_id, id_and_type[0], id_and_type[1]))

        if not has_tss_flag:
            result = ome.execute("INSERT INTO tss(name, position) VALUES ('%s', %i) RETURNING id;"%('start_'+start_strand_bnum[2], start_strand_bnum[0]))
            tss_id = result.fetchone()[0]
            ome.execute("INSERT INTO id2otherid(id, otherID, type, dataset_id) VALUES " + \
                        "(%i, '%s', '%s', %i);"%(tss_id, start_strand_bnum[2], 'gene_start_as_tss', ecocyc_ID))
            ome.execute("INSERT INTO TU_components(tu_id, other_id, type) VALUES (%i, %i, '%s');"%(tu_id, tss_id, 'tss'))
        make_citations(ome, transunits[t], tu_id)


@timing
def load_bindsites(base, components):
    session = base.Session()
    ##First load annotation file containing merger of ecocyc and NCBI
    ecocyc_ID = session.get_or_create(base.DataSource, name="ecocyc").id
    ecocyc_binding_sites = parseEco_dat('dnabindsites.dat')
    ecocyc_regulation = parseEco_dat('regulation.dat')
    
    
    for unique_id,entry in ecocyc_binding_sites.iteritems():
        
        vals = scrub_ecocyc_entry(entry, args=['UNIQUE-ID','TYPES','ABS-CENTER-POS'])
        if vals is None: continue
        
        if 'DNA-Binding-Sites' in vals['TYPES']:
            
            try: centerpos = math.floor(float(vals['ABS-CENTER-POS'][0]))
            except: continue
            
            try:
                length = int(ecocyc_binding_sites[unique_id]['SITE-LENGTH'][0])
                leftpos = math.floor(centerpos-(length/2))
                rightpos = math.floor(centerpos+(length/2))
            except: 
                length = 0
                leftpos = centerpos
                rightpos = centerpos
            
            session.get_or_create(components.DnaBindingSite, name=vals['UNIQUE-ID'][0], leftpos=leftpos,\
                                  rightpos=rightpos, strand='+', centerpos=centerpos, width=length)                  

    
    for unique_id,entry in ecocyc_regulation.iteritems():
        vals = scrub_ecocyc_entry(entry, args=['UNIQUE-ID','TYPES','ASSOCIATED-BINDING-SITE','REGULATOR'])
        if vals is None: continue
        
        if 'Transcription-Factor-Binding' in vals['TYPES']:
            regulator = session.get_or_create(components.Component, name=vals['REGULATOR'][0])
            
            binding_site = session.query(components.DnaBindingSite).filter_by(name=vals['ASSOCIATED-BINDING-SITE'][0]).first()
            if binding_site is None: continue
        
            tf_binding_complex = session.get_or_create(components.Complex, name=vals['UNIQUE-ID'][0])
            
            session.get_or_create(components.ComplexComposition, complex_id=tf_binding_complex.id,\
                                                                     component_id=regulator.id,\
                                                                     stoichiometry=1.)
            
            session.get_or_create(components.ComplexComposition, complex_id=tf_binding_complex.id,\
                                                                     component_id=binding_site.id,\
                                                                     stoichiometry=1.)

                              
@timing
def load_kegg_pathways(session):
    # this table should probably get re-designed properly with a pathways table
    # and a separate pathway genes secondary table - TODO
    kegg_pathways = open(settings.data_directory + '/annotation/KEGG/ecoli_KEGG_pathways.txt','r')
    kegg_dataset = get_or_create(session, Dataset, name="kegg")
    for line in kegg_pathways.readlines():
        vals = line.rstrip('\r\n').split('\t')
        bnums = vals[2].split(',')
        pathway_ID = vals[0]
        pathway_name = vals[1].strip()
        kegg_pathway = Pathway(name=pathway_name)
        for bnum in bnums:
            kegg_pathway.genes.append(
                session.query(Gene).filter_by(bnum=bnum).first())
        session.add(kegg_pathway)
        session.commit()
        id_entry = id2otherid(dataset=kegg_dataset, otherid=pathway_ID)
        id_entry.id = kegg_pathway.id
        session.add(id_entry)
        session.commit()


@timing
def load_regulatory_network(ome):
    sigma_network = open('../../reconstruction/regulonDB/network_sigma_gene.txt','r')
    #sigma_dict = {'Sigma19':'', 'Sigma24':'', 'Sigma38':'b2741'
    for line in sigma_network.readlines():
        vals = line.split('\t')
        if vals[0] != 'Sigma38': continue
        try:
            rpoS = 'rpoS'
            rpoS_bnum = 'b2741'
            #print "INSERT INTO regulatory_network(reg_gene, reg_bnum, regd_gene, regd_bnum, direction) " + \
            #            "VALUES ('%s', '%s', '%s', '%s', '%s');"%(rpoS, rpoS_bnum, vals[1], vals[3], vals[2])
            #ome.execute("INSERT INTO regulatory_network(reg_gene, reg_bnum, regd_gene, regd_bnum, direction) VALUES ('%s', '%s', '%s', '%s', '%s');"%(rpoS, rpoS_bnum, vals[1], vals[3], vals[2]))
            ome.commit()
        except: None
    sigma_network.close()

    regulatory_network = open(settings.trn_directory + 'reconstruction/regulonDB/regDB_regulatory_network.txt','r')
    for line in regulatory_network.readlines():
        vals = line.split('\t')
        if len(vals) < 5: continue
        quality = ''
        evidence = ''
        if len(vals) >= 6:
            evidence, quality = parse_evidence(vals[6].strip())
        ome.execute("INSERT INTO regulatory_network(reg_gene, reg_bnum, regd_gene, regd_bnum, direction, quality, evidence) " + \
                    "VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s');"%(vals[1], vals[2], vals[3], vals[4], vals[5], quality, evidence))
        ome.commit()
    regulatory_network.close()


if __name__ == "__main__":


    # can now use omeORM
    from PrototypeDB.orm import base,components
    session = base.Session()
    # find ecocyc
    ecocyc_ID = session.get_or_create(base.DataSource, name="ecocyc").id


    # components can now be loaded into the database
    #load_genome()
    load_genes(session, ecocyc_ID)
    #load_proteins(session, ecocyc_ID)
    #load_chemicals(session, ecocyc_ID)
    #load_promoters(session, ecocyc_ID)
    #load_transunits(session, ecocyc_ID)
    #load_bindsites(session, ecocyc_ID)
    #load_protein_complexes(session, ecocyc_ID)
    #load_protein_ligand_complexes(session, ecocyc_ID)

    #load_iMC1010(session)
    #load_cobra_model(session)
    #load_kegg_pathways(session)
    #load_regulatory_network(session)

    session.commit()
    session.close()

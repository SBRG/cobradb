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
    file = open(settings.data_directory + '/annotation/15.1/data/'+file_name,'r')
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
    file = open(settings.data_directory + '/annotation/15.1/data/'+file_name,'r')

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
        if session.query(components.Gene).filter(components.Gene.locus_id == bnum).first() is not None:
            continue  # already exists
        gene = session.get_or_create(components.Gene, name=info[1], leftpos=int(vals[3]), rightpos=int(vals[4]), strand=vals[6], locus_id=bnum)
        gene.noncoding = bnum in noncoding_genes
    
        id_entry = session.get_or_create(base.id2otherid,id=gene.id, other_id=bnum,\
                                              type='gene', data_source_id=ecocyc_ID)

    genes = parseEco_dat('genes.dat')

    for g in genes: 
        try:
            ID = genes[g]['UNIQUE-ID'][0]
            bnum = genes[g]['ACCESSION-1'][0]
            start = int(genes[g]['LEFT-END-POSITION'][0])
            stop = int(genes[g]['RIGHT-END-POSITION'][0])
            strand = genes[g]['TRANSCRIPTION-DIRECTION'][0]
        except:
            continue
        try: name = genes[g]['COMMON-NAME'][0]
        except: name = ''
        try: synonyms = genes[g]['SYNONYMS']
        except: synonyms = []
        gene = session.query(components.Gene).filter(components.Gene.locus_id == bnum).first()
        if gene is None:
            continue  # not found
        if name not in synonyms and name != '':
            synonyms.append(name)
        if ID not in synonyms:
            synonyms.append(ID)
        for synonym in synonyms:
            if session.query(base.id2otherid).filter_by(
                        id=gene.id, other_id=synonym).first() is not None:
                continue  # already exists
            id_entry = session.get_or_create(base.id2otherid,id=gene.id, other_id=synonym,\
                                              type='gene', data_source_id=ecocyc_ID)

        
@timing
def load_proteins(session, ecocyc_ID=None):
    if ecocyc_ID is None:
        ecocyc_ID = session.get_or_create(Dataset, name="ecocyc").id
    proteins = parseEco_dat('proteins.dat')
    protein_sequences = parseEco_fsa('protseq.fsa')
    for p in proteins:
        try:
            ID = proteins[p]['UNIQUE-ID'][0]
            types = proteins[p]['TYPES']
            name = proteins[p]['COMMON-NAME'][0].replace('\'','[prime]')
            geneID = proteins[p]['GENE'][0]
        except:
            continue
        short_name = name.split()
        # TODO ensure no protein exists already
        new_protein = Protein()
        new_protein.name = short_name[0]
        new_protein.long_name = name
        try: new_protein.sequence = protein_sequences[ID]
        except: new_protein.sequence = ''     
        session.add(new_protein)
        session.commit()
        # add in the related genes
        new_protein.genes.extend(session.query(Gene).filter(
            id2otherid.otherid==geneID, Gene.id == id2otherid.id).all())
        # create id2otherid entries
        synonyms = []
        if "SYNONYMS" in proteins[p]:
            for entry in proteins[p]["SYNONYMS"]:
                if ";" not in entry:
                    synonyms.append(entry)
        if ID not in synonyms:
            synonyms.append(ID)
        for synonym in synonyms:
            id_entry = id2otherid()
            id_entry.id = new_protein.id
            id_entry.otherid = synonym
            id_entry.type = "protein"
            id_entry.dataset_id = ecocyc_ID
            session.add(id_entry)
        make_citations(session, proteins[p], new_protein.id)
        session.commit()


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
    conn = base.engine.connect()
    
    ##First load annotation file containing merger of ecocyc and NCBI
    ecocyc_ID = session.get_or_create(base.DataSource, name="ecocyc").id
    dna_binding_bound_component_insert = data.dna_binding_bound_component_association.insert()

    binding_sites = parseEco_dat('dnabindsites.dat')
    regulation = parseEco_dat('regulation.dat')
    
    ##commented out code is for adding support for binding sites without a defined position
    #ome.execute("DROP TABLE tu_bind_sites CASCADE;")
    #ome.execute("CREATE TABLE tu_bind_sites ( " + \
    #            "tu_id INT, FOREIGN KEY (tu_id) REFERENCES TU(id) ON DELETE CASCADE, " + \
    #            "binding_site_id INT, FOREIGN KEY (binding_site_id) REFERENCES binding_sites(id) ON DELETE CASCADE, " + \
    #            "PRIMARY KEY(tu_id, binding_site_id));")

    for b in binding_sites:
        try: 
            ID = binding_sites[b]['UNIQUE-ID'][0]
            centerpos = float(binding_sites[b]['ABS-CENTER-POS'][0])
            #component_of = bindsites[b]['COMPONENT-OF']
        except: continue
                    
        """
        Going to try and find what is bound at each binding site, could be a
        protein, protein_complex, protein_ligand_complex, or a 
        protein_complex_ligand_complex for now... All of these complexes will
        be represented by the complex class regardless.
        """
        try: 
            for regulation_id in binding_sites[b]['INVOLVED-IN-REGULATION']:
                try:
                    regulator_id = regulation[regulation_id]['REGULATOR']
                    if ome.query(components.Complex).count() == 0: 
                        load_proteins(base, components)
                        continue
                    db_id = ome.query(base.IDSynonyms).filter(base.IDSynonyms.synonym == regulator_id).id
                    conn.execute(dna_binding_bound_component_insert, chip_experiment_id=experiment_id, chip_peak_analysis_id=peak_analysis.id)
                except: None
        except: None
                        
        try:
            length = int(binding_sites[b]['SITE-LENGTH'][0])
            leftpos = math.floor(centerpos-(length/2))
            rightpos = math.floor(centerpos+(length/2))
        except: 
            leftpos = centerpos
            rightpos = centerpos


        strand = ''
        
        binding_site = ome.get_or_create(components.DnaBindingSite, name=id, leftpos=leftpos, rightpos=rightpos,\
                                         strand=strand, bound_component_id=None)


    ome.commit()


def load_protein_complexes(ome, ecocyc_ID=None):    
    if ecocyc_ID is None:
        ecocyc_ID = session.get_or_create(Dataset, name="ecocyc").id
    complexes = parseEco_col('protcplxs.col')

    for c in complexes:
        try: 
            ID = complexes[c]['UNIQUE-ID'][0]
            name = complexes[c]['NAME'][0].replace('\'','[prime]')
            genes = complexes[c]['GENE-NAME']
            gene_ids = complexes[c]['GENE-ID']
            subunits_comp = complexes[c]['SUBUNIT-COMPOSITION']
        except:
            continue

        result = ome.execute("INSERT INTO complex(name, type) VALUES ('%s', '%s') RETURNING id;"%(name, 'protein_complex'))
        id = result.fetchone()[0]

        subunits = subunits_comp[0].split(',')
        for s in subunits:
            vals = s.split('*')

            id_and_type = find_id_and_type(ome, vals[1])
            if id_and_type == '':
                continue

            ome.execute("INSERT INTO complex_components(complex_id, other_id, type, coefficient) VALUES " + \
                            "(%i, %i, '%s', %i);"%(id, id_and_type[0], id_and_type[1], int(vals[0])))

        ome.execute("INSERT INTO id2otherID(id, otherID, type, dataset_id) VALUES " + \
                        "(%i, '%s', '%s', %i);"%(id, ID, 'protein_complex', ecocyc_ID))


@timing
def load_protein_ligand_complexes(ome, ecocyc_ID=None):
    if ecocyc_ID is None:
        ecocyc_ID = session.get_or_create(Dataset, name="ecocyc").id
    proteins = parseEco_dat('protligandcplxes.dat')

    for p in proteins:
        try:
            ID = proteins[p]['UNIQUE-ID'][0]
            types = proteins[p]['TYPES']
            name = proteins[p]['COMMON-NAME'][0].replace('\'','[prime]')
            components = proteins[p]['COMPONENTS']
        except: continue

        if types[0] != 'Protein-Small-Molecule-Complexes': continue

        result = ome.execute("INSERT INTO complex(name, type) VALUES ('%s', '%s') RETURNING id;"%(name, 'protein_ligand_complex'))
        complex_id = result.fetchone()[0]
        ome.execute("INSERT INTO id2otherid(id, otherid, type, dataset_id) VALUES " + \
                    "(%i, '%s', '%s', %i);"%(complex_id, ID, 'protein_ligand_complex', ecocyc_ID))


        for c in components:
            result = ome.execute("SELECT g.id FROM genes g, proteins p, gene_id_protein_id gpw, id2otherid w2 " + \
                        "WHERE gpw.gene_id = g.id AND gpw.protein_id = p.id AND p.id = w2.id " + \
                        "AND w2.otherid = '%s';"%(c))
            try: 
                gene_id = result.fetchone()[0]
                protein_id = find_id_and_type(ome, c)


                ome.execute("INSERT INTO complex_components(complex_id, other_id, type, coefficient) " + \
                            "VALUES (%i, %i, '%s', %i);"%(complex_id, protein_id, 'protein_complex', 1))

                result = ome.execute("SELECT * from gene_id_protein_id WHERE gene_id = %i AND protein_id = %i;"%(gene_id, protein_id[0]))
                try: result.fetchone()[0]
                except:
                    ome.execute("INSERT INTO gene_id_protein_id(gene_id, protein_id) VALUES " + \
                                "(%i, %i);"%(gene_id, protein_id[0]))

            except: 

                id_and_type = find_id_and_type(ome, c)
                if id_and_type == '': 
                    #print c
                    continue
                ome.execute("INSERT INTO complex_components(complex_id, other_id, type, coefficient) " +
                            "VALUES(%i, %i, '%s', %i);"%(complex_id, id_and_type[0], id_and_type[1], 1))
        make_citations(ome, proteins[p], complex_id)
        ome.commit()  
            
                              
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

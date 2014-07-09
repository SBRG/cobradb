import sys, os, math, re

from sqlalchemy import text, or_

from om import settings, timing


def getAttributes(file_name):
    file = open(settings.data_directory + '/annotation/MetaCyc/17.1/data/'+file_name,'r')
    atts = []
    attflag = 0
    for line in file.readlines():
        if line[0] != '#': return atts
        if line[0:13] == '# Attributes:':
            attflag = 1
            continue

        if attflag == 1 and line != '#\n': atts.append(line.lstrip('#    ').rstrip('\n'))


def parse_metacyc_dat(file_name):
    file = open(settings.data_directory + '/annotation/MetaCyc/17.1/data/'+file_name,'r')
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


def parse_metacyc_col(file_name):
    file = open(settings.data_directory + '/annotation/MetaCyc/17.1/data/'+file_name,'r')

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


def parse_metacyc_fsa(file_name):
    file = open(settings.data_directory + '/annotation/MetaCyc/17.1/data/'+file_name,'r')
    id,seq = '',''
    seq_dict = {}
    for line in file.readlines():
        if line[0] == '>':
            seq_dict[id] = seq
            seq = ''
            id = line.split()[0].lstrip('>')
        else: seq+=line.rstrip('\n')
    return seq_dict


def make_citations(ome, metacyc_entry, id):
    if "CITATIONS" not in metacyc_entry:
        return
    citations = metacyc_entry['CITATIONS']
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


def parse_gpr(gpr):
    gprs = []
    elems = gpr.split(' or ')
    for e in elems: gprs.append(e.replace(')','').replace('(','').split(' and '))
    return gprs


def scrub_metacyc_entry(entry, args=['UNIQUE-ID','COMMON-NAME','TYPES'],extra_args=[]):
    try:
        return {arg: entry[arg] for arg in args+extra_args}
    except:
        return None


def get_gene_with_metacyc(session, base, components, gene_entry):
    vals = scrub_metacyc_entry(gene_entry, args=['UNIQUE-ID','ACCESSION-1','LEFT-END-POSITION',\
                                                          'RIGHT-END-POSITION','TRANSCRIPTION-DIRECTION',\
                                                          'COMMON-NAME'])
    if vals is None: return None

    gene = session.query(components.Gene).filter(or_(components.Gene.locus_id == vals['ACCESSION-1'][0],\
                                                     components.Gene.name == vals['COMMON-NAME'][0])).first()
    if gene is None:
        print vals['ACCESSION-1'][0]
        print 'Exception, MetaCyc gene:'+vals['ACCESSION-1'][0]+' not found in genbank'
        return None
    return gene



def get_protein_with_metacyc(session, base, components, protein_entry):
    vals = scrub_metacyc_entry(protein_entry,extra_args=['GENE'])
    if vals is None: return None

    """First see if the gene exists from genbank, if not log an error message and return"""
    try:
        gene_entry = metacyc_genes[vals['GENE'][0]]
        gene_vals = scrub_metacyc_entry(gene_entry, args=['UNIQUE-ID','ACCESSION-1','LEFT-END-POSITION',\
                                                          'RIGHT-END-POSITION','TRANSCRIPTION-DIRECTION',\
                                                          'COMMON-NAME'])


        gene = session.query(components.Gene).filter(or_(components.Gene.locus_id == vals['ACCESSION-1'][0],\
                                                         components.Gene.name == vals['COMMON-NAME'][0])).first()
        if gene is None:
            print 'Exception, MetaCyc gene:'+vals['ACCESSION-1'][0]+' not found in genbank'
            return None

    except: return None

    """If the gene exists in genbank, then the protein should exist in genbank as well,
       if not log an error message and return
    """
    try:
        return session.query(components.Protein).filter(components.Protein.gene_id == gene.id).one()
    except:
        return None

def update_protein_with_metacyc(session, base, components, protein_entry):
    vals = scrub_metacyc_entry(protein_entry,extra_args=['GENE'])
    if vals is None: return None

    protein = get_protein_with_metacyc(session, base, components, protein_entry)
    if not protein: return None
    protein.long_name = vals['COMMON-NAME']

    session.merge(protein)
    session.commit()


def get_or_create_metacyc_ligand(session, base, components, ligand_entry):
    vals = scrub_metacyc_entry(ligand_entry,extra_args=['SMILES'])
    if vals is None: return None

    name = vals['COMMON-NAME'][0].replace('\'','[prime]')
    markup = re.compile("<.+?>")
    name = markup.sub("", name)

    return session.get_or_create(components.SmallMolecule, name=vals['UNIQUE-ID'][0],\
                                 long_name=name, smiles=vals['SMILES'][0])


def get_or_create_metacyc_protein_complex(session, base, components, protein_complex_entry):
    vals = scrub_metacyc_entry(protein_complex_entry,extra_args=['COMPONENTS'])
    if vals is None: return None

    protein_complex = session.get_or_create(components.Complex, name=vals['UNIQUE-ID'][0], long_name=vals['COMMON-NAME'][0])

    for component in vals['COMPONENTS']:
        try: component_vals = scrub_metacyc_entry(metacyc_proteins[component])
        except:
            try: component_vals = scrub_metacyc_entry(metacyc_ligands[component])
            except: continue

        if component_vals is None: continue

        if 'Protein-Complexes' in component_vals['TYPES']:
            complex_component = get_or_create_metacyc_protein_complex(session, base, components, metacyc_proteins[component])

        elif 'Polypeptides' in component_vals['TYPES']:
            complex_component = get_protein_with_metacyc(session, base, components, metacyc_proteins[component])

        else:
            try:
                complex_component = get_or_create_metacyc_ligand(metacyc_ligands[component])
            except:
                #TODO: Add in the rest of complex type additions here
                continue
        if complex_component is None: continue

        session.get_or_create(components.ComplexComposition, complex_id=protein_complex.id,\
                                                                 component_id=complex_component.id,\
                                                                 stoichiometry=1.)
    return protein_complex


def get_or_create_metacyc_transcription_unit(session, base, components, tu_entry):
    vals = scrub_metacyc_entry(tu_entry,extra_args=['COMPONENTS'])
    if vals is None: return None

    genes = []
    tss = None
    for metacyc_component in vals['COMPONENTS']:
        try:
            gene_entry = metacyc_genes[metacyc_component]
            genes.append(get_gene_with_metacyc(session, base, components, gene_entry))
        except: None
        try:
            promoter_entry = metacyc_promoters[metacyc_component]
            tss_vals = scrub_metacyc_entry(promoter_entry,extra_args=['ABSOLUTE-PLUS-1-POS'])
            tss = tss_vals['ABSOLUTE-PLUS-1-POS'][0]
        except: None
    if len(genes) == 0: return

    try:
        leftpos = min([x.leftpos for x in genes])
        rightpos = max([x.rightpos for x in genes])
        strand = genes[0].strand
    except: return

    try: name = vals['COMMON-NAME'][0]
    except: name = vals['UNIQUE-ID'][0]

    if not tss and strand == '+':
        tss = leftpos
    elif not tss and strand =='-':
        tss = rightpos

    if strand == '+':
        tu = session.get_or_create(components.TU, name=name, leftpos=tss, rightpos=rightpos, strand=strand)
    else:
        tu = session.get_or_create(components.TU, name=name, leftpos=leftpos, rightpos=tss, strand=strand)

    for gene in genes:
        session.get_or_create(components.TUGenes, tu_id=tu.id, gene_id=gene.id)

    return tu


@timing
def load_genbank(genbank_file, base, components):
    """
        Original code by halatif and jtan, modified here to load into om db
    """

    from Bio import SeqIO

    session = base.Session()
    gb_file = SeqIO.read(settings.data_directory+'/annotation/GenBank/'+genbank_file,'gb')

    for feature in gb_file.features:
        om_gene = {'long_name':''}
        om_protein = {'long_name':''}

        if feature.type == 'CDS' and 'product' in feature.qualifiers:  #Only cares for coding sequences which are not pseudogenes

            om_gene['name'] = feature.qualifiers['gene'][0]
            om_gene['leftpos'] = feature.location.start
            om_gene['rightpos'] = feature.location.end

            if feature.strand == 1: om_gene['strand'] = '+'
            elif feature.strand == -1: om_gene['strand'] = '-'

            om_gene['locus_id'] = feature.qualifiers['locus_tag'][0]
            om_gene['info'] = feature.qualifiers['product'][0]

            if 'function' in feature.qualifiers:
                om_gene['info'] = om_gene['info'] + ',' + feature.qualifiers['function'][0]
                om_protein['long_name'] = feature.qualifiers['function'][0]

            gene = components.Gene(**om_gene)
            session.add(gene)
            session.flush()

            om_protein['name'] = feature.qualifiers['protein_id'][0]
            om_protein['gene_id'] = gene.id

            session.add(components.Protein(**om_protein))

    session.commit()
    session.close()



@timing
def load_motifs(base, components):
    motifs = open(settings.data_directory + '/annotation/ec_annotation_from_metacyc_2010July19_wpseudo.gff','r')


@timing
def load_metacyc_proteins(base, components):
    session = base.Session()
    ##First load annotation file containing merger of metacyc and NCBI
    metacyc_ID = session.get_or_create(base.DataSource, name="metacyc").id



    for unique_id,entry in metacyc_proteins.iteritems():

        vals = scrub_metacyc_entry(entry)
        if vals is None: continue

        if 'Protein-Complexes' in vals['TYPES']:
            get_or_create_metacyc_protein_complex(session, base, components, entry)

        elif 'Polypeptides' in vals['TYPES']:
            update_protein_with_metacyc(session, base, components, entry)


    for unique_id,entry in metacyc_protein_cplxs.iteritems():

        vals = scrub_metacyc_entry(entry)
        if vals is None: continue

        if 'Protein-Complexes' in vals['TYPES'] or 'Protein-Small-Molecule-Complexes' in vals['TYPES']:
            get_or_create_metacyc_protein_complex(session, base, components, entry)

        elif 'Polypeptides' in vals['TYPES']:
            get_or_create_metacyc_protein(session, base, components, entry)


@timing
def load_metacyc_transcription_units(base, components):
    session = base.Session()
    ##First load annotation file containing merger of metacyc and NCBI
    metacyc_ID = session.get_or_create(base.DataSource, name="metacyc").id



    for unique_id,entry in metacyc_tus.iteritems():

        vals = scrub_metacyc_entry(entry,extra_args=['COMPONENTS'])
        if vals is None: continue

        if 'Transcription-Units' in vals['TYPES']:
            get_or_create_metacyc_transcription_unit(session, base, components, entry)




        """
        start_strand_bnum = get_start_strand_bnum(ome, components)
        if start_strand_bnum == '': continue

        result = ome.execute("INSERT INTO TU(name, strand) VALUES ('%s', '%s') RETURNING id;"%(name, start_strand_bnum[1]))
        tu_id = result.fetchone()[0]
        ome.execute("INSERT INTO id2otherID(id, otherID, type, dataset_id) VALUES " + \
                    "(%i, '%s', '%s', %i);"%(tu_id, ID, 'TU', metacyc_ID))

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
                        "(%i, '%s', '%s', %i);"%(tss_id, start_strand_bnum[2], 'gene_start_as_tss', metacyc_ID))
            ome.execute("INSERT INTO TU_components(tu_id, other_id, type) VALUES (%i, %i, '%s');"%(tu_id, tss_id, 'tss'))
        make_citations(ome, transunits[t], tu_id)
        """


@timing
def load_metacyc_bindsites(base, components):
    session = base.Session()
    ##First load annotation file containing merger of metacyc and NCBI
    metacyc_ID = session.get_or_create(base.DataSource, name="metacyc").id
    metacyc_binding_sites = parse_metacyc_dat('dnabindsites.dat')
    metacyc_regulation = parse_metacyc_dat('regulation.dat')


    for unique_id,entry in metacyc_binding_sites.iteritems():

        vals = scrub_metacyc_entry(entry, args=['UNIQUE-ID','TYPES','ABS-CENTER-POS'])
        if vals is None: continue

        if 'DNA-Binding-Sites' in vals['TYPES']:

            try: centerpos = math.floor(float(vals['ABS-CENTER-POS'][0]))
            except: continue

            try:
                length = int(metacyc_binding_sites[unique_id]['SITE-LENGTH'][0])
                leftpos = math.floor(centerpos-(length/2))
                rightpos = math.floor(centerpos+(length/2))
            except:
                length = 0
                leftpos = centerpos
                rightpos = centerpos

            session.get_or_create(components.DnaBindingSite, name=vals['UNIQUE-ID'][0], leftpos=leftpos,\
                                  rightpos=rightpos, strand='+', centerpos=centerpos, width=length)


    for unique_id,entry in metacyc_regulation.iteritems():
        vals = scrub_metacyc_entry(entry, args=['UNIQUE-ID','TYPES','ASSOCIATED-BINDING-SITE','REGULATOR'])
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
def load_metacyc_chemicals(ome, metacyc_ID=None):
    if metacyc_ID is None:
        metacyc_ID = session.get_or_create(Dataset, name="metacyc").id
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
                        "(%i, '%s', '%s', %i);"%(id, ID, 'chemical', metacyc_ID))
        make_citations(ome, chemicals[c], id)
        ome.commit()


@timing
def load_metacyc_promoters(ome, metacyc_ID=None):
    if metacyc_ID is None:
        metacyc_ID = session.get_or_create(Dataset, name="metacyc").id
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
                        "(%i, '%s', '%s', %i);"%(id, ID, 'tss', metacyc_ID))
        make_citations(ome, promoters[p], id)
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



metacyc_genes = parse_metacyc_dat('genes.dat')
metacyc_promoters = parse_metacyc_dat('promoters.dat')
metacyc_proteins = parse_metacyc_dat('proteins.dat')
metacyc_ligands = parse_metacyc_dat('compounds.dat')
metacyc_protein_cplxs = parse_metacyc_dat('protligandcplxes.dat')
metacyc_tus = parse_metacyc_dat('transunits.dat')



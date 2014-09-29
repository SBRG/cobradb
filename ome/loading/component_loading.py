import sys, os, math, re

from sqlalchemy import text, or_, and_, func

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


def scrub_metacyc_entry(entry, args=['UNIQUE-ID','TYPES'],extra_args=[]):
    try:
        return {arg: entry[arg] for arg in args+extra_args}
    except:
        return None


def get_gene_with_metacyc(session, base, components, genome, gene_entry):
    vals = scrub_metacyc_entry(gene_entry, args=['UNIQUE-ID','ACCESSION-1','LEFT-END-POSITION',\
                                                          'RIGHT-END-POSITION','TRANSCRIPTION-DIRECTION',\
                                                          'COMMON-NAME'])
    if vals is None: return None

    gene = session.query(components.Gene).filter(and_(components.Gene.genome_id == genome.id,
                                                      or_(components.Gene.locus_id == vals['ACCESSION-1'][0],
                                                          components.Gene.name == vals['COMMON-NAME'][0]))).first()
    if gene is None:
        print 'Exception, MetaCyc gene:'+vals['ACCESSION-1'][0]+' not found in genbank'

    return gene


def update_gene_with_metacyc(session, base, components, gene_entry):
    vals = scrub_metacyc_entry(gene_entry, args=['UNIQUE-ID','ACCESSION-1','LEFT-END-POSITION',\
                                                          'RIGHT-END-POSITION','TRANSCRIPTION-DIRECTION',\
                                                          'COMMON-NAME'])
    if vals is None: return None

    gene = get_gene_with_metacyc(session, base, components, gene_entry)
    if not gene: return None

    gene.long_name = vals['COMMON-NAME'][0]

    session.add(gene)
    session.flush()


def get_protein_with_metacyc(session, base, components, genome, protein_entry):
    vals = scrub_metacyc_entry(protein_entry,extra_args=['GENE'])
    if vals is None: return None

    """First see if the gene exists from genbank, if not log an error message and return"""

    try: gene_entry = metacyc_genes[vals['GENE'][0]]
    except:
        print 'MetaCyc issue, gene: '+vals['GENE'][0]+' does not exist'
        return None

    gene = get_gene_with_metacyc(session, base, components, genome, gene_entry)


    """If the gene exists in genbank, then the protein should exist in genbank as well,
       if not log an error message and return
    """

    try:
        return session.query(components.Protein).filter(components.Protein.gene_id == gene.id).one()
    except:
        return None

def update_protein_with_metacyc(session, base, components, genome, protein_entry):
    vals = scrub_metacyc_entry(protein_entry,extra_args=['GENE', 'COMMON-NAME'])
    if vals is None: return None

    protein = get_protein_with_metacyc(session, base, components, genome, protein_entry)
    if not protein: return None

    protein.long_name = vals['COMMON-NAME'][0]

    session.add(protein)
    session.flush()



def get_or_create_metacyc_ligand(session, base, components, ligand_entry):
    vals = scrub_metacyc_entry(ligand_entry,extra_args=['COMMON-NAME','SMILES'])
    if vals is None: return None

    name = vals['COMMON-NAME'][0].replace('\'','[prime]')
    markup = re.compile("<.+?>")
    name = markup.sub("", name)

    return session.get_or_create(components.Metabolite, name=vals['UNIQUE-ID'][0],\
                                 long_name=name, smiles=vals['SMILES'][0])


def get_or_create_metacyc_protein_complex(session, base, components, genome, protein_complex_entry):
    #if protein_complex_entry['UNIQUE-ID'][0] == 'CPLX0-226': print protein_complex_entry
    #print protein_complex_entry['UNIQUE-ID'][0]
    vals = scrub_metacyc_entry(protein_complex_entry,extra_args=['COMMON-NAME','COMPONENTS'])
    if vals is None: return None

    protein_complex = session.get_or_create(components.Complex, name=vals['UNIQUE-ID'][0], long_name=vals['COMMON-NAME'][0])

    for component in vals['COMPONENTS']:

        component_vals = None
        complex_component = None

        if component in metacyc_proteins:
            component_vals = scrub_metacyc_entry(metacyc_proteins[component])

        elif component in metacyc_ligands:
            component_vals = scrub_metacyc_entry(metacyc_ligands[component])

        elif component in metacyc_protein_cplxs:
            component_vals = scrub_metacyc_entry(metacyc_protein_cplxs[component])

        if component_vals is None: continue

        if 'Protein-Complexes' in component_vals['TYPES']:
            complex_component = get_or_create_metacyc_protein_complex(session, base, components, genome, metacyc_proteins[component])

        elif 'Protein-Small-Molecule-Complexes' in component_vals['TYPES']:
            complex_component = get_or_create_metacyc_protein_complex(session, base, components, genome, metacyc_protein_cplxs[component])

        elif 'Polypeptides' in component_vals['TYPES']:
            complex_component = get_protein_with_metacyc(session, base, components, genome, metacyc_proteins[component])

        elif 'Compounds' in component_vals['TYPES']:
            complex_component = get_or_create_metacyc_ligand(session, base, components, metacyc_ligands[component])

        if complex_component is None: continue

        session.get_or_create(components.ComplexComposition, complex_id=protein_complex.id,\
                                                                 component_id=complex_component.id,\
                                                                 stoichiometry=1.)
    return protein_complex


def get_or_create_metacyc_transcription_unit(session, base, components, genome, tu_entry):

    vals = scrub_metacyc_entry(tu_entry,extra_args=['COMPONENTS', 'COMMON-NAME'])
    if vals is None:
        vals = scrub_metacyc_entry(tu_entry,extra_args=['COMPONENTS'])
    if vals is None: return None

    genes = []
    tss = None
    for metacyc_component in vals['COMPONENTS']:
        try:
            gene_entry = metacyc_genes[metacyc_component]
            genes.append(get_gene_with_metacyc(session, base, components, genome, gene_entry))
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
        tu = session.get_or_create(components.TU, name=name, leftpos=tss, rightpos=rightpos, strand=strand, genome_id=genome.id)
    else:
        tu = session.get_or_create(components.TU, name=name, leftpos=leftpos, rightpos=tss, strand=strand, genome_id=genome.id)

    for gene in genes:
        session.get_or_create(components.TUGenes, tu_id=tu.id, gene_id=gene.id)

    return tu


@timing
def load_genbank(genbank_file, base, components):

    from Bio import SeqIO

    session = base.Session()
    try:
      gb_file = SeqIO.read(settings.data_directory+'/annotation/GenBank/'+genbank_file,'gb')
    except:
      print 'Error, bipython cannot parse %s' % (genbank_file)
      return

    #from IPython import embed; embed()
    bioproject_id = ''
    for value in gb_file.dbxrefs[0].split():
        if 'BioProject' in value:
            bioproject_id = value.split(':')[1]

    if not bioproject_id:
        print 'Invalid genbank file %s: does not contain a BioProject ID' % (genbank_file)


    genome = session.query(base.Genome).filter(base.Genome.bioproject_id == bioproject_id).first()

    if not genome:
        ome_genome = {'bioproject_id' : bioproject_id,
                      'organism'      : gb_file.annotations['organism']}
        genome = base.Genome(**ome_genome)
        session.add(genome)
        session.flush()


    ome_chromosome = {'genome_id': genome.id,
                      'genbank_id': gb_file.annotations['gi'],
                      'ncbi_id': gb_file.id}

    chromosome = base.Chromosome(**ome_chromosome)
    session.add(chromosome)
    session.flush()


    db_xref_data_source_id = {data_source.name:data_source.id for data_source in session.query(base.DataSource).all()}


    for feature in gb_file.features[0:5]:
        ome_gene = {'long_name':''}
        ome_protein = {'long_name':''}


        if feature.type == 'CDS':

            locus_id = ''
            gene_name = ''

            if 'locus_tag' in feature.qualifiers:
                locus_id = feature.qualifiers['locus_tag'][0]

            if 'gene' in feature.qualifiers:
                gene_name = feature.qualifiers['gene'][0]


            if not gene_name and locus_id:
                gene_name = locus_id
            elif gene_name and not locus_id:
                locus_id = gene_name
            elif not gene_name and not locus_id:
                continue



            ome_gene['locus_id'] = locus_id
            ome_gene['name'] = gene_name
            ome_gene['leftpos'] = int(feature.location.start)
            ome_gene['rightpos'] = int(feature.location.end)
            ome_gene['chromosome_id'] = chromosome.id

            if feature.strand == 1: ome_gene['strand'] = '+'
            elif feature.strand == -1: ome_gene['strand'] = '-'

            if 'product' in feature.qualifiers:
                ome_gene['info'] = feature.qualifiers['product'][0][0:300]
            elif 'note' in feature.qualifiers:
                ome_gene['info'] = feature.qualifiers['note'][0][0:300]

            if 'function' in feature.qualifiers:
                ome_gene['info'] = ome_gene['info'] + ',' + feature.qualifiers['function'][0]
                ome_protein['long_name'] = feature.qualifiers['function'][0]

            if len(ome_gene['name']) > 15: continue  #some weird genbank names are too long and misformed



            gene = session.get_or_create(components.Gene, **ome_gene)



            if 'db_xref' in feature.qualifiers:
                for ref in feature.qualifiers['db_xref']:
                    splitrefs = ref.split(':')
                    ome_synonym = {'type':'gene'}
                    ome_synonym['ome_id'] = gene.id
                    ome_synonym['synonym'] = splitrefs[1]

                    try: data_source_id = db_xref_data_source_id[splitrefs[0]]
                    except:
                        data_source = base.DataSource(name=splitrefs[0])
                        session.add(data_source)
                        session.flush()
                        data_source_id = data_source.id
                        db_xref_data_source_id[splitrefs[0]] = data_source_id


                    ome_synonym['synonym_data_source_id'] = data_source_id

                    synonym = base.Synonyms(**ome_synonym)
                    session.add(synonym)


            if 'gene_synonym' in feature.qualifiers:

                for ref in feature.qualifiers['gene_synonym']:
                    syn = ref.split(';')
                    for value in syn:
                        ome_synonym = {'type': 'gene'}
                        ome_synonym['ome_id'] = gene.id
                        ome_synonym['synonym'] = value

                        ome_synonym['synonym_data_source_id'] = None
                        synonym = base.Synonyms(**ome_synonym)
                        session.add(synonym)

            if 'product' in feature.qualifiers and feature.type == 'CDS':

                try: ome_protein['name'] = feature.qualifiers['protein_id'][0]  #if no protein_id
                except: continue                                                #don't make a protein entry
                ome_protein['gene_id'] = gene.id

                session.get_or_create(components.Protein, **ome_protein)

    session.commit()
    session.close()


@timing
def load_genomes(base, components):
    for genbank_file in open(settings.data_directory+'/annotation/genbanklist.txt','r').readlines():
        genbank_file = genbank_file.rstrip('\n')

        #if genbank_file not in ['NC_000913.2.gb']: continue
        load_genbank(genbank_file, base, components)


@timing
def load_motifs(base, components):
    motifs = open(settings.data_directory + '/annotation/ec_annotation_from_metacyc_2010July19_wpseudo.gff','r')


@timing
def load_metacyc_proteins(base, components, genome):
    session = base.Session()

    for unique_id,entry in metacyc_proteins.iteritems():

        vals = scrub_metacyc_entry(entry)
        if vals is None: continue

        if 'Protein-Complexes' in vals['TYPES'] or 'Protein-Small-Molecule-Complexes' in vals['TYPES']:
            get_or_create_metacyc_protein_complex(session, base, components, genome, entry)

        elif 'Polypeptides' in vals['TYPES']:
            update_protein_with_metacyc(session, base, components, genome, entry)

    session.close()

@timing
def load_metacyc_protein_cplxs(base, components, genome):
    session = base.Session()

    for unique_id,entry in metacyc_protein_cplxs.iteritems():

        vals = scrub_metacyc_entry(entry)
        if vals is None: continue

        if 'Protein-Complexes' or 'Protein-Small-Molecule-Complexes' in vals['TYPES']:
            get_or_create_metacyc_protein_complex(session, base, components, entry)



@timing
def load_metacyc_transcription_units(base, components, genome):

    session = base.Session()
    ##First load annotation file containing merger of metacyc and NCBI
    metacyc_ID = session.get_or_create(base.DataSource, name="metacyc").id


    for unique_id,entry in metacyc_tus.iteritems():

        get_or_create_metacyc_transcription_unit(session, base, components, genome, entry)

    session.close()


@timing
def load_metacyc_bindsites(base, components, genome):
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
                                  rightpos=rightpos, strand='+', genome_id=genome.id, centerpos=centerpos, width=length)


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
def load_kegg_pathways(base, components):

    session = base.Session()
    kegg_pathways = open(settings.data_directory+'/annotation/KEGG/ecoli_KEGG_pathways.txt','r')
    kegg_dataset = session.get_or_create(base.DataSource, name="kegg")
    for line in kegg_pathways.readlines():
        vals = line.rstrip('\r\n').split('\t')
        bnums = vals[2].split(',')
        pathway_ID = vals[0]
        pathway_name = vals[1].strip()
        kegg_pathway = session.get_or_create(components.GeneGroup, name = pathway_name)
        for bnum in bnums:
            gene = session.query(components.Gene).filter_by(locus_id=bnum).first()
            if not gene: continue
            session.get_or_create(components.GeneGrouping, group_id = kegg_pathway.id, gene_id = gene.id)
    session.commit()
    session.close()


@timing
def load_regulatory_network(base, components, data, genome):

    data.RegulatoryNetwork.__table__.drop()
    data.RegulatoryNetwork.__table__.create()


    session = base.Session()
    regulatory_network = open(settings.data_directory+'/annotation/regulondb_gene_reg_network.txt','r')
    for line in regulatory_network.readlines():
        vals = line.split('\t')


        reg_gene = session.query(components.Gene).filter(and_(func.lower(components.Gene.name) == vals[0].lower(),
                                                              components.Gene.genome_id == genome.id)).first()

        regd_gene = session.query(components.Gene).filter(and_(func.lower(components.Gene.name) == vals[1].lower(),
                                                               components.Gene.genome_id == genome.id)).first()

        if reg_gene is None or regd_gene is None: continue

        session.get_or_create(data.RegulatoryNetwork, reg_gene_id=reg_gene.id, regd_gene_id=regd_gene.id,
                                                      direction=vals[2],
                                                      evidence=vals[3])


    session.flush()
    session.commit()
    session.close()



@timing
def write_genome_annotation_gff(base, components, genome):
    session = base.Session()

    genbank_fasta_string = 'gi|'+genome.genbank_id+'|ref|'+genome.ncbi_id+'|'

    with open(settings.data_directory+'/annotation/'+genome.ncbi_id+'.gff', 'wb') as gff_file:

        for gene in session.query(components.Gene).filter(components.Gene.genome_id == genome.id).all():

            info_string = 'gene_id "%s"; transcript_id "%s"; gene_name "%s";' % (gene.locus_id, gene.locus_id, gene.name)

            gff_string = '%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s\n' % (genbank_fasta_string, 'ome_db', 'exon', gene.leftpos,
                                                                                                         gene.rightpos,
                                                                                                         gene.strand,
                                                                                                         info_string)
            gff_file.write(gff_string)

    session.close()



metacyc_genes = parse_metacyc_dat('genes.dat')
metacyc_promoters = parse_metacyc_dat('promoters.dat')
metacyc_proteins = parse_metacyc_dat('proteins.dat')
metacyc_ligands = parse_metacyc_dat('compounds.dat')
metacyc_protein_cplxs = parse_metacyc_dat('protligandcplxes.dat')
metacyc_tus = parse_metacyc_dat('transunits.dat')



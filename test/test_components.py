from nose.tools import assert_equal
from nose.tools import assert_not_equal
from nose.tools import assert_raises

import os
import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base

class TestBase:
    
    def setup(self):
        metadata.create_all()
        
    def teardown(self):
        metadata.drop_all()
        
    @classmethod
    def setup_class(cls):
        os.system('createuser -s -e test_user')
        os.system('createdb test_db --u test_user')
        engine = sqlalchemy.create_engine("postgresql://test_user@localhost/test_db")
        Base = declarative_base(bind=engine)
        metadata = sqlalchemy.MetaData(bind=engine, schema='test_schema')
        
        from om.orm import base
        
    @classmethod
    def teardown_class(cls):
        os.system('dropdb test')
        os.system('dropuser test_user')
        
    def test_genome_region_init(self):
        leftpos = 1000
        rightpos = 2000
        strand = '+'
        gr = base.GenomeRegion(leftpos, rightpos, strand)
        
    def test_chip_experiment(self):
        exp_name = 'chipExo-ArcA_ArcA8myc_fructose_NH4Cl_anaerobic_2_anti-myc_S3_L001_R1_001_sorted.bam'
        lab = 'palsson'
        institution = 'UCSD'
        processing_type = ''
        vals = exp_name.split('_')
        exp_type = vals[0].split('-')
        
        strain = ome.get_or_create(data.Strain, name=vals[1])
        data_source = ome.get_or_create(data.DataSource, name=vals[0], lab=lab, institution=institution)
        environment = ome.get_or_create(data.InVivoEnvironment, name='_'.join(vals[2:5]), carbon_source=vals[2],\
                                        nitrogen_source=vals[3], electron_acceptor=vals[4], temperature=37)
        
        ome.get_or_create(data.ChIPExperiment, name='_'.join(vals[0:6])+'_'+processing_type, replicate=vals[5],\
                                       strain=strain, data_source=data_source, environment=environment,\
                                       protocol_type=exp_type[0], antibody=vals[6].rstrip(processing_type+'.gff'),\
                                       target=vals[0].split('-')[1])

        
        
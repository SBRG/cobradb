from nose.tools import assert_equal
from nose.tools import assert_not_equal
from nose.tools import assert_raises

import os
import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base

class TestBase:
    
    def setup(self):
        self.metadata.create_all()
        
    def teardown(self):
        self.metadata.drop_all()
        
    @classmethod
    def setup_class(cls):
        os.system('createuser -s -e test_user')
        os.system('createdb test_db --u test_user')
        engine = sqlalchemy.create_engine("postgresql://test_user@localhost/test_db")
        cls.Base = declarative_base(bind=engine)
        cls.metadata = sqlalchemy.MetaData(bind=engine)
        
        from om import base
        cls.GenomeRegion = base.GenomeRegion
        cls.Session = base.Session
            
    @classmethod
    def teardown_class(cls):
        os.system('dropdb test')
        os.system('dropuser test_user')
        
    def test_genome_region_init(self):
        leftpos = 1000
        rightpos = 2000
        strand = '+'
        gr = self.GenomeRegion(leftpos, rightpos, strand)
        
    def test_chip_experiment(self):
        exp_name = 'chipExo-ArcA_ArcA8myc_fructose_NH4Cl_anaerobic_2_anti-myc_S3_L001_R1_001_sorted.bam'
        lab = 'palsson'
        institution = 'UCSD'
        processing_type = ''
        vals = exp_name.split('_')
        exp_type = vals[0].split('-')
        session = self.Session()
        


        
        
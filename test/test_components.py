from nose.tools import assert_equal
from nose.tools import assert_not_equal
from nose.tools import assert_raises
from om.orm import base

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
        engine = base.make_engine('test_user', '', 'localhost', 'test_db')
        Base = declarative_base(bind=engine)
        cls.metadata = sqlalchemy.MetaData(bind=engine)
        cls.Session = base.sessionmaker(bind=engine, class_=base._Session)
        
    @classmethod
    def teardown_class(cls):
        os.system('dropdb test')
        os.system('dropuser test_user')
        
    def test_genome_region_init(self):
        leftpos = 1000
        rightpos = 2000
        strand = '+'
        gr = base.GenomeRegion(leftpos, rightpos, strand)
        
    def test_genome_region(self):
        ome = self.Session()
        gr = ome.get_or_create(base.GenomeRegion, leftpos=1000, rightpos=2000, strand='+')
        

        
        
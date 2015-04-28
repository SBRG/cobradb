# -*- coding: utf-8 -*-

from ome.loading.map_loading import load_the_map
from ome import base

import pytest

def test_load_the_map(test_db):
    session = base.Session()
    
    assert load_the_map(None, None, None, 'x'*500000) == 1
    with pytest.raises(Exception):
        load_the_map(None, None, None, 'x'*400000)

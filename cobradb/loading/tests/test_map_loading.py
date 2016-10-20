# -*- coding: utf-8 -*-

from cobradb.loading.map_loading import load_the_map
from cobradb import base

import pytest

def test_load_the_map(test_db):
    session = base.Session()

    assert load_the_map(None, None, None, 'x'*1000000) == 1
    with pytest.raises(Exception):
        load_the_map(None, None, None, 'x'*400000)

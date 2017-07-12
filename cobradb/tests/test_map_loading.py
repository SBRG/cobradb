# -*- coding: utf-8 -*-

from cobradb.models import Session
from cobradb.map_loading import load_the_map

import pytest

def test_load_the_map(test_db):
    session = Session()

    assert load_the_map(None, None, None, 'x'*1000000) == 1
    with pytest.raises(Exception):
        load_the_map(None, None, None, 'x'*400000)

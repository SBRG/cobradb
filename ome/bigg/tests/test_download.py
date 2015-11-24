# -*- coding: utf-8 -*-

from ome.bigg.download import *

def test_download_model():
    model = download_model('iJO1366')
    assert model.id == 'iJO1366'
    assert model.optimize().f > 0.1

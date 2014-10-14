#Installation of the OME Framework

## Development version

```
git clone https://github.com/SBRG/ome.git
cd ome
python setup.py develop
```

####After adding data to ome/data
```
python bin/load_db.py
```

####Dependencies
The OME Framework is dependent on [cobrapy](https://github.com/opencobra/cobrapy/blob/master/INSTALL.md) and may require additional installation.

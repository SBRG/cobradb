#Installation of the OME Framework

## Development version

```
git clone https://github.com/SBRG/ome.git
cd ome
pip install -e .
```

####After adding data to ome/data
```
python bin/load_db.py
```

####Dependencies
* [cobrapy](https://github.com/opencobra/cobrapy/blob/master/README.md) for which you may want to refer to the [installation docs](https://github.com/opencobra/cobrapy/blob/master/INSTALL.md).
* [pysam](https://github.com/pysam-developers/pysam) which depends on [samtools](
* 

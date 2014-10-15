#Installation of the OME Framework

## Development version

#### Clone the repository and setup in develop mode
```
git clone https://github.com/SBRG/ome.git
cd ome
sudo python setup.py develop
```

####After adding data to ome/data
```
python bin/load_db.py
```

####Dependencies
While all of the required dependencies *should* install automatically through setup.py above in practice the following generally need to be installed individually.
* [cobrapy](https://github.com/opencobra/cobrapy/blob/master/README.md) for which you may want to refer to the [installation docs](https://github.com/opencobra/cobrapy/blob/master/INSTALL.md).
* [pysam](https://github.com/pysam-developers/pysam) which depends on samtools
* postgres

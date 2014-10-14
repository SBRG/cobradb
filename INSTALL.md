#Installation of the OME Framework

## Development version

#### Clone the repository and setup in develop mode
```
git clone https://github.com/SBRG/ome.git
cd ome
sudo pip install -e .
```

#### Rename the example settings file and edit with appropriate values
```
mv ome/settings.ini.example ome/settings.ini
```
This currently involves a working postgres installation.

####After adding data to ome/data
```
python bin/load_db.py
```

####Dependencies
* [cobrapy](https://github.com/opencobra/cobrapy/blob/master/README.md) for which you may want to refer to the [installation docs](https://github.com/opencobra/cobrapy/blob/master/INSTALL.md).
* [pysam](https://github.com/pysam-developers/pysam) which depends on samtools
* postgres

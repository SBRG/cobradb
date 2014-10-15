#Installation of the OME Framework

## Development version

#### Clone the repository and setup in develop mode
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
All of the required dependencies *should* install automatically through setup.py above.  However, in practice the following generally need to be installed individually.
* [cobrapy](https://github.com/opencobra/cobrapy/blob/master/README.md) for which you may want to refer to the [installation docs](https://github.com/opencobra/cobrapy/blob/master/INSTALL.md).
* [pysam](https://github.com/pysam-developers/pysam) which depends on [samtools](http://samtools.sourceforge.net/)

####Optional Dependencies
If you intend to load more than a few models and any significant quantity of expression or annotation data then you should install postgres and mongodb.

#####On OSX with homebrew http://brew.sh/
```
brew install postgresql
brew install mongodb
```

Setting up the database
```
createdb ome
createuser -d -l -s
```

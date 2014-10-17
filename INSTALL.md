#Installation of the OME Framework

## Development version

####Dependencies
The current version of OME is dependent on postgresql and mongodb. In the future it will also run natively on SQLlite making postgresql and mongodb optional and the overall setup much simpler.

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

All of the rest of the dependencies *should* install automatically through setup.py above.  However, in practice the following generally also need to be installed individually.
* [cobrapy](https://github.com/opencobra/cobrapy/blob/master/README.md) for which you may want to refer to the [installation docs](https://github.com/opencobra/cobrapy/blob/master/INSTALL.md).
* [pysam](https://github.com/pysam-developers/pysam) which depends on [samtools](http://samtools.sourceforge.net/)
* [numpy](http://www.numpy.org/), [scipy](http://www.scipy.org/), and [pandas](http://pandas.pydata.org/) which depending on your OS and configuration can be non-trivial

#####On ubuntu
```
sudo apt-get install python-dev zlib1g-dev samtools g++ libblas-dev liblapack-dev gfortran
sudo pip install cython pysam numpy scipy pandas
```

#### Clone the repository and setup in develop mode
```
git clone https://github.com/SBRG/ome.git
cd ome
python setup.py develop
python bin/load_db.py
```

This will create a directory **ome_data** in your home folder by default.  You can move this folder and change the location by altering the value in **ome/ome/settings.ini**




####Optional Dependencies
If you intend to load more than a few models and any significant quantity of expression or annotation data then you should install postgres and mongodb.


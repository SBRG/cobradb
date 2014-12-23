#Installation of the OME Framework

## Development version

####Dependencies
The current version of OME is dependent on postgresql, mongodb, and a number of python packages. In the future it will run natively on SQLlite making postgresql and mongodb optional and the overall setup much simpler.

######First install and configure the databases.


#####On OSX with homebrew http://brew.sh/
```
brew install postgresql mongodb
```

#####On Ubuntu
```
sudo apt-get install postgresql postgresql-contrib postgresql-server-dev-all mongodb
```
#####Setting up the database

```
sudo -i -u postgres                    
createuser -d -l -P -s <your username>
createdb ome
exit
```
Note: newer versions of OSX only require ```createdb ome```

######Next install additional python packages.

All of the rest of the dependencies *should* install automatically through setup.py below.  However, in practice you may want to install these individually ahead of time.
* [cobrapy](https://github.com/opencobra/cobrapy/blob/master/README.md) for which you may want to refer to the [installation docs](https://github.com/opencobra/cobrapy/blob/master/INSTALL.md).
* [pysam](https://github.com/pysam-developers/pysam) which depends on [samtools](http://samtools.sourceforge.net/)
* [numpy](http://www.numpy.org/), [scipy](http://www.scipy.org/), and [pandas](http://pandas.pydata.org/) which depending on your OS and configuration can be non-trivial


#####On OSX with homebrew http://brew.sh/
```
brew install samtools glpk
```

#####On Ubuntu
```
sudo apt-get install python-dev zlib1g-dev samtools g++ libblas-dev liblapack-dev gfortran
```
#####With pip
```
sudo pip install cython pysam numpy scipy pandas cobra ipython psycopg2 pymongo
```

### Finally, clone the repository and setup in develop mode
```
git clone https://github.com/SBRG/ome.git
cd ome
python setup.py develop
python bin/load_db.py
```

This will create a directory **ome_data** in your home folder by default.  You can move this folder and change the location by altering the value in **ome/ome/settings.ini**



#Installation of the OME Framework

## Dependencies

The current version of OME is dependent on postgresql, mongodb, and a number of
python packages. In the future it will run natively on SQLlite making postgresql
and mongodb optional and the overall setup much simpler.

## On OSX with homebrew http://brew.sh/

```
brew install samtools glpk
brew install postgresql 
```

## On Ubuntu

```
sudo apt-get install postgresql postgresql-contrib postgresql-server-dev-all 
sudo apt-get install python-dev zlib1g-dev samtools g++ libblas-dev liblapack-dev gfortran
```

## Setting up the database

```
sudo -i -u postgres
createuser -d -l -s <your username>
createdb ome
exit
```

## Next install python packages.

All of the rest of the dependencies *should* install automatically through
setup.py below.  However, in practice you may want to install these individually
ahead of time.

* [cobrapy](https://github.com/opencobra/cobrapy/blob/master/README.md) for which you may want to refer to the [installation docs](https://github.com/opencobra/cobrapy/blob/master/INSTALL.md).

For SBML3 output, need to manually install https://github.com/aebrahim/cobrapy/commit/242949bc313b696e1eabf1609730cf20dc7983c0
```
pip install cython
git clone git@github.com:aebrahim/cobrapy.git cobrapy-sbml3
cd cobrapy-sbml3
git checkout 242949bc313b696e1eabf1609730cf20dc7983c0
python setup.py develop
```

* [pysam](https://github.com/pysam-developers/pysam) which depends on [samtools](http://samtools.sourceforge.net/)
* [numpy](http://www.numpy.org/), [scipy](http://www.scipy.org/), and [pandas](http://pandas.pydata.org/) which depending on your OS and configuration can be non-trivial

## Finally, clone the repository and install

```
git clone https://github.com/SBRG/ome.git
cd ome
python setup.py install
# OR
python setup.py develop
python bin/load_db.py --drop-all
```

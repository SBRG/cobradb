# COBRAdb

COBRAdb loads genome-scale metabolic models and genome annotations into a
relational database. It is the software that powers
[BiGG Models](http://bigg.ucsd.edu).

## Installation

COBRAdb requires PostgreSQL, Python 2.7, and a number of Python packages listed
in `setup.py`.

1. Install and set up PostgreSQL. This process varies from platform to platform,
   but the best place to start in the PostgreSQL
   [documentation](https://www.postgresql.org/docs/).

2. Set up your configuration. Copy the file `settings.ini.example` and rename it
   `settings.ini`. The example file includes descriptions of each setting.

3. Install OME and dependencies.

```
python setup.py install
# OR
python setup.py develop
# OR
pip
```

4. Load the database by calling the script `bin/load_db`.

```
cd bin
./load_db --drop-all
```

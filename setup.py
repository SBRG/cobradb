# -*- coding: utf-8 -*-

from os.path import abspath, dirname
from sys import path

# To temporarily modify sys.path
SETUP_DIR = abspath(dirname(__file__))

try:
    from setuptools import setup, find_packages
except ImportError:
    path.insert(0, SETUP_DIR)
    import ez_setup
    path.pop(0)
    ez_setup.use_setuptools()
    from setuptools import setup, find_packages


setup(
    name='cobradb',
    version='0.2.0',
    description="""COBRAdb loads genome-scale metabolic models and genome
                   annotations into a relational database.""",
    url='https://github.com/SBRG/cobradb',
    author='Zachary King',
    author_email='zaking@ucsd.edu',
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
    ],
    keywords='systems biology, genome-scale model',
    packages=find_packages(),
    install_requires=['SQLAlchemy>=1.0.12',
                      'cobra>=0.4.0',
                      'numpy>=1.9.1',
                      'psycopg2>=2.5.4',
                      'biopython>=1.65',
                      'scipy>=0.17.0',
                      'lxml>=3.6.0',
                      'pytest>=2.6.4',
                      'six>=1.10.0',
                      'tornado>=4.4.2',
                      'escher>=1.5.0',
                      'configparser>=3.5.0'],
)

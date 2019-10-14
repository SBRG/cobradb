# -*- coding: utf-8 -*-

from os.path import abspath, dirname
from sys import path

# To temporarily modify sys.path
SETUP_DIR = abspath(dirname(__file__))

from setuptools import setup, find_packages


setup(
    name='cobradb',
    version='0.3.0',
    description="""COBRAdb loads genome-scale metabolic models and genome
                   annotations into a relational database.""",
    url='https://github.com/SBRG/cobradb',
    author='Zachary King',
    author_email='zaking@ucsd.edu',
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='systems biology, genome-scale model',
    packages=find_packages(),
    install_requires=[
        'SQLAlchemy>=1.3.10,<2',
        'cobra>=0.16.0,<0.17',
        'python-libsbml>=5.18.0,<6',
        'numpy>=1.17.2,<2',
        'psycopg2>=2.8.3,<3',
        'biopython>=1.74,<2',
        'scipy>=1.3.1,<2',
        'lxml>=4.4.1,<5',
        'pytest>=4.6.6,<5',
        'six>=1.12.0,<2',
        'tornado>=4.5.3,<5',
        'escher>=1.7.3,<2',
        'configparser>=4.0.2,<5',
    ],
)

# -*- coding: utf-8 -*-

try:
    from setuptools import setup, Command
except:
    from distutils.core import setup, Command

setup(
    name='ome',
    version='0.2.0',
    description="""COBRAdb loads genome-scale metabolic models and genome
                   annotations into a relational database.""",
    url='https://github.com/SBRG/ome',
    author='Zachary King',
    author_email='zaking@ucsd.edu',
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
    keywords='systems biology, genome-scale model',
    packages=['ome'],
    install_requires=['SQLAlchemy>=1.0.12',
                      'cobra>=0.4.0',
                      'numpy>=1.9.1',
                      'psycopg2>=2.5.4',
                      'biopython>=1.65',
                      'scipy>=0.17.0',
                      'lxml>=3.6.0',
                      'pytest>=2.6.4'],
)

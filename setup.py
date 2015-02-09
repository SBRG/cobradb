from codecs import open  # To use a consistent encoding
from os import path

try:
    from setuptools import setup, Command
except:
    from distutils.core import setup, Command

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='ome',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # http://packaging.python.org/en/latest/tutorial.html#version
    version='0.0.1-bigg',

    description='The OME Framework',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/SBRG/ome',

    # Author details
    author='Steve Federowicz, Ali Ebrahim, Zak King, Justin Lu, Joshua Lerman, Edward OBrien',
    author_email='sfederow@gmail.com',

    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Systems Biology :: Data Analysis',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
    keywords='systems biology',
    packages=['ome'],
    install_requires=['SQLAlchemy>=0.9.8',
                      'cobra>=0.3.0',
                      'escher>=1.0.0', # could be optional
                      'numpy>=1.9.1',
                      'pandas>=0.15.2', # could be optional
                      'pysam>=0.8.1', # could be optional
                      'pyzmq>=14.4.1', # could be optional
                      'scipy>=0.15.0', # could be optional
                      'simplejson>=3.6.5',
                      'psycopg2>=2.5.4',
                      'biopython>=1.65'], # could be optional
    extras_require = {'all': ['pymongo>=2.7.2',
                              'ipython>=2.3.1',
                              'pytest>=2.6.4']},
    package_data={'ome': ['data/*']})

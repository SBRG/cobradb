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
    version='0.1.0',

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
    install_requires=['SQLAlchemy>=1.0.12',
                      'cobra>=0.4.0',
                      'numpy>=1.9.1',
                      'psycopg2>=2.5.4',
                      'biopython>=1.65',
                      'scipy>=0.17.0',
                      'lxml>=3.6.0',
                      'pytest>=2.6.4'],
    package_data={'ome': ['data/*']})

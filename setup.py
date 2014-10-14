from distutils.core import setup

setup(
    name='OME Framework',
    version='0.1dev',
    packages=['ome',],
    package_data={'ome': ['data/*']},
    install_requires = ['sqlalchemy','psycopg2','numpy','scipy','cobra','pymongo','simplejson','pysam'],
    license='MIT License',
    long_description=open('README.md').read(),
)

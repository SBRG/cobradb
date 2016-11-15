# -*- coding: utf-8 -*-

from ome.base import *
from ome import settings

import re
import os
import logging


def check_none(v):
    """Return None if v is the empty string or the string 'None'."""
    return None if (v == 'None' or v == '') else v


def get_or_create(session, query_class, **kwargs):
    """Query the query_class, filtering by the given keyword arguments. If no result
    is found, then add a new row to the database.

    Returns a tuple: (result of the query, Boolean: the result already existed)

    Arguments
    ---------

    session: The SQLAlchemy session.

    query_class: The class to query.

    """
    res = session.query(query_class).filter_by(**kwargs).first()
    if res is not None:
        return res, True
    res = query_class(**kwargs)
    session.add(res)
    session.commit()
    return res, False


def load_tsv(filename, required_column_num=None):
    """Try to load a tsv prefs file. Ignore empty lines and lines beginning with #.

    Arguments
    ---------

    filename: A tsv path to load.

    required_column_num: The number of columns to check for.

    """
    if not os.path.exists(filename):
        return []
    with open(filename, 'r') as f:
        # split non-empty rows by tab
        rows = [[check_none(x.strip()) for x in line.split('\t')]
                for line in f.readlines()
                if line.strip() != '' and line[0] != '#']

    # check rows
    if required_column_num is not None:
        def check_row(row):
            if len(row) != required_column_num:
                logging.warn('Line in {} should have {} columns, but found {}: {}'
                             .format(filename, required_column_num, len(row), row))
                return None
            return row
        rows = [x for x in (check_row(r) for r in rows) if x is not None]

    return rows


def _find_data_source_url(bigg_id, url_prefs):
    """Return (bigg_id, name, url prefix) for data source name."""
    name = None; url_prefix = None
    for row in url_prefs:
        if row[0] == bigg_id:
            if len(row) == 1:
                logging.error('Bad row in data-source-prefs: {}'.format(row))
            if len(row) > 1:
                name = check_none(row[1])
            if len(row) > 2:
                url_prefix = check_none(row[2])
            break
        # check synonyms
        elif len(row) == 4 and bigg_id in (x.strip() for x in row[3].split(',')):
            bigg_id, name, url_prefix = row[0], check_none(row[1]), check_none(row[2])
            break
    return bigg_id, name, url_prefix


def get_or_create_data_source(session, bigg_id):
    """Get the data source by name. If it does not exist in the database, then
    add a new row by reading the preference file.

    Arguments
    ---------

    session: The SQLAlchemy session.

    bigg_id: The BiGG ID of the DataSource.

    """
    data_source_db = (session
                      .query(DataSource)
                      .filter(DataSource.bigg_id == bigg_id)
                      .first())
    if not data_source_db:
        # get gene url_prefs
        url_prefs = load_tsv(settings.data_source_preferences)
        bigg_id, name, url_prefix = _find_data_source_url(bigg_id, url_prefs)
        # data source may already exist if this is a synonym
        data_source_db, exists = get_or_create(session, DataSource,
                                               bigg_id=bigg_id,
                                               name=name,
                                               url_prefix=url_prefix)
        if not exists:
            if name is None:
                logging.warn('No name found for data source %s' % bigg_id)
            if url_prefix is None:
                logging.warn('No URL found for data source %s' % bigg_id)
    return data_source_db.id


def increment_id(id, increment_name=''):
    match = re.match(r'(.*)_%s([0-9]+)$' % increment_name, id)
    if match:
        return '%s_%s%d' % (match.group(1), increment_name, int(match.group(2)) + 1)
    else:
        return '%s_%s%d' % (id, increment_name, 1)


def make_reaction_copy_id(bigg_id, copy_number):
    return '{}_copy{}'.format(bigg_id, copy_number)


def check_pseudoreaction(reaction_id):
    patterns = [
        r'^ATPM$',
        r'^EX_.*',
        r'^DM_.*',
        r'^SK_.*',
        r'^BIOMASS_.*' # case insensitive
    ]
    for pattern in patterns:
        if re.match(pattern, reaction_id):
            return True
    return False


def format_formula(formula):
    """Remove unnecessary characters from formula."""
    if formula is None:
        return formula
    formula = formula.strip("'[]").replace('FULLR', 'R')
    if formula.strip() in ['Fe1SO', 'UK', 'HP']:
        return 'R'
    return formula


def scrub_gene_id(the_id):
    """Get a new style gene ID."""
    the_id = re.sub(r'(.*)\.([0-9]{1,2})$', r'\1_AT\2', the_id)
    the_id = re.sub(r'\W', r'_', the_id)
    return the_id


def scrub_name(the_name):
    """Make a nice looking name."""
    if the_name is None:
        return None
    the_name = (the_name
                .replace('_SPACE_SPACE_', ' ')
                .replace('_SPACE_', ' ')
                .replace('_COLON_', ':')
                .replace('_COMMA_', ','))
    the_name = re.sub(r'^[RMG]?_', '', the_name)
    the_name = re.sub(r'_', ' ', the_name)
    # uppercase
    the_name = re.sub('^([a-z])', lambda x: x.group(1).upper(), the_name)
    if the_name.strip() == '':
        return None
    return the_name


def ref_str_to_tuple(ref):
    """String like ' a : b ' to tuple like ('a', 'b')."""
    return tuple(x.strip() for x in ref.split(':'))


def ref_tuple_to_str(key, val):
    """Tuple like ('a', 'b') to string like 'a:b'."""
    return '%s:%s' % (key, val)

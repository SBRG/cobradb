# -*- coding: utf-8 -*-

from ome.base import *
from ome import settings

import re
import os
import logging


def check_none(v):
    return None if (v == 'None' or v == '') else v


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


def find_data_source_url(a_name, url_prefs):
    """Return the url prefix for data source name, or None."""
    for row in url_prefs:
        if row[0] == a_name:
            return row[1]
    return None


def get_or_create_data_source(session, data_source_name):
    data_source_db = (session
                      .query(DataSource)
                      .filter(DataSource.name == data_source_name)
                      .first())
    if not data_source_db:
        # get gene url_prefs
        url_prefs = load_tsv(settings.data_source_preferences)
        url_prefix = find_data_source_url(data_source_name, url_prefs)
        if url_prefix is None:
            logging.warn('No URL found for data source %s' % data_source_name)
        data_source_db = DataSource(name=data_source_name,
                                    url_prefix=url_prefix)
        session.add(data_source_db)
        session.flush()

    return data_source_db.id


def format_formula(formula):
    """Remove unnecessary characters from formula."""
    if formula is None:
        return formula
    return formula.strip("'[]")


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


def ref_str_to_tuple(ref):
    """String like ' a : b ' to tuple like ('a', 'b')."""
    return tuple(x.strip() for x in ref.split(':'))


def ref_tuple_to_str(key, val):
    """Tuple like ('a', 'b') to string like 'a:b'."""
    return '%s:%s' % (key, val)

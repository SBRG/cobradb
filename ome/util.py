# -*- coding: utf-8 -*-

import re
import os
import logging


def increment_id(id, increment_name=''):
    match = re.match(r'(.*)_%s([0-9]+)$' % increment_name, id)
    if match:
        return '%s_%s%d' % (match.group(1), increment_name, int(match.group(2)) + 1)
    else:
        return '%s_%s%d' % (id, increment_name, 1)


def check_pseudoreaction(reaction_id):
    patterns = [
        r'^ATPM$', r'^ATPM_NGAM$',
        r'^EX_.*',
        r'^DM.*',
        r'(?i).*biomass.*' # case insensitive
    ]
    for pattern in patterns:
        if re.match(pattern, reaction_id):
            return True
    return False


def scrub_gene_id(the_id):
    """Get a new style gene ID."""
    the_id = re.sub(r'(.*)\.([0-9]{1,2})$', r'\1_AT\2', the_id)
    the_id = re.sub(r'\W', r'_', the_id)
    return the_id


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
        rows = [[x.strip() for x in line.split('\t')]
                for line in f.readlines()
                if line.strip() != '' and line[0] != '#']
    # check rows
    if required_column_num is not None:
        def check_row(row):
            if len(row) != required_column_num:
                logging.warn('Bad row in gene_reaction_rule_prefs: %s' % row)
                return None
            return row
        rows = [x for x in (check_row(r) for r in rows) if x is not None]

    return rows

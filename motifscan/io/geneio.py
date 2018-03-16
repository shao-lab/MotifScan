"""Module script for gene IO."""

import sys
import logging
import pandas as pd


def load_refseq_gene(input_file):
    """load the RefSeq gene annotation file.

    Args:
        input_file (str): Path of RefSeq gene annotation input file.

    Returns:
        genes (pd.DataFrame): Gene dataframe.

    """

    def define_tss(strand, start, end):
        if strand == '+':
            return start
        elif strand == '-':
            return end
        else:
            logging.ERROR("Invalid gene strand: {}!".format(strand))
            sys.exit(1)

    genes = pd.read_csv(input_file, sep='\t', usecols=[1, 2, 3, 4, 5], names=['id', 'chr', 'strand', 'start', 'end'])
    genes['TSS'] = map(define_tss, genes['strand'], genes['start'], genes['end'])
    return genes

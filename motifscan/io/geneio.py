# Copyright 2014-2017, Hongduo Sun, Jiawei Wang, Zhen Shao
#
# This file is part of MotifScan.
#
# MotifScan is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MotifScan is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MotifScan.  If not, see <http://www.gnu.org/licenses/>.

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

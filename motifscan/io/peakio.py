"""Module script for peak(genomic regions) IO."""

import os
import re
import sys
import logging
import pandas as pd

from motifscan.io.constants import *


def _get_comment_line_num(input_file):
    """Get the number comment and blank line
    
    Args:
        input_file (str): Path of input file.
    
    Returns:
        skip_num (int): Number of comment lines to skip.
        
    """
    skip_num = 0
    with open(input_file) as fin:
        for line in fin:
            tmp_line = line.strip()
            if tmp_line == '' or tmp_line.startswith('#') or not re.match(r'[Cc][Hh][Rr]\S', tmp_line):
                skip_num += 1
            else:
                break
    return skip_num


def load_peak(peak_file, peak_format, peak_length):
    """load peaks from input file.
    
    Args:
        peak_file (str): Path of peak input file.
        peak_format (str): Peak format.
        peak_length (int): Peak length around summit to scan motifs.
    
    Returns:
        peaks (pd.DataFrame): Peak dataframe.
        
    """
    skip_num = _get_comment_line_num(peak_file)  # skip the comment lines
    # load the peaks
    if peak_format == 'macs':
        peaks = pd.read_csv(peak_file, sep='\t', names=MACS_HEADER, skiprows=skip_num)
        peaks['summit'] = peaks['start'] + peaks['summit']
    elif peak_format == 'manorm':
        peaks = pd.read_csv(peak_file, sep='\t', usecols=[0, 1, 2, 3, 4], names=MANORM_HEADER, skiprows=skip_num)
        peaks['summit'] = peaks['start'] + peaks['summit']
    elif peak_format == 'bed3col':
        peaks = pd.read_csv(peak_file, '\t', names=BED3_HEADER, skiprows=skip_num)
        peaks['summit'] = (peaks['start'] + peaks['end']).map(lambda x: int(x / 2))
    elif peak_format == 'bed4col':
        peaks = pd.read_csv(peak_file, '\t', names=BED4_HEADER, skiprows=skip_num)
        peaks['summit'] = peaks['start'] + peaks['summit']
    elif peak_format == 'bed5col':
        peaks = pd.read_csv(peak_file, '\t', names=BED5_HEADER, skiprows=skip_num)
        peaks['summit'] = peaks['start'] + peaks['summit']

    if (peaks['summit'] > peaks['end']).any():
        logging.error("Invalid peak format: 'summit' is greater than 'end'. Please check it!")
        sys.exit(1)

    if peak_length == 0:  # use whole peak region to perform MotifScan
        peaks['seq_start'] = peaks['start']
        peaks['seq_end'] = peaks['end']
    else:  # use specified peak length to perform MotifScan
        peaks['seq_start'] = peaks['summit'] - peak_length / 2
        peaks.loc[peaks['seq_start'] < 0, 'seq_start'] = 0
        peaks['seq_end'] = peaks['summit'] + peak_length / 2
        # TODO chrom_size
    return peaks

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

"""Module script for genome-related operations."""

import re
import logging
import numpy as np
import pandas as pd

ROMAN_CHRS = ('chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX', 'chrXI', 'chrXII',
              'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrXVII', 'chrXVIII', 'chrXIX', 'chrXX', 'chrXXI',
              'chrXXII', 'chrXXIII',
              'ChrI', 'ChrII', 'ChrIII', 'ChrIV', 'ChrV', 'ChrVI', 'ChrVII', 'ChrVIII', 'ChrIX', 'ChrXI', 'ChrXII',
              'ChrXIII', 'ChrXIV', 'ChrXV', 'ChrXVI', 'ChrXVII', 'ChrXVIII', 'ChrXIX', 'ChrXX', 'ChrXXI',
              'ChrXXII', 'ChrXXIII')


class GenomePosError(Exception):
    pass


def compute_bg_freq(genome_file):
    """Compute the genome-wide background frequency of 4 nucleotides, only autosomes are taken into account.
    
    Args:
        genome_file (str): Path of genome sequences input file.
        
    Returns:
        bg_freq (pd.Series): Genome background frequency of nucleotides.
        Example:
            base     frequency
               A    0.29485
               C    0.20491163
               G    0.20500062
               T    0.29523775
        
    """
    nt_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    autosome_flag = True
    roman_name_flag = False
    with open(genome_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if line[0] == '>':
                line = line[1:]
                if re.search(r"[Cc][Hh][Rr][0-9]+\b", line) or line in ROMAN_CHRS or \
                        ((line == 'chrX' or line == 'ChrX') and roman_name_flag):
                    if line in ROMAN_CHRS:
                        roman_name_flag = True
                    logging.info("Processing {}...".format(line))
                    autosome_flag = True
                else:
                    logging.info("Skipping {}...".format(line))
                    autosome_flag = False
            elif autosome_flag:
                for nt in line:
                    nt = nt.upper()
                    try:
                        nt_count[nt] += 1
                    except KeyError:
                        continue
    a, c, g, t = nt_count['A'], nt_count['C'], nt_count['G'], nt_count['T']
    bg_freq_arr = np.array([a, c, g, t]) / float(np.sum([a, c, g, t]))
    bg_freq = pd.Series(bg_freq_arr, index=['A', 'C', 'G', 'T'], name='frequency')
    bg_freq.index.rename('base', inplace=True)
    return bg_freq


def compute_genome_size(genome_file):
    """Compute the size of chromosomes.

     Args:
        genome_file (str): Path of genome sequences input file.

     Returns:
        chr_size (pd.Series):Chromosome size.
        Example of hg19:
             chr       size
            chr1  249250621
            chr2  243199373
            chr3  198022430
            ...

    """
    chrs = []
    chr_size = []
    with open(genome_file, 'r') as fin:
        empty_flag = True
        for line in fin:
            if line[0] == '>':  # chromosome id line
                if empty_flag:
                    empty_flag = False
                else:
                    chrs.append(tmp_chr)
                    chr_size.append(tmp_size)
                tmp_chr = line.strip()[1:]
                tmp_size = 0
            else:  # sequence line
                tmp_size += len(line.strip())
        if not empty_flag:
            chrs.append(tmp_chr)
            chr_size.append(tmp_size)
        chr_size = pd.Series(chr_size, index=chrs, name='size')
        chr_size.index.rename('chr', inplace=True)
        return chr_size


def extract_sequence(genome_dir, chrom, start, end):
    """Extract sequences from the given genomic region.
    Note that the coordinate system is 0-based.

    Args: 
        genome_dir (str): Path of pre-complied genome directory.
        chrom (str): Chromosome name of the genomic region.
        start (int): Start of the genomic region.
        end (int): End of the genomic region.

    Returns:
         Sequence (uppercase) of the given genomic region.

    """
    if start > end:
        raise GenomePosError("Chr:{0}\tStart:{1}\tEnd:{2}".format(chrom, start, end))
    fin = open("{0}/{1}".format(genome_dir, chrom), 'r')
    fin.seek(0, 0)
    fin.readline()  # read the first line; the pointer is at the second line
    length = end - start
    offset = start + int(np.floor(start / 50.0))
    fin.seek(offset, 1)
    tmp_seq = fin.read(length + int(np.ceil(length / 50.0)))
    tmp_seq = tmp_seq.replace('\n', '')
    fin.close()
    return tmp_seq[:length].upper()

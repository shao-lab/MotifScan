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

"""Module script for motif-related operations."""

import ctypes
import logging
from pkg_resources import resource_filename
import numpy as np
import pandas as pd

import motifscan
from motifscan.lib.genome import extract_sequence
from motifscan.lib.matrix import MAT, arr_to_mat, construct_sequence_matrix


def compute_max_score(motifs, bg_freq):
    """Compute the maximal motif score.

    Args:
        motifs (pd.DataFrame): Motif dataframe.
        bg_freq (pd.Series): Genome background frequency of nucleotides.

    Returns:
        motifs (pd.DataFrame): Motif dataframe with maximal motif score.

    """
    log_bg = np.log(bg_freq)
    max_scores = []
    for idx in motifs.index:
        matrix = motifs.loc[idx, 'matrix']
        motif_length = matrix.shape[1]
        max_score = 0
        for tmp_pos in xrange(motif_length):
            tmp_pos_score = np.log(matrix[:, tmp_pos]) - log_bg
            max_score += tmp_pos_score.max()
        max_scores.append(max_score)
        # Bug fix in version 1.1: maximal motif score may not come from the consensus sequence
    motifs['max_score'] = max_scores
    return motifs


def sampling_regions(sampling_times, length, chr_size):
    """Sampling genomic regions from whole genome.

    Args:
        sampling_times (int): Number of sampling.
        length (int): Length of genomic regions.
        chr_size (pd.Series): Chromosome size.

    Returns:
        random_regions (pd.DataFrame): Random sampled regions.

    """
    chrs = []
    starts = []
    ends = []
    chr_weight = chr_size / chr_size.sum()
    #  weighted sampling by chromosome size
    for i in xrange(sampling_times):
        tmp_chr = np.random.choice(chr_size.index, p=chr_weight)
        start = np.random.randint(chr_size[tmp_chr] - length)
        end = start + length
        chrs.append(tmp_chr)
        starts.append(start)
        ends.append(end)
    random_regions = pd.DataFrame({'chr': chrs, 'start': starts, 'end': ends})
    return random_regions


def simulation(motifs, sampling_times, genome_dir, bg_freq, chr_size):
    """Random sampling regions from background and calculate background distribution of each motif score.

    Args:
        motifs (pd.DataFrame): Motif dataframe.
        sampling_times (int): Sampling times.
        genome_dir (str): Path of pre-compiled genome directory.
        bg_freq (pd.Series): Genome background frequency of nucleotides.
        chr_size (pd.Series): Chromosome size.

    Returns:
        motifs (pd.DataFrame): Motif dataframe with motif_score_cutoff.

    """
    score_c = ctypes.CDLL(resource_filename(motifscan.__name__, 'score_c.so'))
    score_c.motif_scan_core_simulation.restype = ctypes.c_double
    score_c.motif_scan_core_simulation.argtypes = [ctypes.POINTER(MAT), ctypes.POINTER(MAT),
                                                   ctypes.POINTER(ctypes.c_double * 4), ctypes.c_double]
    max_length = 0
    for matrix in motifs['matrix']:
        if max_length < matrix.shape[1]:
            max_length = matrix.shape[1]

    try:
        logging.info("Generating coordinates...")
        random_regions = sampling_regions(sampling_times=sampling_times, length=max_length, chr_size=chr_size)
        logging.info("Extracting sequences...")
        genome_dir_iter = [genome_dir] * sampling_times
        random_regions['seq'] = map(extract_sequence, genome_dir_iter, random_regions['chr'],
                                    random_regions['start'], random_regions['end'])
        random_regions['seq_matrix'] = map(construct_sequence_matrix, random_regions['seq'])
        logging.info("Calculating motif score of random sampled regions...")
        bg_freq = np.asarray(bg_freq)
        bg_freq = np.ctypeslib.as_ctypes(bg_freq)
        seq_matrix_col_iloc = random_regions.columns.get_loc('seq_matrix')
        motif_num = len(motifs)
        cnt = 0
        motif_score_all = pd.DataFrame()
        for idx, motif in motifs.iterrows():
            cnt += 1
            logging.info("Scanning motif: {0}\t{1}/{2}...".format(motif['name'], cnt, motif_num))
            max_score = ctypes.c_double(motif['max_score'])
            m_matrix = arr_to_mat(motif['matrix'])
            motif_score = []
            for region_idx in xrange(sampling_times):
                s_matrix = arr_to_mat(random_regions.iat[region_idx, seq_matrix_col_iloc])
                tmp_score = score_c.motif_scan_core_simulation(ctypes.byref(s_matrix), ctypes.byref(m_matrix),
                                                               ctypes.byref(bg_freq), max_score)
                motif_score.append(tmp_score)
            motif_score.sort(reverse=True)
            motif_score_all[motif['id']] = motif_score

        logging.info("Calculating the cutoff of motif score...")
        cnt = 0
        cutoff = sampling_times
        while cutoff > 1:
            cnt += 1
            cutoff = int(cutoff * 0.1)
            motifs["score_cutoff_{}".format(cnt)] = np.array(motif_score_all.iloc[cutoff - 1])
        return motifs
    except:
        logging.error("Sampling failed!")
        raise

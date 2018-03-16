"""Module script for motif IO."""

import re
import logging
import numpy as np
import pandas as pd


class InvalidMotifMatrixFormat(Exception):
    pass


class InvalidMotifMatrixValue(Exception):
    pass


def motif_matrix_normalization(matrix):
    """Normalize the motif PFW and compute the normalized PWM.

    Args:
        matrix (np.ndarray): Motif PFM.

    Returns:
        normalized_matrix (np.ndarray): Normalized motif PWM.

    """
    normalized_matrix = np.zeros_like(matrix)
    motif_length = matrix.shape[1]
    for i in xrange(motif_length):
        nt_count = matrix[:, i]
        zero_pos = nt_count < 0.001
        zero_num = zero_pos.sum()
        if zero_num > 3:
            raise InvalidMotifMatrixValue("Invalid PFM, the elements in a column are all zero!")
        if zero_num > 0:
            nt_count[zero_pos] = nt_count.sum() / (1000 - zero_num)
        nt_count = nt_count / nt_count.sum()
        normalized_matrix[:, i] = nt_count
    return normalized_matrix


def load_motif_pfm(motif_pfm_file):
    """Load motif PFM raw matrix.

    Args: 
        motif_pfm_file (str): Path of motif PFM input file.

    Returns:
        motifs (pd.DataFrame): Motif dataframe.

    """
    motif_id = []
    motif_name = []
    motif_matrix = []
    with open(motif_pfm_file, 'r') as fin:
        tmp_matrix = []
        for line in fin:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':
                # logging.info("Loading motif: {}".format(line))
                tmp_line = re.findall(r"[\w:./()-]+", line)
                if len(tmp_line) < 2:
                    raise InvalidMotifMatrixFormat("Invalid Motif header: {}".format(line))
                motif_id.append(tmp_line[0])
                tmp_name = ' '.join(tmp_line[1:])
                motif_name.append(tmp_name)
                if tmp_matrix:
                    # normalize the matrix
                    tmp_matrix = np.asarray(tmp_matrix, dtype=float)
                    tmp_matrix = motif_matrix_normalization(matrix=tmp_matrix)
                    motif_matrix.append(tmp_matrix)
                    tmp_matrix = []
            else:
                try:
                    tmp_matrix.append([float(i) for i in line.split()])
                except ValueError:
                    tmp_matrix.append([float(i) for i in line.strip("ACGT[] ").split()])
        if tmp_matrix:  # append the last matrix
            tmp_matrix = np.asarray(tmp_matrix, dtype=float)
            tmp_matrix = motif_matrix_normalization(matrix=tmp_matrix)
            motif_matrix.append(tmp_matrix)
        if len(set(motif_id)) != len(motif_id):
            raise InvalidMotifMatrixFormat("Motif ID should be unique!")
    motifs = pd.DataFrame({'id': motif_id, 'name': motif_name, 'matrix': motif_matrix}, index=motif_id)
    return motifs


def write_motif_table(motifs, sampling_times, output_prefix):
    """Write the compiled motif table to output file.

    Args:
        motifs (pd.DataFrame): Motif DataFrame.
        sampling_times (int): Sampling times.
        output_prefix (str): Prefix of output file.  

    """
    digit_num = len(str(sampling_times))
    if digit_num > 7:
        digit_num = 7
    for i in xrange(2, digit_num):
        fout = open("{0}_1e-{1}.txt".format(output_prefix, i), 'w')
        for j, motif in motifs.iterrows():
            fout.write(">{0}\t{1}\tmax_score:{2}\tscore_cutoff:{3}\n".format(motif['id'], motif['name'],
                                                                             motif['max_score'],
                                                                             motif['score_cutoff_{}'.format(i)]))
            for k in xrange(4):
                fout.write("{}\n".format('\t'.join(map(str, motif['matrix'][k]))))
        fout.close()
    return


def read_motif_table(motif_file):
    """Read the motif table from input file to a pandas dataframe.

    Args:
        motif_file (str): Path of motif table input file.

    Returns:
        motifs (pd.DataFrame): Motif dataframe.

    """
    fin = open(motif_file, 'r')
    motif_id = []
    motif_name = []
    max_score = []
    matrix = []
    score_cutoff = []
    cnt = 1
    for line in fin:
        line = line.strip()
        if cnt % 5 == 1:  # motif header line
            tmp_line = line.split('\t')
            motif_id.append(tmp_line[0][1:])
            motif_name.append(tmp_line[1])
            max_score.append(float(tmp_line[2].split(':')[1]))
            score_cutoff.append(float(tmp_line[3].split(':')[1]))
        else:
            tmp_line = line.split('\t')
            tmp_array = np.array(map(float, tmp_line))
            if cnt % 5 == 2:
                a = tmp_array
            elif cnt % 5 == 3:
                c = tmp_array
            elif cnt % 5 == 4:
                g = tmp_array
            elif cnt % 5 == 0:
                t = tmp_array
                matrix.append(np.array([a, c, g, t]))
        cnt += 1
    fin.close()
    motifs = pd.DataFrame({'id': motif_id, 'name': motif_name, 'max_score': max_score,
                           'score_cutoff': score_cutoff, 'matrix': matrix})
    motifs.set_index(motifs['id'], inplace=True)
    return motifs


def load_motif(motif_file, filter_file):
    """Read the motif table from input file to a pandas dataframe, and then filter out motifs not in the given list.

     Args:
        motif_file (str): Path of motif table input file.
        filter_file (str): Path of retained motifs input file.

    Returns:
        motifs (pd.DataFrame): Motif dataframe after filtering (if specified).

    """
    motifs = read_motif_table(motif_file=motif_file)
    if filter_file is not None:
        motif_filter = pd.read_csv(filter_file, header=None, names=['motif'])['motif'].tolist()
        motif_names = motifs['name'].tolist()
        for filter_name in motif_filter:
            if filter_name not in motif_names:
                logging.warn("Pass: motif {} in the filter list is not in whole motif list!".format(filter_name))
        motifs = motifs[motifs['name'].isin(motif_filter)].copy()

    dup_name_idx = []
    for idx, motif in motifs.iterrows():
        if np.sum(motifs['name'] == motif['name']) > 1:
            dup_name_idx.append(idx)
    for idx in dup_name_idx:
        motif = motifs.loc[idx]
        motifs.loc[idx, 'name'] = "{0}_{1}".format(motif['name'], motif['id'])
    return motifs

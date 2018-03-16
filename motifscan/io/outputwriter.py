"""Module script for output."""

import os
import re
import logging
import numpy as np
import pandas as pd
from collections import namedtuple
from motifscan.lib.genome import extract_sequence


def make_output_dir(output_dir, target_site_flag):
    """Make the output directories for to write the results.

     Args:
         output_dir: Path of output directory.
         target_site_flag (bool): Whether make the target_site_dir or not.

    Returns:
        output_stream(dict): A namedtuple which contains the sub-dirs and output file path.

    """
    enrichment_file = output_dir + '/motif_enrichment.csv'
    target_number_file = output_dir + '/peak_motif_target_number.csv'
    motif_score_file = output_dir + '/peak_motif_score.csv'
    plot_dir = output_dir + '/' + 'plot'
    target_site_dir = output_dir + '/' + 'motif_target_sites'
    tmp_dir = output_dir + '/.tmp'

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    if target_site_flag:
        if not os.path.exists(target_site_dir):
            os.mkdir(target_site_dir)
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    OutputFiles = namedtuple('OutputStream', ['enrichment_file', 'target_number_file', 'motif_score_file', 'plot_dir',
                                              'target_site_dir', 'tmp_dir'])
    output_stream = OutputFiles(enrichment_file, target_number_file, motif_score_file, plot_dir, target_site_dir,
                                tmp_dir)
    return output_stream


def output_target_number(output_file, data, columns):
    """Write the motif target number of given regions to output file.
    """
    data.to_csv(output_file, index=False, header=True, columns=columns)
    return


def output_motif_score(output_file, data, columns):
    """Write the motif score of given regions to output file.
    """
    data.to_csv(output_file, index=False, header=True, columns=columns, float_format="%.2f")
    return


def output_target_site(output_dir, peaks, motifs, genome_dir, tmp_dir):
    """Write the motif target site in given regions to output file.
    """
    for idx, motif_record in motifs.iterrows():
        name = motif_record['name']
        motif_len = np.shape([motif_record['matrix']])[2]
        tmp_df = pd.read_pickle("{0}/{1}".format(tmp_dir, motif_record['id']))
        # output details
        output_file = "{0}/{1}_target_site.bed".format(output_dir, re.sub(r'[:\-.]{1,2}', '_', name))
        fout = open(output_file, 'w')
        for peak_idx, target_site in tmp_df['%s.tarsite' % name].iteritems():
            if len(target_site) > 0:
                for i, start in enumerate(target_site):
                    tar_chr = peaks.iloc[peak_idx]['chr']
                    tar_start = int(peaks.iloc[peak_idx]['seq_start'] + start)
                    tar_end = int(tar_start + motif_len)
                    tar_seq = extract_sequence(genome_dir, tar_chr, tar_start, tar_end)
                    tar_ratio = tmp_df["{}.tarratio".format(name)].iloc[peak_idx][i]
                    fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(tar_chr, tar_start, tar_end, tar_seq, tar_ratio))
        fout.close()
    return


def output_enrichment_result(output_file, data):
    """Write the enrichment analysis result.
    """
    data.sort_values(by='enriched_pvalue', inplace=True)
    data.to_csv(output_file, index=False, columns=['name', 'target_number', 'random_target_number', 'fold_change',
                                                   'enriched_pvalue', 'depleted_pvalue', 'pvalue_corrected'])
    return

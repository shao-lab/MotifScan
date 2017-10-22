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

"""Module script for MotifScan backbone."""

import os
import sys
import shutil
import logging
import pandas as pd

import motifscan.core as core
import motifscan.enrich as enrich
import motifscan.plot as plot
import motifscan.io.geneio as geneio
import motifscan.io.peakio as peakio
import motifscan.io.motifio as motifio
import motifscan.io.genomeio as genomeio
from motifscan.io.outputwriter import *
import motifscan.lib.peak as peak
from motifscan.lib.genome import extract_sequence
from motifscan.lib.matrix import construct_sequence_matrix
from motifscan.utils import summarize_result


def motifscan_run(genome_dir, gene_file, motif_file, motif_filter_file, peak_file, peak_format, peak_length,
                  control_file, control_format, location, upstream, downstream, random_times, enrichment_flag,
                  target_site_flag, output_dir):
    output_stream = make_output_dir(output_dir=output_dir, target_site_flag=target_site_flag)
    logging.info("### Loading data ###")
    logging.info("Loading genome...")
    bg_freq, chr_size = genomeio.load_genome(genome_dir=genome_dir)
    if gene_file:
        logging.info("Loading gene annotation...")
        genes = geneio.load_refseq_gene(input_file=gene_file)
    else:
        genes = None
        if location != 'all':
            logging.error("Gene annotation file is missing when MotifScan is performed on {}!".format(location))
            sys.exit(1)

    logging.info("Loading motifs...")
    motifs = motifio.load_motif(motif_file=motif_file, filter_file=motif_filter_file)
    if len(motifs) == 0:
        logging.error("There is no motif in the final motif table!")
        sys.exit(1)
    else:
        logging.info("{} motifs are loaded!".format(len(motifs)))

    logging.info("Loading peaks...")
    peaks = peakio.load_peak(peak_file=peak_file, peak_format=peak_format, peak_length=peak_length)
    logging.info("{} peaks are loaded!".format(len(peaks)))
    if location in ['promoter', 'distal']:
        logging.info("Extracting peaks of {} regions...".format(location))
        peaks = peak.split_peaks_by_location(peaks=peaks, genes=genes, location=location, upstream=upstream,
                                             downstream=downstream)
        logging.info("{0} peaks are extracted!".format(len(peaks)))
    elif location == 'gene':
        logging.info("Extracting gene targeted peaks...".format(location))
        logging.info("Find target gene of peaks...")
        peaks = peak.get_target_gene(peaks=peaks, genes=genes)
        peaks = peaks.loc[peaks['target_gene'].notnull()]
        # pick up the peak that is nearest to the target gene
        targeted_peak_idx = peaks.groupby(['target_gene'])['target_dis'].transform(lambda x: min(abs(x))) == peaks[
            'target_dis'].abs()
        peaks = peaks.loc[targeted_peak_idx].copy()
        peaks.reset_index(inplace=True, drop=True)
        logging.info("{0} peaks are extracted!".format(len(peaks)))

    logging.info("### Scanning motifs on input regions ###")
    logging.info("Extracting sequences...")
    genome_dir_iter = [genome_dir] * len(peaks)
    peaks['seq'] = map(extract_sequence, genome_dir_iter, peaks['chr'], peaks['seq_start'], peaks['seq_end'])
    peaks['seq_matrix'] = map(construct_sequence_matrix, peaks['seq'])

    # Scanning motifs...
    core.scan_motif(peaks=peaks, motifs=motifs, bg_freq=bg_freq, output_dir=output_stream.tmp_dir)
    peak_result, target_number_cols, motif_score_cols = summarize_result(peaks=peaks, motifs=motifs,
                                                                         output_dir=output_stream.tmp_dir)
    logging.info("Saving the target number and motif score...")
    output_target_number(output_file=output_stream.target_number_file, data=peak_result, columns=target_number_cols)
    output_motif_score(output_file=output_stream.motif_score_file, data=peak_result, columns=motif_score_cols)
    if target_site_flag:
        logging.info("Saving the target sites of motifs...")
        output_target_site(output_dir=output_stream.target_site_dir, peaks=peak_result, motifs=motifs,
                           genome_dir=genome_dir, tmp_dir=output_stream.tmp_dir)

    if enrichment_flag:
        if control_file:
            logging.info("### Scanning motifs on control regions ###")
            random_regions = peakio.load_peak(peak_file=control_file, peak_format=control_format,
                                              peak_length=peak_length)
            logging.info("{} control regions are loaded!".format(len(random_regions)))
        else:
            logging.info("### Scanning motifs on random regions ###")
            logging.info("Generating random control regions based on {0} {1} peaks...".format(len(peaks), location))
            if genes is not None:
                random_regions = peak.generate_random_with_ref(peaks=peaks, genes=genes, random_times=random_times)
            else:
                random_regions = peak.generate_random_without_ref(peaks=peaks, chr_size=chr_size,
                                                                  random_times=random_times)
            logging.info("{} random control regions are generated!".format(len(random_regions)))
        logging.info("Extracting sequences...")
        genome_dir_iter = [genome_dir] * len(random_regions)
        if control_file:
            random_regions['seq'] = map(extract_sequence, genome_dir_iter, random_regions['chr'],
                                        random_regions['seq_start'], random_regions['seq_end'])
        else:
            random_regions['seq'] = map(extract_sequence, genome_dir_iter, random_regions['chr'],
                                        random_regions['start'], random_regions['end'])
        random_regions['seq_matrix'] = map(construct_sequence_matrix, random_regions['seq'])
        core.scan_motif(peaks=random_regions, motifs=motifs, bg_freq=bg_freq, output_dir=output_stream.tmp_dir)
        random_result = {}
        for idx, tmp_motif in motifs.iterrows():
            name = tmp_motif['name']
            tmp_df = pd.read_pickle("{0}/{1}".format(output_stream.tmp_dir, idx))
            random_result["{}.number".format(name)] = tmp_df["{}.number".format(name)]
        random_result = pd.DataFrame(random_result)

        logging.info("### Performing enrichment analysis ###")
        enrich_result = enrich.target_enrichment(peak_table=peak_result, rnd_table=random_result, motif_table=motifs)
    if enrichment_flag:
        logging.info("Saving the enrichment analysis results...")
        output_enrichment_result(output_file=output_stream.enrichment_file, data=enrich_result)
    if False:
        logging.info("Plotting the target site distribution of motifs...")
        if (not enrichment_flag or 'value' not in peaks.columns) and peak_length != 0:
            plot.target_site_distribution(peak_df=peak_result, motif_df=motifs, plot_dir=output_stream.plot_dir,
                                          region_radius=peak_length / 2)
        if enrichment_flag and 'value' in peaks.columns and peak_length != 0:
            plot.tarnum_and_tarsite_distribution(peak_table=peak_result, rand_table=random_result, motif_table=enrich_result, plot_out_dir=output_stream.plot_dir,
                                                 region_radius=peak_length / 2)
    logging.info("### MotifScan Finished ###")
    shutil.rmtree(output_stream.tmp_dir)
    if enrichment_flag:
        return peak_result, random_result, enrich_result
    else:
        return peak_result, None, None

# def motifcompare_general(result_a, result_b, output_prefix):
#     rst_a_df = pd.read_csv("{}/motif_enrichment.csv".format(result_a), index_col=0)
#     rst_b_df = pd.read_csv("{}/motif_enrichment.csv".format(result_b), index_col=0)
#     enrich = core.target_enrichment_peak2peak(rst_a_df, rst_b_df)
#     enrich.to_csv(output_prefix)
#
#
# def merge_two_results(peak_result_1, peak_result_2, rnd_result_1, rnd_result_2, motif_table):
#     merged_peak_result = pd.concat([peak_result_1, peak_result_2], ignore_index=True)
#     merged_rnd_result = pd.concat([rnd_result_1, rnd_result_2], ignore_index=True)
#     merged_enrich_result = core.target_enrichment(merged_peak_result, merged_rnd_result, motif_table)
#     return [merged_peak_result, merged_rnd_result, merged_enrich_result]
#
#
# def merge_two_results_no_rnd(peak_result_1, peak_result_2):
#     merged_peak_result = pd.concat([peak_result_1, peak_result_2], ignore_index=True)
#     return merged_peak_result

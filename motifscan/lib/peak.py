"""Module script for peak(genomic regions)-related operations."""

import os
import re
import sys
import random
from collections import defaultdict

import numpy as np
import pandas as pd


def merge_overlap(regions):
    """Merge the overlap regions.
    
    Args:
        regions (dict): Regions to be merged.
    
    Returns:
        merged_regions (dict): Merged regions.
        
    """
    merged_regions = {}
    for tmp_chr in regions.keys():
        tmp_regions = regions[tmp_chr]
        tmp_regions.sort()
        tmp_merged_regions = []
        pre_region = None
        for idx, tmp_region in enumerate(tmp_regions):
            if not pre_region:
                pre_region = tmp_region
                continue
            else:
                if pre_region[1] >= tmp_region[0]:
                    pre_region = (pre_region[0], max(pre_region[1], tmp_region[1]))
                else:
                    tmp_merged_regions.append(pre_region)
                    pre_region = tmp_region
        if pre_region:
            tmp_merged_regions.append(pre_region)
        if tmp_merged_regions:
            merged_regions[tmp_chr] = tmp_merged_regions
        #print tmp_chr, tmp_merged_regions[:10]
    return merged_regions


def overlap_with(regions, chrom, start, end):
    """Returns whether the regions in chrom overlap with (start, end).
    
    Args:
        regions (dict): Reference regions.
        chrom (str): Chromosome name.
        start (int): Start pos.
        end (int): End pos.
    
    Returns:
        True/False: Overlap or not.
        
    """
    if chrom not in regions:
        return False
    ref_regions = regions[chrom]
    left = 0
    right = len(ref_regions) - 1
    while left <= right:
        mid = (left + right) / 2
        head = ref_regions[mid][0]
        tail = ref_regions[mid][1]
        if not (end <= head or start >= tail):
            return True
        elif start >= tail:
            left = mid + 1
        elif end <= head:
            right = mid - 1
    return False


def split_peaks_by_location(peaks, genes, location, upstream, downstream):
    """Extract peaks located in promoter/distal regions.
    
    Args:
        peaks (pd.DataFrame): Peak dataframe.
        genes (pd.DataFrame): Gene dataframe.
        location (str): Location ('promoter' or 'distal') to extract peaks.
        upstream (int): Upstream from gene TSS to define promoter regions.
        downstream (int): Downstream from gene TSS to define promoter regions.

    Returns:
        peaks (pd.DataFrame): Peak dataframe with peaks located at promoter or distal.
    
    """
    promoters = defaultdict(list)
    for idx, tmp_gene in genes.iterrows():
        if tmp_gene['strand'] == '+':
            promoters[tmp_gene['chr']].append((tmp_gene['TSS'] - upstream, tmp_gene['TSS'] + downstream))
        else:
            promoters[tmp_gene['chr']].append((tmp_gene['TSS'] - downstream, tmp_gene['TSS'] + upstream))
    promoters = merge_overlap(promoters)
    loc = []
    for idx, tmp_peak in peaks.iterrows():
        if overlap_with(promoters, tmp_peak['chr'], tmp_peak['start'], tmp_peak['end']):
            loc.append('promoter')
        else:
            loc.append('distal')
    peaks['location'] = loc
    peaks = peaks.loc[peaks['location'] == location].copy()
    peaks.reset_index(inplace=True, drop=True)
    return peaks


def get_target_gene(peaks, genes, distance_cutoff=10000):
    """Find the target gene and distance of peaks.
    
    Args:
        peaks (pd.DataFrame): Peak dataframe.
        genes (pd.DataFrame): Gene dataframe.
        distance_cutoff (int): Distance cutoff of peak to nearest gene.
    
    Returns:
        peaks (pd.DataFrame): Peak dataframe with target gene information.
     
    """
    nearest_gene = []
    nearest_dis = []
    for idx, tmp_peak in peaks.iterrows():
        tmp_chr = tmp_peak['chr']
        tmp_genes = genes.loc[genes['chr'] == tmp_chr]
        tmp_target = np.nan
        tmp_dis = np.nan
        target_dis = tmp_peak['summit'] - tmp_genes['TSS']
        abs_target_dis = target_dis.abs()
        if abs_target_dis.min() <= distance_cutoff:
            idx_min = abs_target_dis.argmin()
            tmp_target = tmp_genes.loc[idx_min, 'id']
            tmp_dis = target_dis[idx_min]
        nearest_gene.append(tmp_target)
        nearest_dis.append(tmp_dis)
    peaks['target_gene'] = nearest_gene
    peaks['target_dis'] = nearest_dis
    return peaks


def generate_random_without_ref(peaks, chr_size, random_times):
    # generate the random control
    chrs = []
    starts = []
    ends = []
    seq_start_col_iloc = peaks.columns.get_loc('seq_start')
    seq_end_col_iloc = peaks.columns.get_loc('seq_end')
    peak_num = len(peaks)
    for i in xrange(peak_num):
        tmp_length = peaks.iat[i, seq_end_col_iloc] - peaks.iat[i, seq_start_col_iloc]
        for j in xrange(random_times):
            tmp_chr = peaks['chr'].iloc[i]
            tmp_start = random.randint(0, chr_size[tmp_chr] - tmp_length)
            tmp_end = tmp_start + tmp_length
            chrs.append(tmp_chr)
            starts.append(tmp_start)
            ends.append(tmp_end)
    random_regions = pd.DataFrame({'chr': chrs, 'start': starts, 'end': ends})
    return random_regions


def generate_random_with_ref(peaks, genes, random_times):
    """generate the random sequence according the peak summit distance relative to the target gene TSS
    """
    merged_df = merge_peak_gene(peaks, genes)
    chr_size = chr_size_margin(peaks, genes)
    dis = []
    chrs = merged_df['chr'].unique()
    for tmp_chr in chrs:
        dis += dis_to_target_gene(merged_df.loc[merged_df['chr'] == tmp_chr].copy())
    random_regions = peak_random(peaks, genes, dis, chr_size, random_times)
    return random_regions


def merge_peak_gene(peak_df, gene_df):
    """merge peak table and gene table and sorting by their position


    Args:
        peak_df: pandas dataframe
        gene_df: pandas dataframe
    Returns:
        merged table: pandas dataframe with 4 columns: chr, position, strand, type
        Note: type "gene" is for record from gene table, type "peak" is for record from peak table.
              positon for gene is TSS; position for peak is peak summit.
              the index of dataframe are reset by the sorting.

    """
    peak_part = peak_df[['chr', 'summit']].copy()
    peak_part['strand'] = 'unknown'
    peak_part['type'] = 'peak'
    peak_part.columns = ['chr', 'position', 'strand', 'type']

    gene_part = gene_df[['chr', 'TSS', 'strand']].copy()
    gene_part['type'] = 'gene'
    gene_part.columns = ['chr', 'position', 'strand', 'type']

    merged_df = peak_part.append(gene_part, ignore_index=True)
    merged_df.sort_values(by=['chr', 'position'], inplace=True)
    merged_df.reset_index(inplace=True, drop=True)
    return merged_df


def chr_size_margin(peak_table, gene_table):
    chr_size_dict = {}
    chrs = gene_table['chr'].unique()

    for chr_id in chrs:
        if chr_id in peak_table['chr']:
            chr_size_dict[chr_id] = max(peak_table[peak_table['chr'] == chr_id]['end'].max(),
                                        gene_table[gene_table['chr'] == chr_id]['end'].max())
        else:
            chr_size_dict[chr_id] = gene_table[gene_table['chr'] == chr_id]['end'].max()
    return chr_size_dict


def dis_to_target_gene(df, distance_cutoff_abs=10000):
    """Compute each peak's distance to its target gene
        Append the column named after distance_to_target_gene in the peak_gene indicating the distance.

    Args:
        df: dataframe  by merge_peak_gene's output
        distance_cutoff_abs: when distance to target gene > 10k, the target gene is considered not exist.

    Returns:

    """
    df['distance_to_target_gene'] = -1
    df.reset_index(inplace=True, drop=True)
    dis = []
    peak_idx = df.loc[df['type'] == 'peak'].index
    for i in peak_idx:
        tmp_peak = df.iloc[i]
        # search for the nearest gene in upstream
        dis_up = float('inf')
        if i > 0:
            for j in xrange(i - 1, -1, -1):
                tmp_item = df.iloc[j]
                if tmp_item['type'] == 'peak':  # skip peak record
                    continue
                dis_up = tmp_peak['position'] - tmp_item['position']
                if tmp_item['strand'] == '-':
                    dis_up = - dis_up
                break
        # search for the nearest gene in downstream
        dis_down = float('-inf')
        if i < len(df):
            for j in np.arange(i + 1, len(df), 1):
                tmp_item = df.iloc[j]
                if tmp_item['type'] == 'peak':  # skip peak record
                    continue
                dis_down = tmp_peak['position'] - tmp_item['position']
                if tmp_item['strand'] == '-':
                    dis_down = -dis_down
                break
        if abs(dis_up) < abs(dis_down):
            tmp_dis = dis_up
        else:
            tmp_dis = dis_down
        dis.append(tmp_dis)

    for idx, tmp_dis in enumerate(dis):
        if tmp_dis == float('inf'):
            dis[idx] = np.random.randint(distance_cutoff_abs, 10 * distance_cutoff_abs)
        elif tmp_dis == float('-inf'):
            dis[idx] = -np.random.randint(distance_cutoff_abs, 10 * distance_cutoff_abs)
        else:
            if abs(tmp_dis) >= distance_cutoff_abs:
                dis[idx] = (tmp_dis / abs(tmp_dis)) * np.random.randint(distance_cutoff_abs, 10 * distance_cutoff_abs)
            else:
                dis[idx] = tmp_dis
    return dis


def peak_random(peak_df, gene_df, dis, chr_size, rnd_times=5):
    """generate rnd peak table

    Args:
        peak_df
        dis: contain all distances
        gene_df: pandas dataframe, gene table
        chr_size: chromosome margin of each chromosome
        rnd_times: random times for each distance

    Returns:
        pandas dataframe containing random times

    """

    gene_part = gene_df[['chr', 'TSS', 'strand']].copy()
    gene_part['type'] = 'gene'
    gene_part.columns = ['chr', 'position', 'strand', 'type']

    rnd_chr = []
    rnd_start = []
    rnd_end = []
    peak_num = len(peak_df)
    for i in xrange(peak_num):
        tmp_peak = peak_df.iloc[i]
        tmp_length = int(tmp_peak['seq_end'] - tmp_peak['seq_start'])
        tmp_dis = dis[i]
        rnd_idx = random.sample(np.arange(len(gene_part)), rnd_times)  # for each peak, generate 5 samples randomly
        for j in rnd_idx:
            tmp_gene = gene_part.iloc[j]
            tmp_chr = tmp_gene['chr']
            tmp_size = chr_size[tmp_chr]
            rnd_chr.append(tmp_gene['chr'])
            if tmp_gene['strand'] == '+':
                tmp_summit = tmp_gene['position'] + tmp_dis
            else:
                tmp_summit = tmp_gene['position'] - tmp_dis
            if tmp_summit - tmp_length / 2 < 0:
                tmp_summit = tmp_length / 2
            if tmp_summit - tmp_length / 2 + tmp_length >= tmp_size:
                tmp_summit = tmp_size + tmp_length / 2 - tmp_length - 1
            tmp_start = tmp_summit - tmp_length / 2
            tmp_end = tmp_start + tmp_length
            # if tmp_end > tmp_size:
            #     tmp_end = tmp_size
            rnd_start.append(tmp_start)
            rnd_end.append(tmp_end)
    rnd_df = pd.DataFrame({'chr': rnd_chr, 'start': rnd_start, 'end': rnd_end})
    return rnd_df

#
# Copyright @ 2014, 2015 Jiawei Wang <jerryeah@gmail.com>, Zhen Shao <shao@enders.tch.harvard.edu>
#
# Licensed under the GPL License, Version 3.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#   http://www.gnu.org/copyleft/gpl.html
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 

"""The core functions of MotifScan.
    
    motif_scan(): the interface of MotifScan core
    target_enrichment(): perform the enrichment analysis between sample and random
    target_enrichment_peak2peak(): perform enrichment analysis between two samples
    fc_tarnum_distribution(): fold change plot graph of each motif ranked by peak mvalue
    target_site_distribution(): motif target distribution graph around peak summit

"""

import os
import re
import sys
import ctypes
import logging
import numpy as np
import pandas as pd
from pkg_resources import resource_filename
import motifscan
from motifscan.lib.matrix import MAT, arr_to_mat


class MotifRes(ctypes.Structure):
    _fields_ = [("number", ctypes.c_int),
                ("ratio", ctypes.c_double),
                ("tarsite", ctypes.POINTER(ctypes.c_int)),
                ("tarratio", ctypes.POINTER(ctypes.c_double))]


def deduplicate_target_site(ts, tr, motif_len):
    """
        if the distance between two neignbor target site is smaller than
        the motif length, we viewed them as duplicate target sites and thus filter out
        the one with the lower ratio.
    """
    if len(ts) >= 2:
        i_pre = float("-Inf")
        j_pre = float("-Inf")
        for (i, j) in zip(ts, tr):
            if i - i_pre < motif_len:
                if j < j_pre:
                    ts.remove(i)
                    tr.remove(j)
                else:
                    ts.remove(i_pre)
                    tr.remove(j_pre)
                    i_pre = i
                    j_pre = j
            else:
                i_pre = i
                j_pre = j
    return ts, tr


def scan_motif(peaks, motifs, bg_freq, output_dir):
    """ The core interface of motif scanning.

    Args:
        peaks (pd.DataFrame): Peak dataframe.
        motifs (pd.DataFrame): Motif dataframe.
        bg_freq (pd.Series): A pd.Series with 4 float values representing the whole genome frequency of nucleotides.
        output_dir (str): temporary motifscan output file for a certain motif

    """
    score_c = ctypes.CDLL(resource_filename(motifscan.__name__, 'score_c.so'))
    score_c.motif_scan_core.restype = ctypes.POINTER(MotifRes)
    score_c.motif_scan_core.argtypes = [ctypes.POINTER(MAT), ctypes.POINTER(MAT), ctypes.POINTER(ctypes.c_double * 4),
                                        ctypes.c_double, ctypes.c_double]
    score_c.freeMOTIF_RES.argtypes = [ctypes.POINTER(MotifRes)]

    motif_num = len(motifs)
    bg_freq = np.asarray(bg_freq)
    bg_freq = np.ctypeslib.as_ctypes(bg_freq)

    peak_num = len(peaks)
    seq_mat_col_iloc = peaks.columns.get_loc('seq_matrix')
    cnt = 0
    for idx, tmp_motif in motifs.iterrows():
        cnt += 1
        logging.info("Scanning motif: {0}\t{1}/{2}...".format(tmp_motif['name'], cnt, motif_num))
        motif_name = tmp_motif['name']
        motif_len = tmp_motif['matrix'].shape[1]
        score_cutoff = ctypes.c_double(tmp_motif['score_cutoff'])
        max_score = ctypes.c_double(tmp_motif['max_score'])
        m_matrix = arr_to_mat(tmp_motif['matrix'])

        ratio, number, tarsite, tarratio = [], [], [], []
        for p_idx in xrange(peak_num):
            s_matrix = arr_to_mat(peaks.iat[p_idx, seq_mat_col_iloc])
            tmp_rst = score_c.motif_scan_core(ctypes.byref(s_matrix), ctypes.byref(m_matrix), ctypes.byref(bg_freq),
                                              max_score, score_cutoff)

            ratio.append(tmp_rst.contents.ratio)
            ctypes_tarsite = (ctypes.c_int * tmp_rst.contents.number).from_address(
                ctypes.addressof(tmp_rst.contents.tarsite.contents))
            tmp_tar_site = []
            for i in np.arange(tmp_rst.contents.number):
                tmp_tar_site.append(ctypes_tarsite[i])

            ctypes_tarratio = (ctypes.c_double * tmp_rst.contents.number).from_address(
                ctypes.addressof(tmp_rst.contents.tarratio.contents))
            tmp_tar_ratio = []
            for i in np.arange(tmp_rst.contents.number):
                tmp_tar_ratio.append(ctypes_tarratio[i])

            # deduplicate overlap motif target sites
            tmp_tar_site, tmp_tar_ratio = deduplicate_target_site(tmp_tar_site, tmp_tar_ratio, motif_len)
            number.append(len(tmp_tar_ratio))
            tarsite.append(tmp_tar_site)
            tarratio.append(tmp_tar_ratio)

            score_c.freeMOTIF_RES(tmp_rst)
        peak_result = pd.DataFrame({"{}.ratio".format(motif_name): ratio, "{}.number".format(motif_name): number,
                                    "{}.tarsite".format(motif_name): tarsite,
                                    "{}.tarratio".format(motif_name): tarratio})
        peak_result.to_pickle("{0}/{1}".format(output_dir, tmp_motif['id']))
    return

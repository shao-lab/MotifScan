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

"""Module script for sequence matrix operations."""

import ctypes
import numpy as np


class MAT(ctypes.Structure):
    _fields_ = [("n", ctypes.c_int),
                ("a_arr", ctypes.POINTER(ctypes.c_double)),
                ("c_arr", ctypes.POINTER(ctypes.c_double)),
                ("g_arr", ctypes.POINTER(ctypes.c_double)),
                ("t_arr", ctypes.POINTER(ctypes.c_double))]


def arr_to_mat(arr):
    n = arr.shape[1]
    a = np.ctypeslib.as_ctypes(arr[0])
    c = np.ctypeslib.as_ctypes(arr[1])
    g = np.ctypeslib.as_ctypes(arr[2])
    t = np.ctypeslib.as_ctypes(arr[3])
    return MAT(n, a, c, g, t)


def construct_sequence_matrix(seq):
    """Construct sequence matrix based on given sequences.

    Args:
        seq (str): Sequences of nucleotides.

    Returns:
        matrix (np.array): Sequence matrix.

    """
    matrix = np.zeros((4, len(seq)))
    base_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for base, idx in zip(seq, xrange(len(seq))):
        try:
            matrix[base_idx[base], idx] = 1
        except KeyError:
            continue
    return matrix

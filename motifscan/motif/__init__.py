"""
motifscan.motif
---------------

Module for motif related classes and functions.
"""

import logging
import os
import re

from motifscan.config import Config
from motifscan.exceptions import PfmsFileNotFoundError, \
    PwmsFileNotFoundError, PfmsJasparFormatError, PwmsMotifScanFormatError
from motifscan.genome import bases
from motifscan.motif.matrix import PositionFrequencyMatrix, \
    PositionWeightMatrix

logger = logging.getLogger(__name__)

pfms_path_fmt = os.path.join("{0}", "{1}_pfms.jaspar")
pwms_path_fmt = os.path.join("{0}", "{1}_{2}_pwms.motifscan")


class MotifMatrices:
    """Generic class for motif matrices."""

    def __init__(self):
        self._matrices = []

    def __iter__(self):
        yield from self._matrices

    def __len__(self):
        return len(self._matrices)

    def append(self, item):
        self._matrices.append(item)

    def extend(self, items):
        self._matrices.extend(items)


class MotifPfms(MotifMatrices):
    """Class for a set of motif PFMs (Position Frequency Matrices).

    Parameters
    ----------
    pfms : list of `PositionFrequencyMatrix`, optional
        PFMs used to construct the class.
    name : str, optional
        The name of the motif PFMs set.

    Attributes
    ----------
    name : str or None
        The name of the motif PFMs set, or None if not specified.
    """

    def __init__(self, pfms=None, name=None):
        super().__init__()
        self.name = name
        if pfms is not None:
            for pfm in list(pfms):
                if isinstance(pfm, PositionFrequencyMatrix):
                    self.append(pfm)
                else:
                    raise ValueError(f"invalid PFM item: {pfm!r}")

    @staticmethod
    def _parse_jaspar_pfms(path):
        """Parse motif PFMs in JASPAR format.

        The regular expression pattern of the header/matrix line is inspired by
        `Biopython`.

        JASPAR PFM Example:

            >MA0006.1	Ahr::Arnt
            A  [     3      0      0      0      0      0 ]
            C  [     8      0     23      0      0      0 ]
            G  [     2     23      0     23      0     24 ]
            T  [    11      1      1      1     24      0 ]

        Raises
        ------
        PfmsJasparFormatError
            If the file does not strictly follow the JASPAR PFMs format.
        """
        header_pattern = re.compile(r"^>\s*(\S+)(\s+(\S+))?")
        matrix_new_pattern = re.compile(r"\s*([ACGT])\s*\[\s*(.+)\s*\]")
        matrix_old_pattern = re.compile(r"\s*(.+)\s*")

        pfms = []
        line_num = 0
        expect_header = True  # expect header (True) or matrix line (False)
        with open(path, 'r') as fin:
            for line in fin:
                line_num += 1
                line = line.strip()
                if not line:  # skip blank lines
                    continue

                m_header = header_pattern.match(line)
                m_matrix_new = matrix_new_pattern.match(line)
                m_matrix_old = matrix_old_pattern.match(line)

                if bool(m_header) != expect_header:  # format check
                    raise PfmsJasparFormatError(line_num, line)

                if m_header:
                    matrix_id = m_header.group(1)
                    name = m_header.group(3)
                    n_matrix = 0
                    values = []
                    expect_header = False
                else:
                    if m_matrix_new:
                        base = m_matrix_new.group(1)
                        if base != bases[n_matrix]:
                            raise PfmsJasparFormatError(line_num, line)
                        tmp_values = m_matrix_new.group(2).split()
                    elif m_matrix_old:
                        tmp_values = m_matrix_old.group(1).split()
                    else:  # neither a header nor a matrix line
                        raise PfmsJasparFormatError(line_num, line)
                    try:
                        values.append(list(map(int, tmp_values)))
                    except (ValueError, TypeError):
                        raise PfmsJasparFormatError(line_num, line)
                    n_matrix += 1
                    if n_matrix == 4:
                        pfm = PositionFrequencyMatrix(values=values, name=name,
                                                      matrix_id=matrix_id)
                        pfms.append(pfm)
                        expect_header = True

            if not expect_header:  # check whether the last matrix is complete
                raise PfmsJasparFormatError(line_num + 1, '')
        return pfms

    def read_pfms(self, path, format='jaspar'):
        """Read motif PFMs.

        Parameters
        ----------
        path : str
            Path to load the PFMs.
        format : {'jaspar'}, optional
            PFMs file format, default='jaspar'.
        """
        if format not in ['jaspar']:
            raise ValueError(f"invalid motif PFMs file format: {format!r}")
        logger.debug(f"Reading motif PFMs from {path} [{format}]")
        pfms = self._parse_jaspar_pfms(path)
        self.extend(pfms)
        logger.debug(f"Found {len(pfms)} motif PFMs")


class MotifPwms(MotifMatrices):
    """Class for a set of motif PWMs (Position Weight Matrices).

    Parameters
    ----------
    pwms : list of `PositionWeightMatrix`, optional
        PWMs used to construct the class.
    name : str, optional
        The name of the motif PWMs set.
    genome : str, optional
        The name of the genome assembly under which these PWMs are built.

    Attributes
    ----------
    name : str or None
        The name of the motif PWMs set, or None if not specified.
    genome : str or None
        The name of the genome assembly under which these PWMs are built, or
        None if not specified.
    """

    def __init__(self, pwms=None, name=None, genome=None):
        super().__init__()
        self.name = name
        self.genome = genome
        if pwms is not None:
            for pwm in list(pwms):
                if isinstance(pwm, PositionWeightMatrix):
                    self.append(pwm)
                else:
                    raise ValueError(f"invalid PWM item: {pwm!r}")

    def save_built_pwms(self):
        """Save built motif PWMs."""
        logger.info(
            f"Saving motif PWMs {self.name!r} under assembly {self.genome!r}")
        motif_dir = Config().get_motif_path(self.name)
        pwms_path = pwms_path_fmt.format(motif_dir, self.name, self.genome)
        self.write_motifscan_pwms(pwms_path)

    def write_motifscan_pwms(self, path):
        """Write motif PWMs in MotifScan format.

        Parameters
        ----------
        path : str
            The file path to write the MotifScan PWMs.
        """
        logger.debug(f"Writing MotifScan PWMs to {path}")
        with open(path, 'w') as f_out:
            for pwm in self:
                f_out.write(f">{pwm.matrix_id}\t{pwm.name}\tPWM\n")
                for idx, base in enumerate(bases):
                    values_str = '\t'.join(
                        map(lambda x: f'{x:8.5f}', pwm.matrix[idx]))
                    f_out.write(f"{base} [{values_str}]\n")
                for p, cutoff in pwm.cutoffs.items():
                    f_out.write(f"Cutoff_p{p}\t{cutoff}\n")

    def read_motifscan_pwms(self, path):
        """Read PWMs in MotifScan format.

         MotifScan PWM Example:

            >MA0006.1   Ahr::Arnt   PWM
            A [-0.85815 -5.68647 -5.68647 -5.68647 -5.68647 -5.68647]
            C [ 0.48657 -5.32257  1.53966 -5.32257 -5.32257 -5.32257]
            G [-0.90016  1.53922 -5.32301  1.53922 -5.32301  1.58174]
            T [ 0.43981 -1.93828 -1.93828 -1.93828  1.21696 -5.68779]
            Cutoff_p1e-3	0.55403
            Cutoff_p1e-4	0.82985
            Cutoff_p1e-5	1.0

        Parameters
        ----------
        path : str
            The file path to read the MotifScan PWMs.

        Raises
        ------
        PwmsMotifScanFormatError
            If the file does not strictly follow the MotifScan PWMs format.
        """
        logger.debug(f"Reading MotifScan PWMs from {path}")
        header_pattern = re.compile(r"^>(\S+)\t(\S+)\tPWM$")
        matrix_pattern = re.compile(r"^([ACGT]) \[(.+)\]$")
        cutoff_pattern = re.compile(r"^Cutoff_p(\S+)\t(\S+)")

        pwms = []
        line_num = 0
        # expect_flag: 0=header, 1=matrix, 2=cutoff, 3=cutoff or header
        # 1 header line + 4 matrix line + at least 1 cutoff line
        expect_flag = 0
        with open(path, 'r') as fin:
            for line in fin:
                line_num += 1
                line = line.strip()
                if not line:  # skip blank lines
                    continue

                m_header = header_pattern.match(line)
                m_matrix = matrix_pattern.match(line)
                m_cutoff = cutoff_pattern.match(line)

                # format checker
                if m_header:
                    if expect_flag != 0 and expect_flag != 3:
                        raise PwmsMotifScanFormatError(line_num, line)
                elif m_matrix:
                    if expect_flag != 1:
                        raise PwmsMotifScanFormatError(line_num, line)
                elif m_cutoff:
                    if expect_flag != 2 and expect_flag != 3:
                        raise PwmsMotifScanFormatError(line_num, line)
                else:  # does not match any pattern
                    raise PwmsMotifScanFormatError(line_num, line)

                if m_header:
                    if expect_flag == 3:  # already got a pwm and save it
                        pwm = PositionWeightMatrix(values=values, name=name,
                                                   matrix_id=matrix_id,
                                                   cutoffs=cutoffs)
                        pwms.append(pwm)
                    matrix_id = m_header.group(1)
                    name = m_header.group(2)
                    n_matrix = 0
                    values = []
                    cutoffs = {}
                    expect_flag = 1  # found header, expect matrix line next
                elif m_matrix:
                    base = m_matrix.group(1)
                    if base != bases[n_matrix]:
                        raise PwmsMotifScanFormatError(line_num, line)
                    tmp_values = m_matrix.group(2).split()
                    try:
                        values.append(list(map(float, tmp_values)))
                    except (ValueError, TypeError):
                        raise PwmsMotifScanFormatError(line_num, line)
                    n_matrix += 1
                    if n_matrix == 4:
                        # got 4 matrix line, expect cutoff line next
                        expect_flag = 2
                elif m_cutoff:
                    p = m_cutoff.group(1)
                    cutoff = m_cutoff.group(2)
                    cutoffs[p] = float(cutoff)
                    # got first cutoff, expect cutoff or header next
                    if expect_flag == 2:
                        expect_flag = 3

            # check whether the last matrix is complete
            if expect_flag == 1 or expect_flag == 2:
                raise PwmsMotifScanFormatError(line_num + 1, '')
            if expect_flag == 3:
                pwm = PositionWeightMatrix(values=values, name=name,
                                           matrix_id=matrix_id,
                                           cutoffs=cutoffs)
                pwms.append(pwm)
        self.extend(pwms)
        logger.debug(f"Found {len(pwms)} MotifScan PWMs")


def load_installed_pfms(name):
    """Load a pre-installed motif PFMs set.

    Parameters
    ----------
    name : str
        Name of the pre-installed motif PFMs set to be loaded.

    Return
    ------
    pfms : `MotifPfms`
        Loaded PFMs of the motif set.

    Raises
    ------
    PfmsFileNotFoundError
        If the motif PFMs file does not exists.
    """
    logger.info(f"Loading motif PFMs set {name!r}")
    motif_dir = Config().get_motif_path(name)
    pfms_path = pfms_path_fmt.format(motif_dir, name)
    if os.path.isfile(pfms_path):
        pfms = MotifPfms(name=name)
        pfms.read_pfms(path=pfms_path, format='jaspar')
    else:
        raise PfmsFileNotFoundError(name)
    return pfms


def load_built_pwms(name, genome):
    """Load built motif PWMs.

    Parameters
    ----------
    name : str
        Name of the built motif PWMs set to be loaded.
    genome : str
        Genome assembly name under which these PWMs are built.

    Raises
    ------
    PwmsFileNotFoundError
        If the motif PWMs file does not exists
    """
    logger.info(
        f"Loading motif PWMs set {name!r} under genome {genome!r}")
    motif_dir = Config().get_motif_path(name)
    pwms_path = pwms_path_fmt.format(motif_dir, name, genome)
    pwms = MotifPwms(name=name, genome=genome)
    if os.path.isfile(pwms_path):
        pwms.read_motifscan_pwms(pwms_path)
    else:
        raise PwmsFileNotFoundError(name, genome)
    return pwms


def get_score_cutoffs(sampling_scores):
    """Get motif score cutoffs given the motif score background distributions.

    Parameters
    ----------
    sampling_scores : array_like
        The sampling scores of motifs in a shape of (n_motifs, n_sampling).
    """
    pwms_cutoffs = []
    n_pwms = len(sampling_scores)
    for i, scores in enumerate(sampling_scores):
        if len(scores) < 100:
            raise ValueError(
                "each motif must have at least 100 sampling scores")
        logger.debug(f"Getting cutoff: {i + 1}/{n_pwms}")
        pwm_cutoffs = {}
        n_scores = len(scores)
        n_bits = min(len(str(n_scores)), 7)
        scores.sort(reverse=True)
        for exponent in range(2, n_bits):
            cutoff = scores[int(n_scores * 0.1 ** exponent) - 1]
            pwm_cutoffs[f'1e-{exponent}'] = cutoff
        pwms_cutoffs.append(pwm_cutoffs)
    return pwms_cutoffs

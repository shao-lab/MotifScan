"""
motifscan.motif.matrix
----------------------

Module for motif position matrix.
"""

import numpy as np

from motifscan.genome import bases


class PositionMatrix:
    """Two-dimensional generic position matrix."""

    def __init__(self, values, name=None, matrix_id=None):
        if len(values) != 4:  # check if values have exactly 4 rows
            raise ValueError("values should have exactly 4 rows for A/C/G/T")

        self.matrix = np.asarray(values)
        if self.matrix.ndim != 2:
            raise ValueError("values should have 2 dimensions in (4 x N)")
        if not (np.issubdtype(self.matrix.dtype, np.integer) or np.issubdtype(
                self.matrix.dtype, np.floating)):
            raise ValueError("values should be integers or floating numbers")

        self._length = self.matrix.shape[1]
        if self._length == 0:
            raise ValueError("values should have at least 1 position per row")

        self.name = name
        self.matrix_id = matrix_id

    @property
    def shape(self):
        """Shape of the position matrix (4 x N)."""
        return self.matrix.shape

    @property
    def length(self):
        """Length of the position matrix (N)."""
        return self._length

    def __len__(self):
        return self._length

    def __str__(self):
        return "A {}\nC {}\nG {}\nT {}\n".format(*self.matrix)


class PositionFrequencyMatrix(PositionMatrix):
    """Two-dimensional position frequency matrix (PFM).

    Parameters
    ----------
    values : array_like
        Non-negative integers which compose the 4 x N position frequency
        matrix. Rows represent the vectors of N positions in the order of
        A, C, G, T.
    name : str, optional
        The name of the PFM.
    matrix_id : str, optional
        The id of the matrix.
    """

    def __init__(self, values, name=None, matrix_id=None):
        super().__init__(values, name, matrix_id)
        if not np.issubdtype(self.matrix.dtype, np.integer) or np.any(
                self.matrix < 0):
            raise ValueError("values in PFM should be non-negative integers")
        if np.any(self.matrix.sum(axis=0) == 0):
            raise ValueError("all values of a PFM position are 0")

    def to_ppm(self, normalize=True, pseudo=0.001):
        """Convert the position frequency matrix (PFM) to the position
        probability matrix (PPM).

        Parameters
        ----------
        normalize : bool, optional
            Whether normalize the PPM or not. Default=True.
        pseudo : float, optional
            The pseudo probability given to positions which contain 0.
            Only valid when `normalize` is True. Default=0.001.
            In this case, frequencies [0, 0, 10, 10] will be normalized to
            probabilities [0.001, 0.001, 0.499, 0.499].

        Returns
        -------
        ppm : PositionProbabilityMatrix
            Converted PPM.
        """
        ppm = PositionProbabilityMatrix(
            values=self.matrix / self.matrix.sum(axis=0), name=self.name,
            matrix_id=self.matrix_id)
        if normalize:
            ppm.normalize(pseudo)
        return ppm


class PositionProbabilityMatrix(PositionMatrix):
    """Two-dimensional position probability matrix (PPM).

    Parameters
    ----------
    values : array_like
        Non-negative probabilities which compose the 4 x N position probability
        matrix. Rows represent the vectors of N positions in the order of
        A, C, G, T.
    name : str, optional
        The name of the PPM.
    matrix_id : str, optional
        The id of the matrix.
    """

    def __init__(self, values, name=None, matrix_id=None):
        super().__init__(values, name, matrix_id)
        if np.any(self.matrix < 0):
            raise ValueError("values in PPM should be non-negative numbers")
        if np.any(self.matrix.sum(axis=0) == 0):
            raise ValueError("all values of a PPM position are 0")
        if not np.allclose(self.matrix.sum(axis=0), 1):
            raise ValueError("the sum probability of a PPM position is not 1")

    def normalize(self, pseudo=0.001):
        """Normalize the PPM, assign a pseudo probability for every zero value
        and normalize the sum probability to 1 for each position.

        Parameters
        ----------
        pseudo : float, optional
            The pseudo probability given to positions which contain 0.
            Only valid when `normalize` is True. Default=0.001.
            In this case, frequencies [0, 0, 10, 10] will be normalized to
            probabilities [0.001, 0.001, 0.499, 0.499].

        Notes
        -----
        This method is guaranteed only when the PPM is valid, which means the
        sum probability for each position should equal to or very close to 1.
        """
        if not 0 < pseudo < 0.25:
            raise ValueError("the range of pseudo should be (0, 0.25)")
        pseudo_count = pseudo / (1 - 4 * pseudo)
        zero_cols = np.any(self.matrix == 0, axis=0)
        self.matrix[:, zero_cols] += pseudo_count
        self.matrix = self.matrix / self.matrix.sum(axis=0)

    def to_pwm(self, bg_freq=None):
        """Convert the position probability matrix (PPM) to the position
        weight matrix (PWM).

        Parameters
        ----------
        bg_freq : dict of {str: float}, optional
            A dict object which stores the nucleotide frequencies used to
            calculate the weight. Default=None and [0.25, 0.25, 0.25, 0.25]
            will be used for A, C, G, T.

        Returns
        -------
        pwm : PositionWeightMatrix
            Converted PWM.
        """
        if bg_freq is None:
            bg_freq = {base: 0.25 for base in bases}
        bg = np.asarray([bg_freq[base] for base in bases]).reshape(4, 1)
        pwm = PositionWeightMatrix(
            values=np.around(np.log(self.matrix / bg), 5),
            name=self.name, matrix_id=self.matrix_id)
        return pwm


class PositionWeightMatrix(PositionMatrix):
    """Two-dimensional position weight matrix (PWM).

    Parameters
    ----------
    values : array_like
        Values which compose the 4 x N position weight matrix. Rows represent
        the vectors of N positions in the order of A, C, G, T.
    name : str, optional
        The name of the PWM.
    matrix_id : str, optional
        The id of the matrix.
    cutoffs : dict of {str: float}, optional
        Score cutoffs for motif occurrences under different p values.
    """

    def __init__(self, values, name=None, matrix_id=None, cutoffs=None):
        super().__init__(values, name, matrix_id)
        self._max_raw_score = None
        self._min_raw_score = None
        self.cutoffs = cutoffs

    def set_cutoff(self, p_value, cutoff):
        """Set motif score cutoff under certain p value."""
        if self.cutoffs is None:
            self.cutoffs = {}
        self.cutoffs[p_value] = cutoff

    @property
    def max_raw_score(self):
        """Calculate the possible maximum raw score for the PWM."""
        if self._max_raw_score is None:
            self._max_raw_score = self.matrix.max(axis=0).sum()
        return self._max_raw_score

    @property
    def min_raw_score(self):
        """Calculate the possible minimum raw score for the PWM."""
        if self._min_raw_score is None:
            self._min_raw_score = self.matrix.min(axis=0).sum()
        return self._min_raw_score

    def score(self, sequence):
        """Calculate the score for certain sequence, which is the ratio of the
        raw score to the possible maximum raw score (raw_score/max_raw_score).

        Parameters
        ----------
        sequence : str
            The sequence to be scored under the PWM.

        Returns
        -------
        score : float
            Score of the specified sequence.
        """
        if len(sequence) != self.length:
            raise ValueError("sequence should have the same length as the PWM")
        row_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        raw_score = 0
        for col_idx, nt in enumerate(sequence.upper()):
            try:
                raw_score += self.matrix[row_idx[nt], col_idx]
            except KeyError:
                continue
        score = raw_score / self.max_raw_score
        return score

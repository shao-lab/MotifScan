"""
motifscan.scanner
-----------------

Module for core motif scanner.
"""

import logging
from collections import namedtuple
from multiprocessing import Pool

from motifscan.motif.score import sliding_motif_score

logger = logging.getLogger(__name__)

MotifSite = namedtuple("MotifSite", ['start', 'score', 'strand'])


class Scanner:
    """
    Motif scanner for a list of genomic regions.

    Parameters
    ----------
    genome : `motifscan.genome.Genome`
        A motifscan genome object to extract sequences from.
    regions : list of `GenomicRegion`
        Specifies the genomic regions for scanning.
    window_size : int, optional
        If `window_size` <= 0, then the whole input region will be scanned.
        Otherwise, only scans the the specified window centered at the summit
        or midpoint of each genomic region.
    strand : {'both', '+', '-'}, optional
        Which strand is used to find motif occurrences, default='both'.
    p_value : str, optional
        Significance level for the motif score cutoff, default='1e-4'.
    remove_dup ：bool, optional
        If True (default), remove adjacent duplicated motif occurrences whose
        distance is less than the motif length.
    n_threads : int, optional
        Number of threads used in parallel, default=1.
    """

    def __init__(self, genome, regions, window_size=0, strand='both',
                 p_value='1e-4', remove_dup=True, n_threads=1):
        if window_size <= 0:
            self.window_size = 0
        else:
            self.window_size = window_size
        self.extend = window_size // 2
        if strand in ['both', '+', '-']:
            self.strand = strand
        else:
            raise ValueError(f"invalid strand option: {strand!r}")
        self.p_value = p_value
        self.remove_dup = remove_dup
        self.n_threads = n_threads
        self.seq_starts = []
        self.seq_ends = []
        self.sequences = []
        self.sequences_rc = None
        self._extract_seq(genome=genome, regions=regions)

    def _extract_seq(self, genome, regions):
        """Extract the sequences on the forward strand and record the
        underlying coordinates.
        """
        logger.debug("Extracting sequences on the forward strand")
        for region in regions:
            if self.window_size <= 0:
                seq_start = region.start
                seq_end = region.end
            else:
                seq_start = max(region.summit - self.extend, 0)
                seq_end = min(region.summit + self.extend,
                              genome.chrom_sizes[region.chrom])
            self.seq_starts.append(seq_start)
            self.seq_ends.append(seq_end)
            self.sequences.append(
                genome.fetch_sequence(region.chrom, seq_start, seq_end))
        if self.strand in ['both', '-']:
            self._rc_seq()

    def _rc_seq(self):
        """Reverse complement the sequences for the reverse strand."""
        logger.debug("Extracting sequences on the reverse strand")
        sequences_rc = []
        table = str.maketrans({'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                               'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'})
        for sequence in self.sequences:
            sequences_rc.append(sequence.translate(table)[::-1])
        self.sequences_rc = sequences_rc

    def _scan_pwm_by_strand(self, pwm, strand):
        """Scan the specified strand for motif occurrences."""
        if strand == '+':
            sequences = self.sequences
        elif strand == '-':
            sequences = self.sequences_rc
        else:
            raise ValueError(f"invalid strand option: {strand!r}")
        score_cutoff = pwm.cutoffs[self.p_value]
        sliding_scores = sliding_motif_score(
            pwm.matrix.tolist(), pwm.length, pwm.max_raw_score, sequences)
        sites = pick_motif_sites(
            sliding_scores, score_cutoff, self.seq_starts, strand=strand)
        if self.remove_dup:
            sites = deduplicate_motif_sites(sites, pwm.length)
        return sites

    def _scan_pwm(self, pwm):
        """Scan with a PWM and returns the detected motif sites."""
        logger.debug(f"Scanning {pwm.name}")
        if self.strand in ['both', '+']:
            sites_fwd = self._scan_pwm_by_strand(pwm=pwm, strand='+')
        if self.strand in ['both', '-']:
            sites_rev = self._scan_pwm_by_strand(pwm=pwm, strand='-')
        if self.strand == '+':
            sites = sites_fwd
        elif self.strand == '-':
            sites = sites_rev
        else:
            sites = [fwd + rev for fwd, rev in zip(sites_fwd, sites_rev)]
        logger.debug(f"Scanned {pwm.name}")
        return sites

    def scan_motifs(self, pwms):
        """Scan for motif occurrences given the motif PWMs.

        Parameters
        ----------
        pwms : `motifscan.motif.MotifPwms`
            Motif Pwms to be scanned.

        Returns
        -------
        motif_sites : nested list, (n_pwms, n_regions, )
            Identified motif occurrences.
        """
        # check if the corresponding motif score cutoff is set
        for pwm in pwms:
            try:
                pwm.cutoffs[self.p_value]
            except (TypeError, KeyError):
                raise ValueError(
                    f"PWM has no motif score cutoff set for P-value "
                    f"{self.p_value!r}")
        if self.n_threads > 1:
            pool = Pool(processes=self.n_threads)
            motif_sites = pool.map(self._scan_pwm, pwms)
        else:
            motif_sites = [self._scan_pwm(pwm) for pwm in pwms]
        return motif_sites


def pick_motif_sites(sliding_scores, cutoff, seq_starts, strand='+'):
    """Given the motif scores of all positions, pick motif occurrences out."""
    if strand not in ['+', '-']:
        raise ValueError(f"expect '+' or '-' for strand, got {strand!r}")
    sites = []
    for region_idx, scores in enumerate(sliding_scores):
        sites_by_region = []
        if strand == '-':
            scores = reversed(scores)
        for pos_idx, score in enumerate(scores):
            if score >= cutoff:
                start = seq_starts[region_idx] + pos_idx
                sites_by_region.append(MotifSite(start, score, strand))
        sites.append(sites_by_region)
    return sites


def deduplicate_motif_sites(sites, length):
    """Deduplicate adjacent motif occurrences with distance less than the motif
    length and only keep the one with higher motif score.
    """
    for sites_by_region in sites:
        idx = 0
        if len(sites_by_region) > 1:
            while idx + 1 < len(sites_by_region):
                site_curr = sites_by_region[idx]
                site_next = sites_by_region[idx + 1]
                if site_next.start - site_curr.start < length:
                    if site_curr.score >= site_next.score:
                        sites_by_region.pop(idx + 1)
                    else:
                        sites_by_region.pop(idx)
                else:
                    idx += 1
    return sites

"""
motifscan.scanner
-----------------

Module for core motif scanner.
"""

import logging
import os
from collections import namedtuple

from motifscan.motif.cscore import c_scan_motif

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
        n_threads = int(n_threads)
        n_cpu = os.cpu_count()
        if n_threads > n_cpu:
            logger.warning(f"Threads number exceed the number of CPUs, "
                           f"using {n_cpu} instead")
            n_threads = n_cpu
        if n_threads < 1:
            n_threads = 1
        self.n_threads = n_threads
        self.seq_starts = []
        self.seq_ends = []
        self.sequences = []
        self._extract_seq(genome=genome, regions=regions)

    def _extract_seq(self, genome, regions):
        """Extract the sequences on the forward strand and record the
        underlying coordinates.
        """
        logger.debug("Extracting sequences")
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

        logger.debug(f"Scanning motif PWMs")
        matrices = []
        cutoffs = []
        for pwm in pwms:
            matrices.append(pwm.matrix.tolist())
            cutoffs.append(pwm.cutoffs[self.p_value])

        if self.strand == '+':
            strand_arg = 1
        elif self.strand == '-':
            strand_arg = 2
        else:
            strand_arg = 3

        sites = c_scan_motif(matrices, cutoffs, self.sequences, strand_arg,
                             self.n_threads)

        motif_sites = make_motif_sites(sites, self.seq_starts)
        if self.remove_dup:
            lengths = [pwm.length for pwm in pwms]
            motif_sites = deduplicate_motif_sites(motif_sites, lengths)
        return motif_sites


def make_motif_sites(sites, seq_starts):
    """Given the pooled motif sites from C extension, rearrange motif sites by
    region(sequence), returns a nested list of shape (n_pwms, n_seqs, n_sites).
    """
    motif_sites = []
    for pwm_idx, sites_pwm in enumerate(sites):
        motif_sites.append([])
        for _ in seq_starts:
            motif_sites[pwm_idx].append([])
        for site in sites_pwm:
            seq_idx, pos_idx, score, strand = site
            start = seq_starts[seq_idx] + pos_idx
            if strand == 1:
                strand = '+'
            else:
                strand = '-'
            motif_sites[pwm_idx][seq_idx].append(
                MotifSite(start=start, score=score, strand=strand))
    return motif_sites


def _deduplicate_sites(sites, length):
    idx = 0
    if len(sites) > 1:
        while idx + 1 < len(sites):
            site_curr = sites[idx]
            site_next = sites[idx + 1]
            if site_next.start - site_curr.start < length:
                if site_curr.score >= site_next.score:
                    sites.pop(idx + 1)
                else:
                    sites.pop(idx)
            else:
                idx += 1


def deduplicate_motif_sites(motif_sites, lengths):
    """Deduplicate adjacent motif occurrences with distance less than the motif
    length and only keep the one with higher motif score.
    """
    motif_sites_dedup = []
    for sites_pwm, length in zip(motif_sites, lengths):
        sites_pwm_dedup = []
        for sites in sites_pwm:
            # split fwd sites and rev sites
            sites_fwd = []
            sites_rev = []
            for site in sites:
                if site.strand == '+':
                    sites_fwd.append(site)
                else:
                    sites_rev.append(site)
            _deduplicate_sites(sites_fwd, length)
            _deduplicate_sites(sites_rev, length)
            sites_dedup = sites_fwd + sites_rev
            sites_dedup.sort(key=lambda x: x.start)
            sites_pwm_dedup.append(sites_dedup)
        motif_sites_dedup.append(sites_pwm_dedup)
    return motif_sites_dedup

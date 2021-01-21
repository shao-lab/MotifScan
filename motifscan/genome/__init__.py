"""
motifscan.genome
----------------

This module contains classes and functions implemented for genome assemblies.
"""

import logging
import os

import numpy as np
import pysam

from motifscan.config import Config
from motifscan.exceptions import GenomeFileNotFoundError, BackgroundFormatError
from motifscan.genome.annotation import read_gene_annotation

logger = logging.getLogger(__name__)

bases = 'ACGT'

fasta_path_fmt = os.path.join("{0}", "{1}.fa")
bg_freq_path_fmt = os.path.join("{0}", "{1}_bg_freq.txt")
gene_path_fmt = os.path.join("{0}", "{1}_gene_annotation.txt")


class Genome:
    """Load a pre-installed genome assembly.

    This class handles basic genome properties like chromosome sizes,
    nucleotide frequencies and other advanced functions (fetching sequences,
    sampling regions etc.).

    Parameters
    ----------
    name : str
        The name of a pre-installed genome assembly to be loaded.
    path : str, optional
        Directory where the genome data files are located. If not specified,
        the path will be retrieved from the config file.

    Attributes
    ----------
    name : str
        The name of the genome assembly.
    path : str
        Directory where the genome data files are located.
    fa : `pysam.FastaFile`
        Fasta file handler.
    bg_freq : dict of {str: float}
        Nucleotide frequencies of whole genome background.
    genes : `motifscan.genome.annotation.Genes` or None
        Annotated genes, or None if the gene annotation file is not available.

    Raises
    ------
    GenomeFileNotFoundError
        If required genome data files are missing or broken.
    """

    def __init__(self, name, path=None):
        logger.info(f"Loading genome {name!r}")
        self.name = name
        self.path = path or Config().get_genome_path(self.name)
        self._fasta_path = fasta_path_fmt.format(self.path, self.name)
        self._bg_freq_path = bg_freq_path_fmt.format(self.path, self.name)
        self._gene_path = gene_path_fmt.format(self.path, self.name)
        if os.path.isfile(self._fasta_path):
            self.fa = pysam.FastaFile(self._fasta_path)
        else:
            raise GenomeFileNotFoundError(self.name, 'sequence')
        if os.path.isfile(self._bg_freq_path):
            self.bg_freq = read_bg_freq(self._bg_freq_path)
        else:
            raise GenomeFileNotFoundError(self.name, 'background frequency')
        if os.path.isfile(self._gene_path):
            self.genes = read_gene_annotation(self._gene_path)
        else:
            logger.warning("No gene annotation file found")
            self.genes = None

        self._chroms = None
        self._chrom_sizes = None

    def close(self):
        """Close the fasta handler."""
        self.fa.close()

    @property
    def chroms(self):
        """Returns sorted chromosome names of the genome assembly.

        Returns
        -------
        list of str
            Chromosome names (sorted) of the genome assembly.
        """
        if self._chroms is None:
            self._chroms = sorted(self.fa.references)
        return self._chroms

    @property
    def chrom_sizes(self):
        """Returns the sizes of chromosomes.

        Returns
        -------
        dict of {str: int}
            A dict object with key as the chromosome name and value as the
            size of that chromosome.
        """
        if self._chrom_sizes is None:
            self._chrom_sizes = {chrom: self.fa.get_reference_length(chrom) for
                                 chrom in self.chroms}
        return self._chrom_sizes

    def fetch_sequence(self, chrom, start, end):
        """Returns the sequence of given genomic region (0-based coordinates).
        Wrapper of `pysam.FastaFile.fetch`.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            Start position.
        end : int
            End position.

        Returns
        -------
        str
            The sequence of the specified region.
        """
        return self.fa.fetch(chrom, start, end)

    def random_sequences(self, n_times, length, max_n=0, random_seed=None):
        """Random sampling sequences from whole genome. The weight for each
        chromosome to be selected equals to the proportion of its chromosome
        size.

        Parameters
        ----------
        n_times: int
            Sampling times.
        length: int
            The length of sequences to generate.
        max_n: int, optional
            The maximal number of `N` base allowed in a sampled sequence,
            default=0.
        random_seed: int, optional
            random_seed: The seed used to set random state.

        Yields
        ------
        str
            Random sampled sequence of length `length`.
        """
        if random_seed is not None:
            logger.debug(f"Setting random seed: {random_seed}")
            np.random.seed(random_seed)
        chrom_sizes_sum = sum(self.chrom_sizes.values())
        chrom_weight = [self.chrom_sizes[chrom] / chrom_sizes_sum for chrom in
                        self.chroms]
        random_chroms = np.random.choice(self.chroms, size=n_times,
                                         p=chrom_weight)
        n_seq = 0
        n_loop = 0
        while n_seq < n_times:
            chrom = random_chroms[n_loop % n_times]
            start = np.random.randint(self.chrom_sizes[chrom] - length)
            seq = self.fetch_sequence(chrom, start, start + length)
            if seq.count('N') + seq.count('n') <= max_n:
                yield seq
                n_seq += 1
            n_loop += 1


def cal_bg_freq(path, skip_non_autosomes=True):
    """Calculate the nucleotide frequencies of whole genome background.

    Parameters
    ----------
    path : str
        A fasta file which contains the sequences of whole genome background.
    skip_non_autosomes : bool, optional
        If True (default), only count bases on autosomes. Sex chromosomes
        (chrX, chrY), mitochondria DNA (chrM) and other unconfirmed sequences
        (chrUn_*, *_alt*, *_hap*, *_random*) are excluded.

    Returns
    -------
    bg_freq : dict of {str: float}
        A dict object which stores the nucleotide frequencies of A, C, G, T.
    """
    logger.debug(f"Calculating nucleotide frequencies: {path}")
    bg_count = {base: 0 for base in bases}
    keywords_to_skip = ['chrX', 'chrY', 'chrM', 'chrUn_',
                        '_random', '_hap', '_alt']

    fa = pysam.FastaFile(path)
    for chrom in fa.references:
        if skip_non_autosomes:
            skip = False
            for keyword in keywords_to_skip:
                if keyword in chrom:
                    skip = True
                    break
            if skip:
                logger.debug(f"Skipped: {chrom}")
                continue
        logger.debug(f"Processing: {chrom}")
        sequence = fa.fetch(chrom, 0, fa.get_reference_length(chrom)).upper()
        for base in bases:
            bg_count[base] += sequence.count(base)
    fa.close()

    total_count = sum(bg_count.values())
    bg_freq = {base: round(bg_count[base] / total_count, 5) for base in bases}
    return bg_freq


def write_bg_freq(path, bg_freq):
    """Write the nucleotide frequencies of whole genome background.

    Parameters
    ----------
    path : str
        The file path to write the nucleotide frequencies.
    bg_freq: dict of {str: float}
        A dict object which stores the nucleotide frequencies of A, C, G, T.
    """
    logger.debug(f"Writing nucleotide frequencies to {path}")
    with open(path, 'w') as f_out:
        for base in bases:
            f_out.write(f"{base}\t{bg_freq[base]}\n")


def read_bg_freq(path):
    """Read the nucleotide frequencies of whole genome background.

    Parameters
    ----------
    path : str
        The file path to read the nucleotide frequencies.

    Returns
    -------
    bg_freq : dict of {str: float}
        A dict object which stores the nucleotide frequencies of A, C, G, T.

    Raises
    ------
    BackgroundFormatError
        If the frequency file not strictly follows the format.
    """
    logger.debug(f"Reading nucleotide frequencies from {path}")
    bg_freq = {}
    with open(path, 'r') as f_in:
        for idx, expected in enumerate(bases):
            line = f_in.readline().strip()
            base, freq = line.split('\t')
            if base != expected:
                raise BackgroundFormatError(idx + 1, line)
            try:
                bg_freq[base] = float(freq)
            except (ValueError, TypeError):
                raise BackgroundFormatError(idx + 1, line)
    return bg_freq

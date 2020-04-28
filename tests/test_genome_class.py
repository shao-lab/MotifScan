import os

import pytest

from motifscan.genome import Genome, cal_bg_freq, read_bg_freq, write_bg_freq
from motifscan.exceptions import GenomeFileNotFoundError, \
    BackgroundFormatError


def test_genome_init(genome_root):
    genome_path = os.path.join(genome_root, 'test')
    genome = Genome(name='test', path=genome_path)
    assert genome.name == 'test'
    assert genome.path == genome_path
    assert genome.chroms == ['chr1', 'chr2', 'chrM', 'chrX']
    assert genome.chrom_sizes['chr1'] == 10
    assert genome.chrom_sizes['chr2'] == 17
    assert genome.chrom_sizes['chrM'] == 15
    assert genome.chrom_sizes['chrX'] == 16
    assert genome.fetch_sequence('chr1', 0, 10) == 'AaTtCcGgNn'
    assert genome.fetch_sequence('chr2', 0, 17) == 'AAAaCCccTTtGNNNNN'
    assert genome.fetch_sequence('chrM', 0, 15) == 'AaaaaAAAAAAAnnn'
    assert genome.fetch_sequence('chrX', 0, 16) == 'AAACCTACNNTnggAC'
    seq_random = list(genome.random_sequences(n_times=3, length=5))
    assert len(seq_random) == 3
    seq_random = list(
        genome.random_sequences(n_times=3, length=5, random_seed=1))
    assert seq_random == ['AAAAA', 'AaTtC', 'AAAaC']
    genome.close()
    assert genome.fa.closed
    with pytest.raises(GenomeFileNotFoundError):
        Genome(name='unknown_genome', path=genome_path)


def test_genome_cal_bg_freq(genome_root):
    fasta_path = os.path.join(genome_root, 'test', 'test.fa')
    bg_freq = cal_bg_freq(fasta_path, skip_non_autosomes=True)
    assert sorted(bg_freq.keys()) == ['A', 'C', 'G', 'T']
    assert bg_freq['A'] == pytest.approx(0.3)
    assert bg_freq['C'] == pytest.approx(0.3)
    assert bg_freq['G'] == pytest.approx(0.15)
    assert bg_freq['T'] == pytest.approx(0.25)
    bg_freq = cal_bg_freq(fasta_path, skip_non_autosomes=False)
    assert bg_freq['A'] == pytest.approx(0.51111)
    assert bg_freq['C'] == pytest.approx(0.22222)
    assert bg_freq['G'] == pytest.approx(0.11111)
    assert bg_freq['T'] == pytest.approx(0.15556)


def test_read_bg_freq(genome_root):
    bg_freq_path = os.path.join(genome_root, 'test', 'test_bg_freq.txt')
    bg_freq = read_bg_freq(bg_freq_path)
    assert sorted(bg_freq.keys()) == ['A', 'C', 'G', 'T']
    assert bg_freq['A'] == 0.3
    assert bg_freq['C'] == 0.3
    assert bg_freq['G'] == 0.15
    assert bg_freq['T'] == 0.25


def test_bg_freq_format(genome_root):
    genome_path = os.path.join(genome_root, 'bad')
    with pytest.raises(BackgroundFormatError):
        read_bg_freq(os.path.join(genome_path, 'test_bg_freq_bad1.txt'))
    with pytest.raises(BackgroundFormatError):
        read_bg_freq(os.path.join(genome_path, 'test_bg_freq_bad2.txt'))


def test_write_bg_freq(genome_root, tmp_dir):
    fasta_path = os.path.join(genome_root, 'test', 'test.fa')
    bg_freq = cal_bg_freq(fasta_path, skip_non_autosomes=True)
    bg_freq_path = os.path.join(tmp_dir, 'test_bg_freq.txt')
    write_bg_freq(bg_freq_path, bg_freq)
    bg_freq = read_bg_freq(bg_freq_path)
    assert bg_freq['A'] == 0.3
    assert bg_freq['C'] == 0.3
    assert bg_freq['G'] == 0.15
    assert bg_freq['T'] == 0.25

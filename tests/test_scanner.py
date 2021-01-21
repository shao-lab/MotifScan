import os
import pytest

from motifscan.genome import Genome
from motifscan.region import GenomicRegion
from motifscan.scanner import Scanner, MotifSite, deduplicate_motif_sites
from motifscan.motif import MotifPwms
from motifscan.motif.matrix import PositionWeightMatrix


def test_scanner_init(genome_root):
    genome_path = os.path.join(genome_root, 'test')
    genome = Genome(name='test', path=genome_path)
    regions = [GenomicRegion(chrom='chr1', start=2, end=5)]
    scanner = Scanner(genome=genome, regions=regions, window_size=0)
    assert scanner.window_size == 0
    assert scanner.sequences == ['TtC']
    assert scanner.seq_starts == [2]
    assert scanner.seq_ends == [5]
    scanner = Scanner(genome=genome, regions=regions, window_size=4,
                      strand='+')
    assert scanner.sequences == ['aTtC']
    assert scanner.seq_starts == [1]
    assert scanner.seq_ends == [5]
    with pytest.raises(ValueError):
        Scanner(genome=genome, regions=regions, window_size=0, strand='*')


def test_scan_motif(genome_root):
    genome_path = os.path.join(genome_root, 'test')
    genome = Genome(name='test', path=genome_path)
    regions = [GenomicRegion(chrom='chr1', start=2, end=5)]
    pwms = MotifPwms()
    pwm = PositionWeightMatrix([[1, 0], [0, 1], [0, 0], [1, 0]],
                               cutoffs={'1e-3': 0.5, '1e-4': 1})
    pwms.append(pwm)
    scanner = Scanner(genome=genome, regions=regions, window_size=4,
                      p_value='1e-4')
    sites = scanner.scan_motifs(pwms)
    assert len(sites[0]) == 1
    assert len(sites[0][0]) == 1
    scanner = Scanner(genome=genome, regions=regions, window_size=4,
                      p_value='1e-3')
    sites = scanner.scan_motifs(pwms)
    assert len(sites[0][0]) == 3
    print(sites[0][0])
    scanner = Scanner(genome=genome, regions=regions, window_size=4,
                      p_value='1e-2')
    with pytest.raises(ValueError):
        scanner.scan_motifs(pwms)
    scanner = Scanner(genome=genome, regions=regions, window_size=4,
                      p_value='1e-3', remove_dup=False)
    sites = scanner.scan_motifs(pwms)
    assert len(sites[0][0]) == 5


def test_deduplicate_motif_sites():
    site1 = MotifSite(1, 1, '+')
    site2 = MotifSite(3, 0.8, '+')
    site3 = MotifSite(1, 1, '-')
    site4 = MotifSite(2, 3, '-')
    site5 = MotifSite(5, 1, '+')
    motif_sites = \
        deduplicate_motif_sites([[[site1, site2, site3, site4, site5]]], [3])
    assert len(motif_sites) == 1
    assert len(motif_sites[0]) == 1
    assert len(motif_sites[0][0]) == 3
    assert motif_sites[0][0][0].start == 1
    assert motif_sites[0][0][0].strand == '+'
    assert motif_sites[0][0][1].start == 2
    assert motif_sites[0][0][1].strand == '-'
    assert motif_sites[0][0][2].start == 5
    assert motif_sites[0][0][2].strand == '+'

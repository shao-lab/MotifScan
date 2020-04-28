import os

from motifscan.region import GenomicRegion
from motifscan.region.utils import overlap_with, subset_by_location, \
    generate_control_regions
from motifscan.genome.annotation import read_gene_annotation


def test_overlap_with():
    assert not overlap_with([], 1, 100)
    intervals = [[1, 5], [3, 8], [10, 12]]
    assert overlap_with(intervals, 1, 3)
    assert overlap_with(intervals, 3, 6)
    assert not overlap_with(intervals, 0, 1)
    assert not overlap_with(intervals, 8, 10)
    assert not overlap_with(intervals, 15, 20)


def test_subset_by_location(genome_root):
    genes = read_gene_annotation(
        os.path.join(genome_root, 'test', 'test_gene_annotation.txt'))
    regions = [GenomicRegion(chrom='chr1', start=9868, end=13868)]
    regions_subset = subset_by_location(regions=regions, genes=genes,
                                        location='promoter')
    assert len(regions_subset) == 1
    regions_subset = subset_by_location(regions=regions, genes=genes,
                                        location='distal')
    assert len(regions_subset) == 0
    regions = [GenomicRegion(chrom='chr1', start=9868, end=10868)]
    regions_subset = subset_by_location(regions=regions, genes=genes,
                                        location='promoter', upstream=1000)
    assert len(regions_subset) == 0


def test_generate_control_regions(genome_root):
    genes = read_gene_annotation(
        os.path.join(genome_root, 'test', 'test_gene_annotation.txt'))
    regions = [GenomicRegion(chrom='chr1', start=9868, end=13868),
               GenomicRegion(chrom='chr1', start=50000, end=51000),
               GenomicRegion(chrom='chr1', start=17200, end=17500)]
    regions_random = generate_control_regions(
        n_random=2, regions=regions, chrom_size={'chr1': 1000000},
        genes=genes)
    assert len(regions_random) == 6
    regions_random = generate_control_regions(
        n_random=2, regions=regions, chrom_size={'chr1': 1000000},
        genes=genes, random_seed=1)
    assert len(regions_random) == 6
    regions_random = generate_control_regions(
        n_random=2, regions=regions, chrom_size={'chr1': 1000000})
    assert len(regions_random) == 6

"""
motifscan.region.utils
----------------------

Module for genomic regions utilities.
"""

import logging
import random

from motifscan.region import GenomicRegion

logger = logging.getLogger(__name__)


def overlap_with(intervals, start, end):
    """Returns whether the interval [start, end) overlaps with the intervals.
    The input intervals have been sorted in advance for binary search.

    Parameters
    ----------
    intervals : list
        Intervals.
    start : int
        Start coordinate.
    end : int
        End coordinate.

    Returns
    -------
    bool
        Overlaps or not.
    """
    if not intervals:
        return False
    left = 0
    right = len(intervals) - 1
    while left <= right:
        mid = (left + right) // 2
        start_ref = intervals[mid][0]
        end_ref = intervals[mid][1]
        if not (end <= start_ref or start >= end_ref):
            return True
        elif start >= end_ref:
            left = mid + 1
        elif end <= start_ref:
            right = mid - 1
    return False


def subset_by_location(regions, genes, location, upstream=2000,
                       downstream=2000):
    """Filter genomic regions by location (promoter or distal).

    Parameters
    ----------
    regions: list
        Genomic regions to be filtered.
    genes : `motifscan.genome.annotation.Genes`
        Gene annotations.
    location : {'promoter', 'distal'}
        Only keep regions located at promoters or distal regions.
    upstream : int, optional
        TSS upstream distance to define promoters.
    downstream : int, optional
        TSS downstream distance to define promoters.

    Returns
    -------
        Filtered genomic regions subset.
    """
    filtered_regions = []
    promoters = {}
    for region in regions:
        if region.chrom not in promoters:
            # get promoter intervals of that chromosome
            promoters[region.chrom] = []
            for gene in genes.fetch(region.chrom):
                promoters[region.chrom].append(
                    gene.promoter(upstream, downstream))
            promoters[region.chrom].sort()
        overlap = overlap_with(promoters[region.chrom], region.start,
                               region.end)
        if not overlap ^ (location == 'promoter'):  # xor
            filtered_regions.append(region)
    return filtered_regions


def generate_control_regions(n_random, regions, chrom_size, genes=None,
                             random_seed=None):
    """Generate random control regions for the given reference regions.
    The length and chromosome distribution for each region are controlled to
    match the reference regions. If `genes` and location is specified,

    Parameters
    ----------
    n_random : int
        How many random control regions is generated for each reference region.
    regions : list
        Reference regions.
    chrom_size : dict
        Chromosome sizes.
    genes : `motifscan.genome.annotation.Genes`
        Gene annotations, controls the distances of regions to nearest genes
        when generating control regions. In this way, the control regions
        generated for a promoter region is also picked from promoters.
        If not specified, ignore annotations and random sampling from the
        whole background.
    random_seed: int, optional
        random_seed: The seed used to set random state.
    """
    if random_seed is not None:
        logger.debug(f"Setting random seed: {random_seed}")
        random.seed(random_seed)
    regions_control = []
    for region in regions:
        length = region.end - region.start
        if genes is None:
            for _ in range(n_random):
                start = random.randint(0, chrom_size[region.chrom] - length)
                regions_control.append(
                    GenomicRegion(chrom=region.chrom, start=start,
                                  end=start + length))
        else:
            genes_chrom = genes.fetch(region.chrom)
            if not genes_chrom:
                continue
            distance = dis_to_nearest_gene(region, genes_chrom)
            n_generated = 0
            while n_generated < n_random:
                if distance is None:
                    # randomize a distance if it's too far away gene TSSs
                    distance = random.randint(10000, 100000)
                gene_random = random.choice(genes_chrom)
                if gene_random.strand == '+':
                    start = gene_random.tss + distance
                else:
                    start = gene_random.tss - distance
                if (start >= 0) and (
                        start + length <= chrom_size[region.chrom]):
                    regions_control.append(
                        GenomicRegion(chrom=region.chrom, start=start,
                                      end=start + length))
                    n_generated += 1
    return regions_control


def dis_to_nearest_gene(region, genes, distance_cutoff=10000):
    """Compute the distance of  a specified region to the nearest gene TSS.

    Parameters
    ----------
    region : GenomicRegion
        GenomicRegion object.
    genes : list
        List of genes on that chromosome.
    distance_cutoff : int
        Maximal distance cutoff.

    Returns
    -------
    int or None
        A signed distance, positive if the region is located at the downstream
        of the nearest gene and vice versa. Or None if the distance is beyond
        the distance cutoff.
    """
    dis_min = distance_cutoff
    target_gene = None
    for gene in genes:
        tmp_dis = region.start - gene.tss
        if abs(tmp_dis) < dis_min:
            dis_min = tmp_dis
            target_gene = gene
    if target_gene is None:
        return None
    else:
        if target_gene.strand == '-':
            # signed value, positive for downstream and negative for upstream
            dis_min = -dis_min
        return dis_min

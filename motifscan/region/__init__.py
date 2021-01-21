"""
motifscan.region
----------------

Module for MotifScan genomic regions.
"""

import logging
from motifscan.region.parsers import get_region_parser

logger = logging.getLogger(__name__)

REGION_FORMATS = ['bed', 'bed3-summit', 'macs', 'macs2', 'narrowpeak',
                  'broadpeak', 'manorm']


class GenomicRegion:
    """Class for a genomic region.

    Parameters
    ----------
    chrom : str
        The chromosome name of the region.
    start : int
        The start coordinate of the region.
    end : int
        The end coordinate of the region.
    summit : int, optional
        The summit coordinate of the region. If not specified, the middle point
        will be taken as the summit.
    score : float, optional
        The score of the region.

    Attributes
    ----------
    chrom : str
        The chromosome name of the region.
    start : int
        The start coordinate of the region.
    end : int
        The end coordinate of the region.
    summit : int
        The summit coordinate of the region.
    score : float or None
        The score of the region or None if not specified.

    Notes
    -----
    The coordinates are 0-based, which means the region range is [start, end)
    and `start` <= `summit` < `end` is strictly inspected.
    """

    def __init__(self, chrom, start, end, summit=None, score=None):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        if self.start >= self.end:
            raise ValueError(
                f"expect start < end, got: start={start} end={end}")
        if summit is not None:
            self.summit = int(summit)
        else:
            self.summit = (self.start + self.end) // 2
        if not self.start <= self.summit < self.end:
            logger.warning(f"expect start <= summit < end, got chrom={chrom} "
                           f"start={start} summit={summit} end={end}")
        self.score = score

    def __repr__(self):
        return f"GenomicRegion({self.chrom}:{self.start}-{self.end})"


def load_motifscan_regions(path, format='bed'):
    """Load genomic regions from the specified path.

    Parameters
    ----------
    path : str
        Path to load the genomic regions.
    format : str, optional
        File format, default='bed'.

    Returns
    -------
    regions : list of `GenomicRegion`
        Loaded genomic regions.
    """
    logger.info(f"Loading genomic regions from {path} [{format}]")
    parser = get_region_parser(format)()
    regions = []
    for chrom, start, end, summit, score in parser.parse(path):
        regions.append(GenomicRegion(chrom, start, end, summit, score))
    logger.info(f"Loaded {len(regions)} genomic regions")
    return regions

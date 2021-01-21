import pytest
from motifscan.region import GenomicRegion


def test_region_init():
    region = GenomicRegion(chrom='chr1', start=1, end=100, score=10)
    assert region.summit == 50
    assert region.score == 10
    region = GenomicRegion(chrom='chr1', start=1, end=100, summit=51)
    assert region.summit == 51
    with pytest.raises(ValueError):
        GenomicRegion(chrom='chr1', start=1, end=1)
    with pytest.raises(ValueError):
        GenomicRegion(chrom='chr1', start=2, end=1)

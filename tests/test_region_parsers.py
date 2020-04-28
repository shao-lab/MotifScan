import os

import pytest

from motifscan.exceptions import RegionFileFormatError
from motifscan.region import load_motifscan_regions


def test_unsupported_format(region_root):
    with pytest.raises(ValueError):
        load_motifscan_regions(os.path.join(region_root, 'test_regions.bed'),
                               format='unknown_format')


def test_invalid_format(region_root):
    with pytest.raises(RegionFileFormatError):
        load_motifscan_regions(os.path.join(region_root, 'test_regions.bed'),
                               format='bed3-summit')


def test_bed_parser(region_root):
    regions = load_motifscan_regions(
        os.path.join(region_root, 'test_regions.bed'),
        format='bed')
    assert len(regions) == 4
    assert regions[0].score == 100
    assert regions[-1].score == None


def test_bed3_summit_parser(region_root):
    regions = load_motifscan_regions(
        os.path.join(region_root, 'test_regions_summit.bed'),
        format='bed3-summit')
    assert len(regions) == 4


def test_macs_parser(region_root):
    regions = load_motifscan_regions(
        os.path.join(region_root, 'test_regions_macs.xls'),
        format='macs')
    assert len(regions) == 9
    assert regions[0].score == 100.54


def test_macs2_parser(region_root):
    regions = load_motifscan_regions(
        os.path.join(region_root, 'test_regions_macs2.xls'),
        format='macs2')
    assert len(regions) == 10
    assert regions[0].score == 33.13990


def test_narrowpeak_parser(region_root):
    regions = load_motifscan_regions(
        os.path.join(region_root, 'test_regions.narrowPeak'),
        format='narrowpeak')
    assert len(regions) == 10
    assert regions[0].score == 331


def test_broadpeak_parser(region_root):
    regions = load_motifscan_regions(
        os.path.join(region_root, 'test_regions.broadPeak'),
        format='broadpeak')
    assert len(regions) == 10
    assert regions[0].score == 331


def test_manorm_parser(region_root):
    regions = load_motifscan_regions(
        os.path.join(region_root, 'test_regions_manorm.xls'),
        format='manorm')
    assert len(regions) == 4
    assert regions[0].score == 6.77003

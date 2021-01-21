import os

import pytest

from motifscan.exceptions import MotifSetNotFoundError, \
    PfmsJasparFormatError, PwmsMotifScanFormatError
from motifscan.motif import MotifPfms, MotifPwms, load_installed_pfms, \
    load_built_pwms
from motifscan.motif.matrix import PositionFrequencyMatrix, \
    PositionWeightMatrix
from motifscan.motif import get_score_cutoffs


def test_pfms_init(motif_root):
    with pytest.raises(ValueError):
        MotifPfms(
            pfms=[PositionFrequencyMatrix([[1], [2], [3], [4]]), 's', True])
    name = 'motif_set'
    pfm = PositionFrequencyMatrix([[1], [2], [3], [4]])
    pfms = MotifPfms(pfms=[pfm], name=name)
    assert pfms.name == name
    assert list(pfms) == [pfm]
    assert len(pfms) == 1
    with pytest.raises(MotifSetNotFoundError):
        load_installed_pfms('uninstalled_motif_set')


def test_read_jaspar_pfms(motif_root):
    pfms = MotifPfms()
    pfms.read_pfms(os.path.join(motif_root, 'test', 'test_pfms.jaspar'),
                   format='jaspar')
    assert len(pfms) == 2
    with pytest.raises(ValueError):
        pfms.read_pfms(os.path.join(motif_root, 'test', 'test_pfms.jaspar'),
                       format='jas')


def test_pfms_jaspar_format(motif_root):
    pfms = MotifPfms()
    with pytest.raises(PfmsJasparFormatError):
        pfms.read_pfms(
            os.path.join(motif_root, 'bad', 'test_pfms_bad1.jaspar'),
            format='jaspar')
    with pytest.raises(PfmsJasparFormatError):
        pfms.read_pfms(
            os.path.join(motif_root, 'bad', 'test_pfms_bad2.jaspar'),
            format='jaspar')
    with pytest.raises(PfmsJasparFormatError):
        pfms.read_pfms(
            os.path.join(motif_root, 'bad', 'test_pfms_bad3.jaspar'),
            format='jaspar')
    with pytest.raises(PfmsJasparFormatError):
        pfms.read_pfms(
            os.path.join(motif_root, 'bad', 'test_pfms_bad4.jaspar'),
            format='jaspar')


def test_pwms_init(motif_root):
    with pytest.raises(ValueError):
        MotifPwms(
            pwms=[PositionWeightMatrix([[1], [2], [3], [4]]), 's', True])
    pwm = PositionWeightMatrix([[1], [2], [3], [4]])
    name = 'motif_set'
    pwms = MotifPwms(pwms=[pwm], name=name)
    assert len(pwms) == 1
    assert pwms.name == name
    pwms = MotifPwms(name=name)
    pwms.read_motifscan_pwms(
        os.path.join(motif_root, 'test', 'test_pwms.motifscan'))
    assert len(pwms) == 2
    with pytest.raises(MotifSetNotFoundError):
        load_built_pwms('uninstalled_motif_set', 'hg19')


def test_pwms_motifscan_format(motif_root):
    pwms = MotifPwms()
    with pytest.raises(PwmsMotifScanFormatError):
        pwms.read_motifscan_pwms(
            os.path.join(motif_root, 'bad', 'test_pwms_bad1.motifscan'))
    with pytest.raises(PwmsMotifScanFormatError):
        pwms.read_motifscan_pwms(
            os.path.join(motif_root, 'bad', 'test_pwms_bad2.motifscan'))
    with pytest.raises(PwmsMotifScanFormatError):
        pwms.read_motifscan_pwms(
            os.path.join(motif_root, 'bad', 'test_pwms_bad3.motifscan'))
    with pytest.raises(PwmsMotifScanFormatError):
        pwms.read_motifscan_pwms(
            os.path.join(motif_root, 'bad', 'test_pwms_bad4.motifscan'))
    with pytest.raises(PwmsMotifScanFormatError):
        pwms.read_motifscan_pwms(
            os.path.join(motif_root, 'bad', 'test_pwms_bad5.motifscan'))
    with pytest.raises(PwmsMotifScanFormatError):
        pwms.read_motifscan_pwms(
            os.path.join(motif_root, 'bad', 'test_pwms_bad6.motifscan'))


def test_pwms_write_motifscan_pwms(tmp_dir):
    pwm = PositionWeightMatrix([[1], [2], [3], [4]], name='motif1',
                               matrix_id='id1',
                               cutoffs={'1e-2': 0.1, '1e-3': 0.4})
    name = 'motif_set'
    pwms = MotifPwms(pwms=[pwm], name=name)
    pwms_path = os.path.join(tmp_dir, 'test_pwm.motifscan')
    pwms.write_motifscan_pwms(pwms_path)
    assert os.path.isfile(pwms_path)


def test_get_score_cutoffs():
    cutoffs = get_score_cutoffs([list(range(0, 1000000))])
    assert len(cutoffs) == 1
    assert len(cutoffs[0]) == 5
    assert cutoffs[0]['1e-2'] == 990000
    assert cutoffs[0]['1e-3'] == 999000
    assert cutoffs[0]['1e-4'] == 999900
    assert cutoffs[0]['1e-5'] == 999990
    assert cutoffs[0]['1e-6'] == 999999
    with pytest.raises(ValueError):
        get_score_cutoffs([[1], [2]])
    with pytest.raises(TypeError):
        get_score_cutoffs([1, 2, 3])
    with pytest.raises(ValueError):
        get_score_cutoffs([[1, 2, 3]])

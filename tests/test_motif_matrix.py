import pytest

from motifscan.motif.matrix import PositionMatrix, PositionFrequencyMatrix, \
    PositionProbabilityMatrix, PositionWeightMatrix


def test_matrix_init():
    matrix = PositionMatrix(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]], name='matrix1',
        matrix_id='matrix_1')
    assert matrix.name == 'matrix1'
    assert matrix.matrix_id == 'matrix_1'
    assert matrix.shape == (4, 3)
    assert matrix.length == 3
    assert len(matrix) == 3
    with pytest.raises(ValueError):
        PositionMatrix([])
    with pytest.raises(ValueError):
        PositionMatrix([1, 2, 3, 4])
    with pytest.raises(ValueError):
        PositionMatrix([[], [], [], []])
    with pytest.raises(ValueError):
        PositionMatrix(['s', 'a', 1, 2])
    with pytest.raises(ValueError):
        PositionMatrix([['s'], ['a'], ['o'], ['e']])


def test_pfm_init():
    with pytest.raises(ValueError):
        PositionFrequencyMatrix([[1], [0.4], [7], [10]])
    with pytest.raises(ValueError):
        PositionFrequencyMatrix([[-1], [4], [7], [10]])
    with pytest.raises(ValueError):
        PositionFrequencyMatrix([[0], [0], [0], [0]])


def test_pfm_to_ppm():
    pfm = PositionFrequencyMatrix([[1, 1], [1, 2], [1, 2], [1, 0]])
    ppm = pfm.to_ppm(normalize=False)
    assert ppm.matrix[0][0] == pytest.approx(0.25)
    assert ppm.matrix[0][1] == pytest.approx(0.2)
    assert ppm.matrix[1][0] == pytest.approx(0.25)
    assert ppm.matrix[1][1] == pytest.approx(0.4)
    assert ppm.matrix[3][1] == pytest.approx(0)
    ppm = pfm.to_ppm(normalize=True, pseudo=0.001)
    assert ppm.matrix[0][0] == pytest.approx(0.25)
    assert ppm.matrix[0][1] == pytest.approx(0.2002)
    assert ppm.matrix[1][0] == pytest.approx(0.25)
    assert ppm.matrix[1][1] == pytest.approx(0.3994)
    assert ppm.matrix[3][1] == pytest.approx(0.001)


def test_ppm_init():
    with pytest.raises(ValueError):
        PositionProbabilityMatrix([[0], [0], [0], [0]])
    with pytest.raises(ValueError):
        PositionProbabilityMatrix([[0], [0.2], [-0.1], [0.9]])
    with pytest.raises(ValueError):
        PositionProbabilityMatrix([[0.3], [0.2], [0.5], [0.3]])


def test_ppm_normalize():
    ppm = PositionProbabilityMatrix(
        [[0.2, 0.2], [0.2, 0.2], [0.3, 0.6], [0.3, 0]])
    with pytest.raises(ValueError):
        ppm.normalize(pseudo=1)
    ppm.normalize(pseudo=0.001)
    assert ppm.matrix[0][0] == pytest.approx(0.2)
    assert ppm.matrix[0][1] == pytest.approx(0.2002)
    assert ppm.matrix[2][0] == pytest.approx(0.3)
    assert ppm.matrix[2][1] == pytest.approx(0.5986)
    assert ppm.matrix[3][1] == pytest.approx(0.001)


def test_ppm_to_pwm():
    ppm = PositionProbabilityMatrix(
        [[0.2, 0.2], [0.2, 0.2], [0.3, 0.6], [0.3, 0]])
    ppm.normalize(pseudo=0.001)
    pwm = ppm.to_pwm()
    assert pwm.matrix[0][0] == pytest.approx(-0.22314)
    assert pwm.matrix[0][1] == pytest.approx(-0.22214)
    assert pwm.matrix[2][0] == pytest.approx(0.18232)
    assert pwm.matrix[2][1] == pytest.approx(0.87313)
    assert pwm.matrix[3][1] == pytest.approx(-5.52146)
    pwm = ppm.to_pwm(bg_freq={'A': 0.22, 'C': 0.23, 'G': 0.28, 'T': 0.27})
    assert pwm.matrix[0][0] == pytest.approx(-0.09531)
    assert pwm.matrix[0][1] == pytest.approx(-0.09431)
    assert pwm.matrix[1][0] == pytest.approx(-0.13976)
    assert pwm.matrix[1][1] == pytest.approx(-0.13876)
    assert pwm.matrix[2][0] == pytest.approx(0.06899)
    assert pwm.matrix[2][1] == pytest.approx(0.75980)
    assert pwm.matrix[3][1] == pytest.approx(-5.59842)


def test_pwm_raw_score():
    pwm = PositionWeightMatrix(
        [[1.35, 0.21, -5.23], [0.07, -0.21, 0.6], [2.15, 2.22, -0.84],
         [-2.64, -1.89, 5.47]])
    assert pwm.max_raw_score == pytest.approx(9.84)
    assert pwm.max_raw_score == pytest.approx(9.84)
    assert pwm.min_raw_score == pytest.approx(-9.76)
    assert pwm.min_raw_score == pytest.approx(-9.76)


def test_pwm_score():
    pwm = PositionWeightMatrix(
        [[1.35, 0.21, -5.23], [0.07, -0.21, 0.6], [2.15, 2.22, -0.84],
         [-2.64, -1.89, 5.47]])
    with pytest.raises(ValueError):
        pwm.score("")
    with pytest.raises(ValueError):
        pwm.score("NNNN")
    assert pwm.score("NNN") == 0
    assert pwm.score("AGT") == pytest.approx(0.918699)
    assert pwm.score("ANT") == pytest.approx(0.693089)
    assert pwm.score("CTA") == pytest.approx(-0.716463)

import pytest

from motifscan.motif.score import motif_score, sliding_motif_score


def test_motif_score():
    matrix = [[1.35, 0.21, -5.23], [0.07, -0.21, 0.6], [2.15, 2.22, -0.84],
              [-2.64, -1.89, 5.47]]
    max_raw_score = 1
    scores = motif_score(matrix, 3, max_raw_score,
                         ['NNN', 'AGT', 'ANT', 'CTA'])
    assert scores == pytest.approx([0, 9.04, 6.82, -7.05])


def test_sliding_motif_score():
    matrix = [[1.35, 0.21, -5.23], [0.07, -0.21, 0.6], [2.15, 2.22, -0.84],
              [-2.64, -1.89, 5.47]]
    max_raw_score = 1
    sliding_scores = sliding_motif_score(matrix, 3, max_raw_score,
                                         ['NNNAG', 'TANTCTA'])
    assert sliding_scores[0] == pytest.approx([0, -5.23, -0.63])
    assert sliding_scores[1] == pytest.approx(
        [-2.43, 6.82, -1.29, 2.62, -7.05])

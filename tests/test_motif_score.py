import pytest

from motifscan.motif.cscore import c_score, c_scan_motif


def test_c_score():
    matrix = [[[1.35, 0.21, -5.23], [0.07, -0.21, 0.6], [2.15, 2.22, -0.84],
               [-2.64, -1.89, 5.47]]]
    scores = c_score(matrix, ['NNN', 'AGT', 'ANT', 'CTA'], 1, 1)
    assert len(scores) == 1
    assert scores[0] == pytest.approx(
        [0.0, 0.9186991869918698, 0.693089430894309, -0.7164634146341464])
    scores = c_score(matrix, ['NNN', 'AGT', 'ANT', 'CTA'], 2, 1)
    assert len(scores) == 1
    assert scores[0] == pytest.approx(
        [0.0, 0.6717479674796748, 0.693089430894309, -0.3323170731707317])
    scores = c_score(matrix, ['NNN', 'AGT', 'ANT', 'CTA'], 3, 1)
    assert len(scores) == 1
    assert scores[0] == pytest.approx(
        [0.0, 0.9186991869918698, 0.693089430894309, -0.3323170731707317])


def test_c_scan_motif():
    matrix = [[[1.35, 0.21, -5.23], [0.07, -0.21, 0.6], [2.15, 2.22, -0.84],
               [-2.64, -1.89, 5.47]]]
    motif_sites = c_scan_motif(matrix, [0.2], ['NNNAG', 'TANTCTA'], 3, 1)
    assert len(motif_sites) == 1
    assert len(motif_sites[0]) == 4
    assert motif_sites[0][0] == pytest.approx([1, 1, 0.693089430894309, 1])
    assert motif_sites[0][1] == pytest.approx([1, 1, 0.693089430894309, 2])
    assert motif_sites[0][2] == pytest.approx([1, 2, 0.23983739837398374, 2])
    assert motif_sites[0][3] == pytest.approx([1, 3, 0.266260162601626, 1])

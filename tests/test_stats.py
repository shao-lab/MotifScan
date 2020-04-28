import numpy as np
import pytest

from motifscan.stats import motif_enrichment


def test_stats_motif_enrichment():
    # n_motifs = 2, n_total_input = 3, n_total_control = 5
    motif_sites = [[[True], [True, True], []],
                   [[], [], []],
                   [[], [], [True]]]
    motif_sites_control = [[[True], [], [True, True, True], [], [True]],
                           [[], [True], [], [True], []],
                           [[], [], [], [], []]]
    enrichment_results = motif_enrichment(motif_sites, motif_sites_control)
    assert len(enrichment_results) == 3
    assert enrichment_results[0].n_input == 2
    assert enrichment_results[0].n_control == 3
    assert enrichment_results[0].fold_change == pytest.approx(2 * 5 / 3 / 3)
    assert enrichment_results[1].n_input == 0
    assert enrichment_results[1].n_control == 2
    assert enrichment_results[1].fold_change == 0
    assert enrichment_results[2].n_input == 1
    assert enrichment_results[2].n_control == 0

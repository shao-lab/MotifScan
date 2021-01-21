"""
motifscan.stats
---------------

Statistics used in MotifScan.
"""

from collections import namedtuple

import numpy as np
import scipy.stats as stats

EnrichmentResult = namedtuple('EnrichmentResult',
                              ['name', 'n_input', 'n_control', 'fold_change',
                               'p_enriched', 'p_depleted', 'p_corrected'])


def motif_enrichment(pwms, motif_sites, motif_sites_control):
    """Perform the motif enrichment analysis between input and control regions.
    """
    enrichment_results = []
    n_motifs = len(motif_sites)
    for pwm, sites, sites_control in zip(pwms, motif_sites,
                                         motif_sites_control):
        # number of total input/control regions
        n_input_total = len(sites)
        n_control_total = len(sites_control)
        # number of input/control regions with >= 1 motif site
        n_input = sum([len(sites_by_region) > 0 for sites_by_region in sites])
        n_control = sum(
            [len(sites_by_region) > 0 for sites_by_region in sites_control])
        if (n_input_total > 0) and (n_control > 0):
            fold_change = n_input * n_control_total / n_control / n_input_total
        else:
            fold_change = np.nan
        table = [[n_input, n_input_total - n_input],
                 [n_control, n_control_total - n_control]]
        odds_ratio, p_enriched = stats.fisher_exact(table, 'greater')
        odds_ratio, p_depleted = stats.fisher_exact(table, 'less')
        p_corrected = min(min(p_enriched, p_depleted) * n_motifs, 1)
        enrichment_results.append(
            EnrichmentResult(pwm.matrix_id + ',' + pwm.name, n_input,
                             n_control, fold_change, p_enriched, p_depleted,
                             p_corrected))
    return enrichment_results

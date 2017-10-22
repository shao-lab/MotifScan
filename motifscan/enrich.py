import re
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from motifscan.stats import fisher_exact_custom


def target_enrichment(peak_table, rnd_table, motif_table):
    """Perform the enrichment analysis between sample and random.

    Args:
        peak_table: pandas dataframe, motifscan result table on sample
        rnd_table: pandas dataframe, motifscan result table on random
        motif_table: pandas dataframe, motif information table

    Returns:
        motif_table: pandas dataframe, table containing both motif information and
                                       fisher exact tests statistics
    """
    n_motif = len(motif_table)
    n_peak = len(peak_table)
    n_rand = len(rnd_table)
    n_samp = int(n_rand/n_peak)

    fold_change = np.zeros(n_motif)
    enriched_pvalue = np.zeros(n_motif)
    depleted_pvalue = np.zeros(n_motif)
    oddsratio = np.ones(n_motif)
    pvalue_corrected = np.ones(n_motif)

    peak_tarnum = np.zeros(n_motif)
    rand_tarnum = np.zeros(n_motif)
    peak_tarnum_table = peak_table[
        pd.Index([i for i in peak_table.columns if re.search(r'\.number', i)])]
    rnd_tarnum_table = rnd_table[
        pd.Index([i for i in rnd_table.columns if re.search(r'\.number', i)])]

    for mti, motif_name in zip(range(n_motif), motif_table['name']):
        peak_tarnum[mti] = len(
                [i for i in peak_tarnum_table['%s.number' % motif_name] if i > 0])
        rand_tarnum[mti] = len(
                [i for i in rnd_tarnum_table['%s.number' % motif_name] if i > 0])

        if peak_tarnum[mti] != 0 and rand_tarnum[mti] != 0:
            fold_change[mti] = float(peak_tarnum[mti] * n_rand) / (rand_tarnum[mti] * n_peak)
        else:
            fold_change[mti] = 'NaN'
        table = [[peak_tarnum[mti], n_peak - peak_tarnum[mti]],
                 [rand_tarnum[mti], n_rand - rand_tarnum[mti]]]
        oddsratio[mti], enriched_pvalue[mti] = fisher_exact(table, 'greater')
        oddsratio[mti], depleted_pvalue[mti] = fisher_exact(table, 'less')
        pvalue_corrected[mti] = min(min(depleted_pvalue[mti],
                                        enriched_pvalue[mti]) * n_motif, 1)
    motif_table['target_number'] = peak_tarnum
    motif_table['random_target_number'] = rand_tarnum
    motif_table['fold_change'] = fold_change
    motif_table['enriched_pvalue'] = enriched_pvalue
    motif_table['depleted_pvalue'] = depleted_pvalue
    motif_table['pvalue_corrected'] = pvalue_corrected
    return motif_table


def target_enrichment_peak2peak(res1_tarnum, res2_tarnum):
    merged_table = pd.merge(res1_tarnum, res2_tarnum, left_index=True, right_index=True)
    oddsratio = merged_table.apply(fisher_exact_custom('oddsratio'), 1)
    enriched_pvalue = merged_table.apply(fisher_exact_custom('enrich_pvalue'), 1)
    depleted_pvalue = merged_table.apply(fisher_exact_custom('deplete_pvalue'), 1)
    pvalue_corrected = [min(min(i) * len(merged_table), 1) for i in zip(enriched_pvalue, depleted_pvalue)]
    fold_change = merged_table.apply(lambda r: (r[0] * (r[6] + r[7])) / (r[6] * (r[0] + r[1])), 1)

    merged_table['fold_change'] = fold_change
    merged_table['enriched_pvalue'] = enriched_pvalue
    merged_table['depleted_pvalue'] = depleted_pvalue
    merged_table['pvalue_corrected'] = pvalue_corrected
    merged_table.sort('enrich_pvalue', inplace=True)
    return merged_table




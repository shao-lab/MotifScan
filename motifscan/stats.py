from scipy import stats


def fisher_exact_custom(v):
    def test(r):
        table = [[r[0], r[1]],
                 [r[6], r[7]]]
        oddsratio, enrich_pvalue = stats.fisher_exact(table, 'greater')
        oddsratio, deplete_pvalue = stats.fisher_exact(table, 'less')
        if v == 'oddsratio':
            return oddsratio
        elif v == 'enrich_pvalue':
            return enrich_pvalue
        elif v == 'deplete_pvalue':
            return deplete_pvalue

    return test

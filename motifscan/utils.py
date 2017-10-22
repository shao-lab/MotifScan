import re
import pandas as pd


def extract_motif_name_from_peak_result_table(peak_result_table):
    """
    extract the motif name from the peak_result_table
    """
    motif_name = []
    for col_name in peak_result_table.columns:
        if re.search(r'.*\.number', col_name):
            motif_name.append(col_name[:-7])
    return motif_name


def summarize_result(peaks, motifs, output_dir):
    """ """
    # motifscan result output
    peak_result = {'chr': peaks['chr'], 'start': peaks['start'], 'end': peaks['end'],
                   'seq_start': peaks['seq_start'], 'seq_end': peaks['seq_end'],
                   'summit': peaks['summit'] - peaks['start']}
    target_number_cols = ['chr', 'start', 'end', 'summit']
    motif_score_cols = ['chr', 'start', 'end', 'summit']
    if 'target_gene' in peaks.columns:
        peak_result['target_gene'] = peaks['target_gene']
        peak_result['target_dis'] = peaks['target_dis']
        target_number_cols.extend(['target_gene', 'target_dis'])
        motif_score_cols.extend(['target_gene', 'target_dis'])

    if 'value' in peaks.columns:
        peak_result['value'] = peaks['value']
        target_number_cols.append('value')
        motif_score_cols.append('value')
    for idx, tmp_motif in motifs.iterrows():
        name = tmp_motif['name']
        tmp_df = pd.read_pickle("{0}/{1}".format(output_dir, tmp_motif['id']))
        peak_result['{}.number'.format(name)] = tmp_df["{}.number".format(name)]
        peak_result['{}.ratio'.format(name)] = tmp_df["{}.ratio".format(name)]
        peak_result['{}.taridx'.format(name)] = tmp_df["{}.tarsite".format(name)]
        peak_result['{}.tarsite'.format(name)] = tmp_df["{}.tarsite".format(name)]
        peak_result['{}.tarratio'.format(name)] = tmp_df["{}.tarratio".format(name)]
        target_number_cols.append("{}.number".format(name))
        motif_score_cols.append("{}.ratio".format(name))
    peak_result = pd.DataFrame(peak_result)
    return peak_result, target_number_cols, motif_score_cols

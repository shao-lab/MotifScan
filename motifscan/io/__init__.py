"""
motifscan.io
------------
I/O functions for MotifScan.
"""

import os

from motifscan.io.utils import replace_special_char


def write_sites_table(output_dir, pwms, regions, motif_sites):
    """Write the motif sites number/score summary table."""
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    path_num = os.path.join(output_dir, 'motif_sites_number.xls')
    path_score = os.path.join(output_dir, 'motif_sites_score.xls')
    with open(path_num, 'w') as f_num, open(path_score, 'w') as f_score:
        name_fields = '\t'.join(
            [pwm.matrix_id + ',' + pwm.name for pwm in pwms])
        f_num.write(f"chr\tstart\tend\t{name_fields}\n")
        f_score.write(f"chr\tstart\tend\t{name_fields}\n")
        for idx, region in enumerate(regions):
            n_sites = []
            scores = []
            for sites in motif_sites:  # sites_by_motif
                num = len(sites[idx])  # sites_by_region
                n_sites.append(num)
                if num == 0:
                    scores.append('NA')
                else:
                    scores.append(max([site.score for site in sites[idx]]))
            num_fields = '\t'.join(map(str, n_sites))
            score_fields = '\t'.join(map(str, scores))
            f_num.write(f"{region.chrom}\t{region.start + 1}\t{region.end}\t"
                        f"{num_fields}\n")
            f_score.write(f"{region.chrom}\t{region.start + 1}\t{region.end}\t"
                          f"{score_fields}\n")


def write_sites_bed(output_dir, pwms, regions, motif_sites):
    """Write motif sites bed file(s)."""
    output_dir = os.path.join(output_dir, 'motif_sites')
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    for pwm, sites in zip(pwms, motif_sites):
        name = replace_special_char(pwm.matrix_id + '_' + pwm.name)
        path = os.path.join(output_dir, f"{name}_sites.bed")
        with open(path, 'w') as f_out:
            for idx, region in enumerate(regions):
                for site in sites[idx]:
                    f_out.write(f"{region.chrom}\t{site.start}\t"
                                f"{site.start + pwm.length}\t.\t"
                                f"{site.score}\t{site.strand}\n")


def write_enrich_table(output_dir, enrichment_results):
    """Write the motif enrichment results table."""
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    path = os.path.join(output_dir, 'motif_enrichment.xls')
    # sort by enriched p value and then fold change if p value ties
    enrichment_results.sort(key=lambda x: (x.p_enriched, -x.fold_change))
    with open(path, 'w') as f_out:
        f_out.write("Motif\tNum_input_regions\tNum_control_regions\t"
                    "Fold_change\tEnriched_P_value\tDepleted_P_value\t"
                    "Corrected_P_value\n")
        for res in enrichment_results:
            f_out.write(f"{res.name}\t{res.n_input}\t{res.n_control}\t"
                        f"{res.fold_change}\t{res.p_enriched}\t"
                        f"{res.p_depleted}\t{res.p_corrected}\n")

import logging
import os

import matplotlib as mpl
import numpy as np

mpl.use('Agg')

import matplotlib.pyplot as plt

from motifscan.io.utils import replace_special_char

logger = logging.getLogger(__name__)


def have_same_region_length(regions):
    length = None
    for region in regions:
        if not length:
            length = region.end - region.start
        else:
            if region.end - region.start != length:
                return False
    return True


def have_value_attr(regions):
    for region in regions:
        if region.score is None:
            return False
    return True


def smooth(x, window_len=11):
    if len(x) <= window_len:
        return x
    s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    w = np.hanning(window_len)
    y = np.convolve(w / w.sum(), s, mode='same')
    return y[window_len - 1:-window_len + 1]


def plot_motif_sites_dist(output_dir, regions, pwms, motif_sites,
                          window_size):
    if window_size <= 0:
        if len(regions) == 0:
            logger.error("No regions found for plotting")
            return
        if not have_same_region_length(regions):
            logger.error("Unable to plot when the scanning length is "
                         "different across regions")
            return

    output_dir = os.path.join(output_dir, 'plots')
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if window_size <= 0:
        window_size = regions[0].end - regions[0].start
    extend = window_size // 2

    for pwm, sites in zip(pwms, motif_sites):
        logger.debug(f"Plotting for {pwm.matrix_id + ',' + pwm.name}")
        distances = []
        for idx, region in enumerate(regions):
            for site in sites[idx]:
                distances.append(site.start + pwm.length / 2 - region.summit)
        bin_edges = np.arange(-extend - 5, extend + 6, 10)
        freq, _ = np.histogram(distances, bins=bin_edges)
        if len(distances) > 0:
            freq = smooth(freq / len(distances))
        x = []
        for i in range(len(freq)):
            x.append((bin_edges[i] + bin_edges[i + 1]) // 2)
        fig = plt.figure(figsize=(4, 3.5))
        ax = fig.gca()
        ax.bar(x, freq, width=10, color='#4169E1',
               label=pwm.matrix_id + ',' + pwm.name)
        ax.legend(loc='upper right', fontsize=8, frameon=False)
        ax.set_xlabel('Distance to Center/Summit', fontsize=8)
        ax.set_ylabel('Fraction', fontsize=8)
        ax.set_xlim(-extend - 5, extend + 5)
        if len(distances) > 0:
            ax.set_ylim(0, 1.2 * max(freq))
        else:
            ax.set_ylim(0, 0.1)
        ax.tick_params(axis='both', which='major', labelsize=8)
        fig.subplots_adjust(left=0.15, right=0.98, bottom=0.15, top=0.95)
        name = replace_special_char(pwm.matrix_id + '_' + pwm.name)
        path = os.path.join(output_dir, f"{name}_sites_distributions.pdf")
        fig.savefig(path)
        plt.close()


def plot_motif_sites_enrich(output_dir, regions, pwms, motif_sites,
                            motif_sites_control):
    if not have_value_attr(regions):
        logger.error("Unable to plot when some regions have no scores set for "
                     "sorting")
        return
    n_regions_input = len(regions)
    if len(str(n_regions_input)) < 2:
        logger.error(
            f"Too few regions to plot: {n_regions_input}")
        return

    output_dir = os.path.join(output_dir, 'plots')
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # sort regions by score
    ranked_idx = sorted(range(n_regions_input),
                        key=lambda x: regions[x].score, reverse=True)
    flanking_size = n_regions_input // 100

    for pwm, sites_input, sites_control in zip(pwms, motif_sites,
                                               motif_sites_control):
        logger.debug(f"Plotting for {pwm.matrix_id + ',' + pwm.name}")
        n_regions_control = len(sites_control)
        n_control = sum([len(sites) > 0 for sites in sites_control])
        ratio_control = n_control / n_regions_control
        if ratio_control == 0:
            ratio_control = 1

        sites_input = np.asarray(sites_input)[ranked_idx]
        fold_changes = []
        has_site_flag = [len(sites) > 0 for sites in sites_input]
        for idx in range(n_regions_input):
            head = max(0, idx - flanking_size)
            tail = min(idx + flanking_size, n_regions_input)
            ratio_input = sum(has_site_flag[head:tail]) / (tail - head)
            fold_changes.append(ratio_input / ratio_control)
        fold_changes = smooth(fold_changes)

        fig = plt.figure(figsize=(4, 3.5))
        ax = fig.gca()
        ax.bar(range(1, n_regions_input + 1), fold_changes, width=1,
               color='#4169E1', label=pwm.matrix_id + ',' + pwm.name)
        ax.legend(loc='upper right', fontsize=8, frameon=False)
        ax.set_xlabel('Regions Ranked by Score (Descending)', fontsize=8)
        ax.set_ylabel('Fold Change', fontsize=8)
        ax.set_xlim(0, n_regions_input)
        y_max = max(fold_changes)
        if y_max > 0:
            ax.set_ylim(0, 1.2 * y_max)
        else:
            ax.set_ylim(0, 0.1)
        ax.tick_params(axis='both', which='major', labelsize=8)
        fig.subplots_adjust(left=0.15, right=0.98, bottom=0.15, top=0.95)
        name = replace_special_char(pwm.matrix_id + '_' + pwm.name)
        path = os.path.join(output_dir, f"{name}_sites_enrichment.pdf")
        fig.savefig(path)
        plt.close()

import re
import sys
import math
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from motifscan.utils import extract_motif_name_from_peak_result_table


def _get_highest_digit(n):
    if n / 100000000:
        n /= 100000000
    if n / 10000:
        n /= 10000
    if n / 100:
        n /= 100
    if n / 10:
        n /= 10
    return n


def smooth(x, window_len=5, window='hanning'):
    if isinstance(x, list):
        x = np.array(x)
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len < 3:
        return x
    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett',  blackman'"

    s = np.r_[2 * x[0] - x[window_len:1:-1], x, 2 * x[-1] - x[-1:-window_len:-1]]
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='same')
    # return the smoothed signal,  chopping off the ends so that it has the previous size.
    return y[window_len - 1:-window_len + 1]


def target_site_distribution(peak_df, motif_df, plot_dir, region_radius):
    # win_size = region_radius * 2

    win_size = region_radius * 2
    half_win_size = region_radius
    bin_size = 10

    bin_center = np.arange(-half_win_size, win_size - half_win_size + bin_size, bin_size)
    bin_edge = bin_center - round(bin_size / 2)
    bin_edge = np.append(bin_edge, bin_center[-1] + round(bin_size / 2))
    bin_edge = map(float, bin_edge)
    motif_num = len(motif_df)
    cnt = 0
    for motif_idx, motif_record in motif_df.iterrows():
        cnt += 1
        motif_len = np.shape(motif_record['matrix'])[1]
        motif_name = motif_record['name']
        motif_tarsite = []
        for tar_idx, i in enumerate(peak_df['%s.tarsite' % motif_name]):
            if len(i) > 0:
                for j in i:
                    motif_tarsite.append(j + motif_len / 2 - region_radius)

        motif_tarsite_freq = np.histogram(motif_tarsite, bin_edge)
        motif_marker = motif_tarsite_freq[0]
        motif_marker = motif_marker / float(sum(motif_marker))
        n_motif_marker = len(motif_marker)
        # circular smooth
        motif_marker_long = np.concatenate([motif_marker, motif_marker, motif_marker])
        motif_marker_long = smooth(motif_marker_long, 20)
        motif_marker = motif_marker_long[n_motif_marker:(2 * n_motif_marker)]

        plt.cla()
        # plt.plot(bin_edge[:-1], motif_marker, lw = 2, color='#4169E1', label=motif_name)
        plt.bar(bin_edge[:-1], motif_marker, width=bin_size - float(bin_size) / 10 + float(bin_size) / 20,
                color='#4169E1', linewidth=0, label=motif_name)
        plt.legend()
        ax = plt.gca()
        ax.set_xlabel('Distance to Peak Summit', weight='bold')
        ax.set_ylabel('Fraction', weight='bold')
        ax.set_xlim([min(bin_edge), max(bin_edge)])
        ax.set_ylim([0, 1.5 * max(motif_marker)])
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontweight('bold')
        ax.tick_params(axis='x', which='both', top='off')
        ax.tick_params(axis='y', which='both', right='off')
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(2)
        plt.savefig('%s/%s_target_site.png' % (plot_dir, re.sub(r'[:\-.]{1,2}', '_', motif_name)), dpi=600)
        plt.close()

    return


def tarnum_and_tarsite_distribution(peak_table, rand_table, motif_table, plot_out_dir, bin_size=10, region_radius=500):
    # fold change plot parameter
    motif_name_list = extract_motif_name_from_peak_result_table(peak_table)
    npeak = len(peak_table)
    nrand = len(rand_table)
    nsamp = round(nrand / npeak)
    bin = 1000
    half_bin = bin / 2
    peak_table.sort_values(by='value', ascending=False, inplace=True)
    n_motif = len(motif_name_list)

    # target site paramater
    win_size = region_radius * 2
    half_win_size = win_size / 2

    bin_center = np.arange(-half_win_size, win_size - half_win_size + 1, bin_size)
    bin_edge = bin_center - round(bin_size / 2)
    bin_edge = np.append(bin_edge, bin_center[-1] + round(bin_size / 2))
    bin_edge = map(int, bin_edge)
    n_motif = len(motif_table)

    cnt = 0
    for midx, motif_record in motif_table.iterrows():
        cnt += 1
        motif_len = np.shape(motif_record['matrix'])[1]
        motif_name = motif_record['name']
        # fold change data preparation
        tarnum_fc_smooth = np.zeros(npeak)
        tarnum_smooth = np.zeros(npeak)
        for pi in np.arange(npeak):
            peak_start_idx = max(0, pi - half_bin)
            peak_end_idx = min(npeak, pi + half_bin)
            peak_tarnum = peak_table['%s.number' %
                                     motif_name].iloc[peak_start_idx:peak_end_idx]

            peak_tarnum_smooth = float(
                len(peak_tarnum[peak_tarnum > 0])) / len(peak_tarnum)
            rnd_idx = peak_table.iloc[peak_start_idx:peak_end_idx].index * int(nsamp)
            rand_tarnum = rand_table['%s.number' %
                                     motif_name].ix[rnd_idx]
            rand_tarnum_smooth = float(
                len(rand_tarnum[rand_tarnum > 0])) / len(rand_tarnum)

            tarnum_smooth[pi] = peak_tarnum_smooth
            if rand_tarnum_smooth == 0:
                tarnum_fc_smooth[pi] = peak_tarnum_smooth
            else:
                tarnum_fc_smooth[pi] = peak_tarnum_smooth / rand_tarnum_smooth

        tarnum_fc_smooth = pd.Series(tarnum_fc_smooth)
        tarnum_smooth = pd.Series(tarnum_smooth)

        # target site data preparation
        motif_tarsite = []
        for idx, i in enumerate(peak_table['%s.tarsite' % motif_name]):
            if len(i) > 0:
                for j in i:
                    motif_tarsite.append(j + motif_len / 2 - half_win_size)
        motif_tarsite_freq = np.histogram(motif_tarsite, bin_edge)
        motif_marker = motif_tarsite_freq[0]
        motif_marker = motif_marker / float(sum(motif_marker))
        n_motif_marker = len(motif_marker)
        motif_marker_long = np.concatenate([motif_marker, motif_marker, motif_marker])
        motif_marker_long = smooth(motif_marker_long, 20)
        motif_marker = motif_marker_long[n_motif_marker:(2 * n_motif_marker)]

        # plot
        plt.cla()
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        # fc
        fig.suptitle(motif_name, weight='black')
        ax1.bar(np.arange(npeak), tarnum_fc_smooth, 1,
                color='#4169E1', lw=0, label='fc of %s' % motif_name)
        xtick_step = _get_highest_digit(npeak / 6) * 10 ** int(math.log10(npeak / 6))
        ax1.set_xticks(np.arange(0, npeak, xtick_step))
        ax1.set_xlim(xmin=1, xmax=npeak)
        ax1.set_ylim(ymin=0, ymax=1.3 * max(tarnum_fc_smooth))
        for axis in ['top', 'bottom', 'left', 'right']:
            ax1.spines[axis].set_linewidth(2)
        ax1.tick_params(axis='x', which='both', top='off')
        ax1.tick_params(axis='y', which='both', right='off')
        ax1.set_ylabel('Fold Change', weight='bold')
        ax1.set_xlabel('Peak Rank (Sorted by Value in Descending Order)', weight='bold')
        for tick in ax1.xaxis.get_major_ticks():
            tick.label.set_fontweight('bold')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontweight('bold')

        # targetsite
        ax2.bar(bin_edge[:-1], motif_marker, width=bin_size - 1.5,
                color='#4169E1', linewidth=0, label=motif_name)
        ax2.set_xlabel('Distance to Peak Summit', weight='bold')
        ax2.set_ylabel('Fraction', weight='bold')
        ax2.set_xlim([min(bin_edge), max(bin_edge)])
        ax2.set_ylim([0, 1.5 * max(motif_marker)])
        for tick in ax2.xaxis.get_major_ticks():
            tick.label.set_fontweight('bold')
        for tick in ax2.yaxis.get_major_ticks():
            tick.label.set_fontweight('bold')
        ax2.tick_params(axis='x', which='both', top='off')
        ax2.tick_params(axis='y', which='both', right='off')
        for axis in ['top', 'bottom', 'left', 'right']:
            ax2.spines[axis].set_linewidth(2)
        fig.savefig('%s/%s_%s_tarsite_fc_dist.png' % (plot_out_dir, cnt, re.sub(r'[:\-.]{1,2}', '_', motif_name)),
                    dpi=600)
        plt.close()
    return

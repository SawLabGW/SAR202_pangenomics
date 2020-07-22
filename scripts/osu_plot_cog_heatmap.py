#!/usr/bin/env python
__author__ = "Jimmy Saw"

"""

Script used to generate Figure 2A

usage:
cd ~/cgrb_workspace/SAR202/refined_final_set/rpsblast/results

osu_plot_cog_heatmap.py \
    -c top50_normalized_COGs_counts_high_precision_bin_sorted.txt \
    -g ../../bins/new_groups_050918.txt -s Purples \
    -d ../../bins/depths.txt -o SAR202_cog_abundances_zscore.pdf \
    -a 0.8 -p custom -z yes -m euclidean -t median

osu_plot_cog_heatmap.py \
    -c top50_normalized_COGs_counts_high_precision_bin_sorted.txt \
    -g ../../bins/new_groups_050918.txt -s Purples \
    -d ../../bins/depths.txt -o SAR202_cog_abundances_pabund.pdf \
    -a 0.8 -p custom -z no -m euclidean -t single

## best combinations of column clustering for metric (-m) and method (-t) (with z_score=0):
## 1. euclidean, median
## 2. sqeuclidean, complete
## 3. braycurtis, average
## 4. euclidean, single

## best combinations of column clustering for metric (-m) and method (-t) (without z_score):
## 1. euclidean, single
## 2. euclidean, weighted
## 3. cityblock, average

## final plot for the manuscript
cd ~/cgrb_workspace/SAR202/refined_final_set/rpsblast/results
osu_plot_cog_heatmap.py \
    -c top50_normalized_COGs_counts_high_precision_bin_sorted.txt \
    -g ../../bins/new_groups_081018.txt \
    -s Greys \
    -d ../../bins/depths.txt \
    -a 0.8 -p custom -z yes -m euclidean -t median \
    -o SAR202_cog_abundances_zscore.pdf

"""

import sys
import argparse
import operator
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

def plot_heatmap(cogs, group_mapping, palette, outfile, alph, shade, depth, zscore, metric, method):
    """
    This function plots seaborn heatmap
    :param cogs:
    :param group_mapping:
    :param palette:
    :param outfile:
    :return: None
    """
    group_dict = {}
    with open(group_mapping) as gf:
        lines = gf.readlines()
        for line in lines:
            c = line.split("\t")
            group_dict[c[0]] = c[1].strip()
    unique_groups = list(set([v for k, v in group_dict.iteritems()]))
    group_colors = {}

    if palette == "custom":
        x = {'Group1a': '#EED2EE', 'Group1b': '#CDB6CD', 'Group1c': '#CD96CD', 'Group2': '#607B8B',
             'Group3a': '#EECFA1', 'Group3b': '#CDB38B', 'Group3c': '#8B795E', 'Group4': '#8FBC8F',
             'Group5': '#FA8072', 'Group6': '#8DB6CD'}
        for i in unique_groups:
            group_colors[i] = x[i]
    else:
        gcolors = sns.color_palette(palette, len(unique_groups))
        for i, j in enumerate(unique_groups):
            group_colors[j] = gcolors[i]

    df = pd.read_csv(cogs, sep="\t", header=0)
    df = df.set_index('bins')
    colors = []
    for i in df.index:
        if i in group_dict:
            colors.append(group_colors[group_dict[i]])

    depth_dict = {}
    depth_colors = {}
    depth_bins = ['A: 0-200m', 'B: 200-1000m', 'C: 1000-4000m', 'D: 4000-6000m', 'E: 6000-12000m', 'F: unknown']

    if palette == "custom":
        dcolors = sns.color_palette("ocean_r", len(depth_bins))
    else:
        dcolors = sns.color_palette(palette, len(depth_bins))
    for a, b in zip(depth_bins, dcolors):
        depth_colors[a] = b
    depth_colors['F: unknown'] = "#CDC1C5"
    depth_df = pd.read_csv(depth, sep="\t", header=None, names=['orig_id','bins','depths'])
    for a, b in zip(depth_df['bins'], depth_df['depths']):
        bcat = 'F: unknown'
        if 0 <= int(b) <= 200:
            #epipelagic
            bcat = 'A: 0-200m'
        elif 200 < int(b) <= 1000:
            #mesopelagic
            bcat = 'B: 200-1000m'
        elif 1000 < int(b) <= 4000:
            #bathypelagic
            bcat = 'C: 1000-4000m'
        elif 4000 < int(b) <= 6000:
            #abyssopelagic
            bcat = 'D: 4000-6000m'
        elif 6000 < int(b) <= 12000:
            #hadopelagic
            bcat = 'E: 6000-12000m'
        else:
            bcat = 'F: unknown'
        depth_dict[a] = bcat
    dd_colors = []
    for i in df.index:
        bin_category = 'F: unknown'
        if i in depth_dict:
            bin_category = depth_dict[i]
        dd_colors.append(depth_colors[bin_category])

    sns.set_style('white')
    #sns.clustermap(df, z_score=1, method='ward', cmap='OrRd', row_colors=colors)
    #grid = sns.clustermap(df, cmap='OrRd', row_colors=colors)

    alpha = float(alph)
    df.to_csv("data_frame.tab", sep="\t") #save the dataframe

    if zscore in ['Yes', 'yes']:
        grid = sns.clustermap(df, figsize=(10, 8), cmap=shade, row_cluster=False, row_colors=[colors, dd_colors],
                              alpha=alpha, xticklabels=True, yticklabels=False, z_score=0, metric=metric,
                              method=method, cbar_kws={'label': 'z score'})
    else:
        grid = sns.clustermap(df, figsize=(10, 8), cmap=shade, row_cluster=False, row_colors=[colors, dd_colors],
                              alpha=alpha, xticklabels=True, yticklabels=False, metric=metric,
                              method=method, cbar_kws={'label': '% abundance'})
    grid.ax_col_dendrogram.set_visible(False)
    grid.cax.set_position([0.15, 0.40, 0.02, 0.15])

    sorted_groups = sorted(group_colors.items(), key=operator.itemgetter(0))
    cs = []
    ls = []
    for i in sorted_groups:
        cs.append(i[1])
        ls.append(i[0])
    gpatches = [mpatches.Patch(color=c, label=l, alpha=alpha) for c, l in zip(cs, ls)]

    sorted_depths = sorted(depth_colors.items(), key=operator.itemgetter(0))
    dcs = []
    dls = []
    for i in sorted_depths:
        dcs.append(i[1])
        dls.append(i[0])
    dpatches = [mpatches.Patch(color=c, label=l, alpha=alpha) for c, l in zip(dcs, dls)]

    plt.setp(grid.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(grid.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.setp(grid.ax_heatmap.yaxis.get_majorticklabels(), fontsize=8)
    plt.setp(grid.ax_heatmap.xaxis.get_majorticklabels(), fontsize=8)
    grid.fig.legend(handles=gpatches, loc='center left', ncol=1, prop={'size': 9})
    grid.fig.legend(handles=dpatches, loc='lower left', ncol=1, prop={'size': 9})


    grid.savefig(outfile)
    #grid.savefig(outfile, dpi=600)

    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script plots heatmap of COG categories")
    parser.add_argument("-c", "--cog_counts", required=True, help="COG counts file")
    parser.add_argument("-g", "--group_mapping", required=True, help="group mapping file")
    parser.add_argument("-d", "--depth", required=True, help="Depth information for SAR202 bins")
    parser.add_argument("-p", "--palette", required=True, help="color palette to use")
    parser.add_argument("-a", "--alpha", required=True, help="Alpha value for both heatmap and legends: 0-1")
    parser.add_argument("-s", "--shade", required=True, help="Shade color of heatmap: use palette colors")
    parser.add_argument("-z", "--zscore", required=True, help="To plot using zscore or not: yes or no")
    parser.add_argument("-o", "--outfile", required=True, help="outfile to save")
    parser.add_argument("-m", "--metric", required=True, help="Clustering metric")
    parser.add_argument("-t", "--method", required=True, help="Clustering method")
    args = parser.parse_args()

    plot_heatmap(args.cog_counts, args.group_mapping, args.palette, args.outfile,
                 args.alpha, args.shade, args.depth, args.zscore, args.metric, args.method)


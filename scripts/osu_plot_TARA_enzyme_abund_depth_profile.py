#!/usr/bin/env python
__author__ = "Jimmy Saw"

"""
Script used to plot Figure 7A, 7B, and S2

usage example:
cd ~/cgrb_workspace/TARAOCEANS/enzyme_abundances/SAR202

osu_plot_TARA_enzyme_abund_depth_profile.py \
    -a sorted_SAR202_fmno_counts.txt \
    -t ../../OM.CompanionTables.txt \
    -e FMNO -o SAR202_fmnos_B109_abundances_depth

osu_plot_TARA_enzyme_abund_depth_profile.py \
    -a sorted_SAR202_enol_counts.txt \
    -t ../../OM.CompanionTables.txt \
    -e Enolases -o SAR202_enols_B109_abundances_depth
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import operator
import re
from pandas import Series
import matplotlib.patches as mpatches

def plot_ars_red(comp_table, abund, enz, out, dataset):
    """
    Plots relative abundances of a given enzyme in question
    :param comp_table:
    :param abund:
    :param enz:
    :return:
    """
    sns.set_style('white')
    comp_df = pd.read_csv(comp_table, sep="\t", header=0)
    unique_biomes = Series.unique(
        comp_df['Ocean and sea regions (IHO General Sea Areas 1953) [MRGID registered at www.marineregions.com]'])
    biome_colors = {}

    custom = ['#333333', '#CD96CD', '#607B8B', '#CDBE70', '#EE7621', '#4F94CD', '#8B8B83', '#548B54']

    colors = sns.color_palette(custom)

    for i, b in enumerate(sorted(unique_biomes)):
        bn = re.sub("\[.*", "", b)
        biome_colors[bn] = colors[i]
    abund_df = pd.read_csv(abund, sep="\t", header=0)

    comp_dict = {}
    for a, b, c, d, e in zip(comp_df['Sample label [TARA_station#_environmental-feature_size-fraction]'],
                          comp_df['Latitude [degrees North]'], comp_df['Longitude [degrees East]'],
                          comp_df['Sampling depth [m]'],
                          comp_df['Ocean and sea regions (IHO General Sea Areas 1953) [MRGID registered at www.marineregions.com]']):
        comp_dict[a] = (b, c, d, e)

    counts = []

    for i, j, k in zip(abund_df['bins'], abund_df[enz], abund_df['genes']):
        if i in comp_dict:
            depth =  comp_dict[i][2]
            percent = 0
            if not k == 0:
                percent = (j / float(k)) * 100
            #print i, percent
            biome = re.sub("\[.*", "", comp_dict[i][3])
            color = biome_colors[biome]
            counts.append([i, depth, percent, biome, color])

    #M = 20
    df2 = pd.DataFrame(data=dict(x=[n for n, v in enumerate(counts)],
                                 y=[i[1] for i in counts],
                                 a2=[i[2] for i in counts],
                                 names=[i[0] for i in counts],
                                 colors=[biome_colors[i[3]] for i in counts]
                                 ))
    #bins = np.linspace(df2.a2.min(), df2.a2.max(), M)
    print "min", df2.a2.min()
    print "maximum relative abundance =", df2.a2.max()
    #print "bins", bins

    size_factor = 10
    if df2.a2.max() >= 10:
        size_factor = 10
    elif df2.a2.max() >= 1:
        size_factor = 100
    elif df2.a2.max() >= 0.5:
        size_factor = 500
    elif df2.a2.max() >= 0.1:
        size_factor = 1000
    elif df2.a2.max() >= 0.05:
        size_factor = 1500
    elif df2.a2.max() >= 0.01:
        size_factor = 4000
    else:
        size_factor = 10000

    fig, ax = plt.subplots(1, 1, figsize=(14, 7))

    for index, j in enumerate(sorted(unique_biomes)):
        x = []
        y = []
        bc = []
        s = []
        for k, v in enumerate(counts):
            if v[3] == re.sub("\[.*", "", j):
                x.append(k)
                y.append(v[1])
                bc.append(biome_colors[v[3]])
                s.append(v[2] * size_factor)
        ax.scatter(x, y, s=s, c=bc, alpha=0.5)

    ## draw legend patches for biomes
    sorted_biomes = sorted(biome_colors.items(), key=operator.itemgetter(0))
    cs = []
    ls = []
    for i in sorted_biomes:
        cs.append(i[1])
        ls.append(i[0])
    patches = [mpatches.Patch(color=c, label=l, alpha=0.8) for c, l in zip(cs, ls)]
    leg = plt.legend(handles=patches, labels=[i for i in ls], bbox_to_anchor=(1.25, 1), ncol=1, prop={'size': 10})
    ax.add_artist(leg)

    ax.grid(True)

    ## scatter point legend
    slabels = []
    spoints = []
    points = np.linspace(0, df2.a2.max(), 6)
    for i in points:
        x = plt.scatter([], [], s=i * size_factor, c="#104E8B", alpha=0.5)
        spoints.append(x)
        slabels.append('{0:.3f}'.format(i))
    plt.legend(handles=spoints[1:], labels=slabels[1:], ncol=1, fontsize=10, bbox_to_anchor=(1.11, 0.5),
               prop={'size': 10}, labelspacing=1)

    plt.xlabel("Samples")
    plt.ylabel("depth (m)")
    plt.title(dataset + " " + enz + " abundances in the TARA Oceans metagenomes")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(out + ".pdf")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script plots SAR202 16S abundances in EBI TARA assemblies")
    parser.add_argument("-a", "--abund", required=True, help="Abundances file from miTAG data")
    parser.add_argument("-t", "--table", required=True, help="Companion table file")
    parser.add_argument("-e", "--enz", required=True, help="Enzyme to check")
    parser.add_argument("-o", "--out", required=True, help="Out prefix")
    parser.add_argument("-s", "--set", required=True, help="Choose: SAR202 or TARA")
    args = parser.parse_args()

    plot_ars_red(args.table, args.abund, args.enz, args.out, args.set)

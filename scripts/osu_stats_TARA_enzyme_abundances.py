#!/usr/bin/env python

__author__ = "Jimmy Saw"

"""
Script used to make Figure 7C

Usage example:
cd ~/cgrb_workspace/TARAOCEANS/enzyme_abundances/SAR202
osu_stats_TARA_enzyme_abundances.py \
    -t sorted_SAR202_fmno_counts.txt \
    -c FMNO -s pearsonr \
    -o TARA_SAR202_fmnos_B109_abundances_cor_pearson.pdf

osu_stats_TARA_enzyme_abundances.py \
    -t sorted_SAR202_fmno_counts.txt \
    -c FMNO -s spearmanr \
    -o TARA_SAR202_fmnos_B109_abundances_cor_spearman.pdf


"""

import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr, kendalltau, weightedtau

def plot(table, output, enzyme, stat):
    """
    Plots enzyme abundances
    :param table:
    :param output:
    :param enzyme:
    :return: plot
    """
    df = pd.read_csv(table, sep="\t", header=0)
    bins = df['bins']
    biomes = df['biome']
    percents = []
    for i, j in zip(df[enzyme], df['genes']):
        if j != 0:
            percents.append((float(i) / float(j))*100)
        else:
            percents.append(0)
    depths = [float(i) for i in df['depth']]
    newdf = pd.DataFrame({'bin': bins, 'biome': biomes, 'depth': depths, 'abundance': percents})
    sns.set(style="ticks", color_codes=True)
    if stat == "pearsonr":
        g = sns.jointplot("abundance", "depth", data=newdf, kind="reg", color="#000000", size=6, robust=True,
                          scatter_kws={"s": 10, "color": "#228B22"}, stat_func=pearsonr)
    elif stat == "spearmanr":
        g = sns.jointplot("abundance", "depth", data=newdf, kind="reg", color="#000000", size=6, robust=True,
                          scatter_kws={"s": 10, "color": "#228B22"}, stat_func=spearmanr)
    elif stat == "kendalltau":
        g = sns.jointplot("abundance", "depth", data=newdf, kind="reg", color="#000000", size=6, robust=True,
                          scatter_kws={"s": 10, "color": "#228B22"}, stat_func=kendalltau)
    else:
        print "Choose either pearsonr, spearmanr, or kendalltau for option -s"
        sys.exit()

    g.fig.axes[0].invert_yaxis()
    g.fig.suptitle("% abundance of " + enzyme + " vs. depth", size=12)
    plt.xlabel("normalized % abundance")
    plt.ylabel("depth (m)")
    plt.subplots_adjust(top=0.90)
    plt.savefig(output, format='pdf', dpi=1000)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script plots statistical plots for: fmnos vs depth, enolases vs depth, etc")
    parser.add_argument("-t", "--table", required=True, help="Table of COG counts")
    parser.add_argument("-c", "--cog", required=True, help="COG to check")
    parser.add_argument("-s", "--stats", required=True, help="Choose: pearsonr, spearmanr, or kendalltau")
    parser.add_argument("-o", "--outfile", required=True, help="outfile prefix for saving output")
    args = parser.parse_args()
    plot(args.table, args.outfile, args.cog, args.stats)


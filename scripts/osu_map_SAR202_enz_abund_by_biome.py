#!/usr/bin/env python

__author__ = "Jimmy Saw"

"""
Script used to make Figures 7A and 7B

usage example:

cd ~/cgrb_workspace/SAR202/refined_final_set/bins
osu_map_SAR202_enz_abund_by_biome.py -p moll -t TableS1.txt -e rhds -r l -o RHD_bin_distribution
osu_map_SAR202_enz_abund_by_biome.py -p moll -t TableS1.txt -e fmnos -r l -o FMNO_bin_distribution
"""

import argparse
import operator
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.basemap import Basemap
import matplotlib.patches as mpatches

def plot_distribution(table, enz, proj, res, out):
    df = pd.read_csv(table, sep="\t", header=0)
    bins = df['bin']
    lats = df['lat']
    lons = df['lon']
    depths = df['depth(m)']
    genes = df['genes']
    ecounts = df[enz]

    m = Basemap(projection=proj, lon_0=0, lat_0=0, resolution=res)
    pcts = []

    fig = plt.figure(figsize=(12, 6))

    factor = 50

    for bin, lon, lat, d, g, ec in zip(bins, lons, lats, depths, genes, ecounts):
        if lon.startswith("mixed") or lat.startswith("mixed"):
            pass
        else:
            p = (float(ec) / float(g)) * 100
            pcts.append(p)
            x, y = m(float(lon), float(lat))
            m.scatter(x, y, marker='o', c='#FF6103', s=p * factor, alpha=0.5, zorder=3)

    m.drawcoastlines(linewidth=0.5, color='#9E9E9E', antialiased=1)
    m.drawlsmask(land_color='#CCCCCC', ocean_color='#B0E2FF')

    slabels = []
    spoints = []

    points = np.linspace(0, max(pcts), 6)
    for i in points[1:]:
        x = plt.scatter([], [], s=i * factor, c="#FF6103",
                        alpha=0.8)  # multiple by size factor to get similar size as actual data points
        spoints.append(x)
        slabels.append('{0:.2f}'.format(i))
    leg1 = plt.legend(spoints, slabels, ncol=5, frameon=False, fontsize=8, handlelength=2, loc=8,
                      bbox_to_anchor=(0.5, -0.07))
    ax1 = plt.gca().add_artist(leg1)

    enzymes = {'fmnos': 'FMNO', 'enolases': 'Enolase', 'rhds': 'RHD', 'dehydro': 'Dehydrogenases'}

    plt.title("SAR202 SAG/MAG isolation sources and " + enzymes[enz] + " abundances at all depths")
    plt.savefig(out + "_map.pdf", bbox_extra_artists=(leg1,), bbox_inches='tight')
    #plt.show()

def plot_depth_profile(table, enz, out):
    sns.set_style('white')
    df = pd.read_csv(table, sep="\t", header=0)
    bins = df['bin']
    #lats = df['lat']
    #lons = df['lon']
    depths = df['depth(m)']
    genes = df['genes']
    ecounts = df[enz]
    groups = df['subgroup']

    group_colors = {'Group1a': '#EED2EE', 'Group1b': '#CDB6CD', 'Group1c': '#CD96CD', 'Group2': '#607B8B',
             'Group3a': '#EECFA1', 'Group3b': '#CDB38B', 'Group3c': '#8B795E', 'Group4': '#8FBC8F',
             'Group5': '#FA8072', 'Group6': '#8DB6CD'}

    enzymes = {'fmnos': 'FMNO', 'enolases': 'Enolase', 'rhds': 'RHD', 'dehydro': 'Dehydrogenases'}

    counts = []

    for i, j, k, l, m in zip(bins, depths, genes, ecounts, groups):
        pct = (float(l) / float(k)) * 100
        color = group_colors[m]
        if not j.startswith("mixed"):
            counts.append([i, int(j), pct, m, color])

    sorted_counts = sorted(counts, key=operator.itemgetter(3))

    df2 = pd.DataFrame(data=dict(x=[n for n, v in enumerate(sorted_counts)],
                                 y=[i[1] for i in sorted_counts],
                                 a2=[i[2] for i in sorted_counts],
                                 names=[i[3] for i in sorted_counts],
                                 colors=[i[4] for i in sorted_counts]
                                 ))
    print "min", df2.a2.min()
    print "maximum relative abundance =", df2.a2.max()
    print "maximum depth = ", df2.y.max()

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

    for index, j in enumerate(sorted(group_colors.keys())):
        x = []
        y = []
        bc = []
        s = []
        for k, v in enumerate(sorted_counts):
            if j == v[3]:
                x.append(k)
                y.append(v[1])
                bc.append(v[4])
                s.append(v[2] * size_factor)
        ax.scatter(x, y, s=s, c=bc, alpha=0.8)

    ## draw legend patches for groups
    sorted_groups = sorted(group_colors.items(), key=operator.itemgetter(0))
    cs = []
    ls = []
    for i in sorted_groups:
        cs.append(i[1])
        ls.append(i[0])
    patches = [mpatches.Patch(color=c, label=l, alpha=0.8) for c, l in zip(cs, ls)]
    leg = plt.legend(handles=patches, labels=[i for i in ls], bbox_to_anchor=(1.13, 1), ncol=1, prop={'size': 10})
    ax.add_artist(leg)

    ax.grid(True)

    ## scatter point legend
    slabels = []
    spoints = []
    points = np.linspace(0, df2.a2.max(), 6)
    for i in points:
        x = plt.scatter([], [], s=i * size_factor, c="white", edgecolors="#000000", alpha=0.8)
        spoints.append(x)
        slabels.append('{0:.3f}'.format(i))
    plt.legend(handles=spoints[1:], labels=slabels[1:], ncol=1, fontsize=10, bbox_to_anchor=(1.11, 0.5),
               prop={'size': 10}, labelspacing=1)

    plt.xlabel("bins")
    plt.ylabel("depth (m)")
    plt.ylim(-100, 2500)
    #plt.xticks(np.arange(0, len(xtlabels) + 1, 1.0))
    plt.title("Normalized SAR202 " + enzymes[enz] + " abundances in each genome")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(out + "_profile.pdf")
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script plots locations of SAR202 enzyme abundances on the "
                                                 "World Map.")
    parser.add_argument("-p", "--projection", required=True, help="Choose projection types (examples: kav7, gall, eck4,"
                                                                  "robin, mbtfpq, cea)")
    parser.add_argument("-t", "--table", required=True, help="Table containing all data")
    parser.add_argument("-e", "--enz", required=True, help="Enzyme to check")
    parser.add_argument("-r", "--res", required=True, help="Resolution for map boundaries: c, l, i, h, f")
    parser.add_argument("-o", "--outprefix", required=True, help="Out prefix")
    args = parser.parse_args()

    plot_distribution(args.table, args.enz, args.projection, args.res, args.outprefix)
    plot_depth_profile(args.table, args.enz, args.outprefix)

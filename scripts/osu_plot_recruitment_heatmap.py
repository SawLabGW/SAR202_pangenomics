#!/usr/bin/env python

__author__ = "Jimmy Saw"

"""
Script used to make Figures 6, S5, S6

example usage: 
cd ~/cgrb_workspace/SAR202/refined_final_set/recruitment/trenches
osu_plot_recruitment_heatmap.py \
    -f 3samples_recr_bp_sorted_by_depth.txt \
    -b bins.txt -c ../bats/fr-hit/groups.txt \
    -s Purples -g OGT/fr-hit/genome_sizes.txt \
    -d trench_depths.txt -o trenches_recruitment_heatmap \
    -y correlation \
    -t weighted \
    -z [yes or no] \
    -l [yes or no] \
    -p [tab20|Paired_r|Accent|custom] \
    -m [trenches|bats|tara]

Example use cases:
cd ~/cgrb_workspace/SAR202/refined_final_set/recruitment/trenches

osu_plot_recruitment_heatmap.py -f 3samples_recr_bp_sorted_by_depth_and_newgroup.txt -b bins.txt -c new_groups_050918.txt -s Purples -g genome_sizes.txt -d trench_depths.txt -p custom -m trenches -o test -z yes -y cityblock -t weighted -l no -x no
osu_plot_recruitment_heatmap.py -f 3samples_recr_bp_sorted_by_depth_and_newgroup_062718.txt -b bins.txt -c new_groups_062718.txt -s Purples -g genome_sizes.txt -d trench_depths.txt -p custom -m trenches -o test -z yes -y cityblock -t weighted -l no -x yes

cd ~/cgrb_workspace/SAR202/refined_final_set/recruitment/tara/fr-hit
osu_plot_recruitment_heatmap.py -f all_tara_bps_sorted_depths_groups.txt -b bins.txt -c new_groups_050918.txt -s Purples -g genome_sizes.txt -d sorted_samples_depths.txt -o test -p custom -m tara -z yes -y cityblock -t weighted -l yes
osu_plot_recruitment_heatmap.py -f all_tara_bps_sorted_depths_groups.txt -b bins.txt -c new_groups_050918.txt -s Purples -g genome_sizes.txt -d sorted_samples_depths.txt -o test -p custom -m tara -z yes -y cityblock -t weighted -l no

cd ~/cgrb_workspace/SAR202/refined_final_set/recruitment/bats/fr-hit
osu_plot_recruitment_heatmap.py -f bats_frhit_bp_counts_sorted_depths_new_groups.txt -b bins.txt -c new_groups_050918.txt -s Purples -g genome_sizes.txt -d depths.txt -o test -p custom -m bats -z yes -y cityblock -t weighted -l yes
osu_plot_recruitment_heatmap.py -f bats_frhit_bp_counts_sorted_depths_new_groups.txt -b bins.txt -c new_groups_050918.txt -s Purples -g genome_sizes.txt -d depths.txt -o test -p custom -m bats -z yes -y cityblock -t weighted -l no

## best options for clustering so far:
1. -l no (clustering using known subgroups)
2. -l yes -y cityblock -t weighted (de novo clustering)
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import operator
import seaborn as sns
from matplotlib.patches import Rectangle

def plot_heatmap(table, groups, depths, shade, palette, out, metaset, zscore, metric, method, cluster, label):
    ## parse count table
    df = pd.read_csv(table, sep="\t", header=0)
    ndf = df.set_index('bins')
    tdf = ndf.T

    ## parse depth info
    ddf = pd.read_csv(depths, sep=" ", header=None, names=['sample','source','depth'])
    depth_dict = {}
    for a, b, c in zip(ddf['sample'], ddf['source'], ddf['depth']):
        depth_dict[a] = (b, c)

    sns.set_style('white')
    param = {'xtick.labelsize': 6, 'ytick.labelsize': 6}
    sns.set(rc = param)
    group_dict = {}
    with open(groups) as gf:
        lines = gf.readlines()
        for line in lines:
            c = line.split("\t")
            group_dict[c[0]] = c[1].strip()
    unique_groups = list(set([v for k, v in group_dict.iteritems()]))
    group_colors = {}
    #gcolors = []
    if palette == "custom":
        ## make custom group colors
        gcolors = {'Group1a': '#EED2EE', 'Group1b': '#CDB6CD', 'Group1c': '#CD96CD', 'Group2': '#607B8B',
                   'Group3a': '#EECFA1', 'Group3b': '#CDB38B', 'Group3c': '#8B795E', 'Group4': '#8FBC8F', 
                   'Group5': '#FA8072', 'Group6': '#8DB6CD', 'Group7': '#4682B4'}
        for i in unique_groups:
            group_colors[i] = gcolors[i]
    else:
        gcolors = sns.color_palette(palette, len(unique_groups))
        for i, j in enumerate(unique_groups):
            group_colors[j] = gcolors[i]

    colors = []
    for i in tdf.columns:
        if i in group_dict:
            colors.append(group_colors[group_dict[i]])

    renamed = []
    row_colors = []
    samples = {'MRT': 'Mariana Trench', 'OGT': 'Ogasawara Trench', 'JPT': 'Japan Trench', 'BATS': 'BATS', 
        'IO': 'Indian Ocean', 'MS': 'Mediterranean Sea', 'NAO': 'North Atlantic Ocean', 
        'NPO': 'North Pacific Ocean',  'RS': 'Red Sea', 'SAO': 'South Atlantic Ocean', 
        'SO': 'Southern Ocean', 'SPO': 'South Pacific Ocean'}
    rc = {'MRT': '#8B8B83', 'OGT': '#CDCDC1', 'JPT': '#EEEEE0', 'BATS':'#CDCDC1',
        'IO': '#CDAF95', 'MS': '#CD96CD', 'NAO': '#6CA6CD',
        'NPO': '#B4EEB4',  'RS': '#FF7256', 'SAO': '#7EC0EE', 
        'SO': '#C1CDC1', 'SPO': '#9BCD9B'}
    for i in tdf.index:
        if i in depth_dict:
            x = ''
            y = depth_dict[i][0]
            if y in samples:
                x = samples[y]
            if y in rc:
                row_colors.append(rc[y])
            #renamed.append(x + " (" + str(depth_dict[i][1]) + ")")
            renamed.append(str(depth_dict[i][1]))
    tdf.index = renamed

    figsize = (14, 10)
    tlabels = False
    if label == "yes":
        tlabels = True
    if zscore in ['Yes', 'yes']:
        if cluster == "yes":
            grid = sns.clustermap(tdf, figsize=figsize, cmap=shade, col_cluster=True, row_colors=row_colors,
                                  col_colors=colors, row_cluster=False, alpha=0.85, xticklabels=tlabels, yticklabels=True,
                                  metric=metric, method=method, z_score=0, cbar_kws={'label': 'z score'})
            grid.cax.set_position([0.2, 0.40, 0.02, 0.15])
        else:
            grid = sns.clustermap(tdf, figsize=figsize, cmap=shade, col_cluster=False, row_colors=row_colors,
                                  col_colors=colors, row_cluster=False, alpha=0.85, xticklabels=tlabels, yticklabels=True,
                                  metric=metric, method=method, z_score=0, cbar_kws={'label': 'z score'})
            grid.cax.set_position([0.2, 0.40, 0.02, 0.15])
    else:
        if cluster == "yes":
            grid = sns.clustermap(tdf, figsize=figsize, cmap=shade, col_cluster=True, row_colors=row_colors,
                                  col_colors=colors, row_cluster=False, alpha=0.95, xticklabels=tlabels, yticklabels=True,
                                  metric=metric, method=method, cbar_kws={'label': 'recruitment'})
            grid.cax.set_position([0.2, 0.40, 0.02, 0.15])
        else:
            grid = sns.clustermap(tdf, figsize=figsize, cmap=shade, col_cluster=False, row_colors=row_colors,
                                  col_colors=colors, row_cluster=False, alpha=0.95, xticklabels=tlabels, yticklabels=True,
                                  metric=metric, method=method, cbar_kws={'label': 'recruitment'})
            grid.cax.set_position([0.2, 0.40, 0.02, 0.15])

    ax = grid.ax_heatmap
    if metaset == "trenches":
        ax.add_patch(Rectangle((0, 0), 122, 6, fill=False, edgecolor='grey', lw=1, alpha=0.5))
        ax.add_patch(Rectangle((0, 6), 122, 9, fill=False, edgecolor='grey', lw=1, alpha=0.5))
        ax.add_patch(Rectangle((0, 15), 122, 7, fill=False, edgecolor='grey', lw=1, alpha=0.5))
    elif metaset == "bats":
        ax.add_patch(Rectangle((0, 0), 122, 8, fill=False, edgecolor='grey', lw=1, alpha=0.5))
        ax.add_patch(Rectangle((0, 8), 122, 4, fill=False, edgecolor='grey', lw=1, alpha=0.5))
        ax.add_patch(Rectangle((0, 12), 122, 5, fill=False, edgecolor='grey', lw=1, alpha=0.5))
    elif metaset == "tara":
        ax.add_patch(Rectangle((0, 0), 122, 9, fill=False, edgecolor='grey', lw=1, alpha=0.5))
        ax.add_patch(Rectangle((0, 9), 122, 2, fill=False, edgecolor='grey', lw=1, alpha=0.5))
        ax.add_patch(Rectangle((0, 11), 122, 3, fill=False, edgecolor='grey', lw=1, alpha=0.5))
        ax.add_patch(Rectangle((0, 14), 122, 5, fill=False, edgecolor='grey', lw=1, alpha=0.5))
        ax.add_patch(Rectangle((0, 19), 122, 8, fill=False, edgecolor='grey', lw=1, alpha=0.5))
        ax.add_patch(Rectangle((0, 27), 122, 5, fill=False, edgecolor='grey', lw=1, alpha=0.5))
        ax.add_patch(Rectangle((0, 32), 122, 10, fill=False, edgecolor='grey', lw=1, alpha=0.5))
        ax.add_patch(Rectangle((0, 42), 122, 1, fill=False, edgecolor='grey', lw=1, alpha=0.5))
    else:
        print "No boxes drawn on the heatmap" 

    sorted_groups = sorted(group_colors.items(), key=operator.itemgetter(0))

    ## draw legend patches for groups
    cs = []
    ls = []
    for i in sorted_groups:
        cs.append(i[1])
        ls.append(i[0])
    patches = [mpatches.Patch(color=c, label=l, alpha=0.8) for c, l in zip(cs, ls)]

    #grid.fig.legend(handles=patches, labels=[i for i in ls], loc='lower left', ncol=1, prop={'size': 8})
    grid.fig.legend(handles=patches, labels=[i for i in ls], loc='upper right', ncol=1, prop={'size': 8})

    #grid.fig.legend(handles=patches, labels=[i for i in ls], loc='center left', ncol=1, prop={'size': 9} )
    #plt.legend(bbox_to_anchor=(0, 0), loc='upper left', handles=patches, prop={'size': 9}) 

    ## draw legend patches for samples
    sdict = {'bats': ['BATS'], 'tara': ['IO', 'RS', 'MS', 'NAO', 'SAO', 'NPO', 'SPO', 'SO'], 'trenches': ['JPT', 'OGT', 'MRT']}
    scs = []
    sls = []
    ## draw patches only for samples in a given metagenomic dataset
    targets = sdict[metaset]
    newdict = {}
    for t in targets:
        newdict[t] = rc[t]
        scs.append(rc[t])
        sls.append(samples[t])
    sample_patches = [mpatches.Patch(color=c, label=l, alpha=0.8) for c, l in zip(scs, sls)]
    #grid.fig.legend(handles=sample_patches, labels=[i for i in sls], loc='center left', ncol=1, prop={'size': 8})
    grid.fig.legend(handles=sample_patches, labels=[i for i in sls], loc='upper center', ncol=1, prop={'size': 8})
    
    plt.setp(grid.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.setp(grid.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    #grid.fig.subplots_adjust(top=0.95)

    grid.savefig(out + ".pdf")
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script calculates fragment recruitment")
    parser.add_argument("-f", "--frhit", required=True, help="FR-Hit output file summary")
    parser.add_argument("-b", "--bins", required=True, help="list of shortened bin names")
    parser.add_argument("-c", "--classification", required=True, help="group classes")
    parser.add_argument("-s", "--shade", required=True, help="eg: Blues, Reds")
    parser.add_argument("-p", "--palette", required=True, help="eg: terrain, Set3, etc")
    parser.add_argument("-d", "--depths", required=True, help="Trench depths")
    parser.add_argument("-g", "--genome_sizes", required=True, help="Genome sizes")
    parser.add_argument("-m", "--meta", required=True, help="Name of metagenomic dataset. Choose either: trenches, bats, or tara")
    parser.add_argument("-z", "--zscore", required=True, help="Whether to calculate z-score or not")
    parser.add_argument("-y", "--metric", required=True, help="Metric for calculation of distance matrix: correlation, euclidean, braycurtis, canberra, etc")
    parser.add_argument("-t", "--method", required=True, help="Method for linkage clustering: weighted, average, centroid, single, complete, etc")
    parser.add_argument("-l", "--cluster", required=True, help="Whether to cluster columns or not: if 'yes' then will ignore -y and -t options")
    parser.add_argument("-o", "--outfile", required=True, help="outfile prefix for saving output")
    parser.add_argument("-x", "--label", required=True, help="Label the x ticks (bin names): yes or no")
    args = parser.parse_args()
    plot_heatmap(args.frhit, args.classification, args.depths, args.shade, args.palette, args.outfile, 
        args.meta, args.zscore, args.metric, args.method, args.cluster, args.label)

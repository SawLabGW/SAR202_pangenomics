#!/usr/bin/env python
__author__ = "Jimmy Saw"

"""
usage example:

cd ~/cgrb_workspace/SAR202/refined_final_set/phylogenies/phylogenomics/chloNOGs/analyses/single_gene_trees
osu_SAR202_single_gene_trees.py \
    -t tree.contree \
    -m group.mappings \
    -c group.colors \
    -a ../fasta/chloNOG44_annotations.txt \
    -n mono.out \
    -o test.pdf

for i in `ls *.contree`;do
    osu_SAR202_single_gene_trees.py \
    -t $i -m group.mappings -c group.colors \
    -a ../fasta/chloNOG44_annotations.txt -o $i.pdf \
    -n $i.monophyly
done

cd ~/cgrb_workspace/SAR202/refined_final_set/phylogenies/phylogenomics/protein_ortho108/single_gene_trees
for i in `ls *.contree`;do
    osu_SAR202_single_gene_trees.py \
    -t $i -m group.mappings \
    -c group.colors -a ortholog_annotations.txt \
    -n $i.monophyly -o $i.pdf
done
"""

import operator
import argparse
import pandas as pd
import ete3
from ete3 import Tree, TreeStyle, NodeStyle, faces, TextFace

def tree_layout(node):
    """
    side function for tree layouts (following ete3 manual)
    :param node:
    :return:
    """
    if node.is_leaf():
        color = get_color(node.name)
        N = ete3.faces.TextFace(node.name, fsize=10, fgcolor=color)
        faces.add_face_to_node(N, node, 0)

def get_color(node_name):
    """
    Color the leaves according to known groups
    :param node_name:
    :return: color hex code
    """
    cdf = pd.read_csv(color_mapping, sep="\t", header=None, names=['group', 'color'])
    color_dict = dict(zip(cdf['group'], cdf['color']))
    gdf = pd.read_csv(group_mapping, sep="\t", header=None, names=['taxon', 'group'])
    group_dict = dict(zip(gdf['taxon'], gdf['group']))
    color = "black"
    if node_name in group_dict:
        group = group_dict[node_name]
        if group in color_dict:
            color = color_dict[group]
    return color

def consider_other_outgroups(members_found):
    """
    Consider alternative outgroups if cyanobacteria are NOT monophyletic
    :param members_found:
    :return: alternative outgroup members
    """
    members_gte2 = []
    for k, v in members_found.iteritems():
        if not k == "Cyanobacteria":
            if len(members_found[k]) >= 1:
                members_gte2.append(k)
    if len(members_gte2) > 0:
        return members_gte2[0], members_found[members_gte2[0]], len(members_found[members_gte2[0]])
    else:
        return "None", "None", 0

def determine_outgroup(tree, group_dict):
    """
    function to determine outgroup to use
    :param tree:
    :param group_dict:
    :return: outgroup to use, members in the group, and count
    """
    og_to_check = ['Cyanobacteria', 'Anaerolineae', 'Caldilineae', 'Thermoflexia', 'Chloroflexia', 'Dehalococcoidia',
                   'Ktedonobacteria']
    members_found = {}
    for i in og_to_check:
        members_found[i] = []
    for og in og_to_check:
        for n in tree.traverse():
            if n.name in group_dict[og]:
                members_found[og].append(n.name)

    if len(members_found['Cyanobacteria']) == 1:
        return 'Cyanobacteria', members_found['Cyanobacteria'], len(members_found['Cyanobacteria'])
    elif len(members_found['Cyanobacteria']) == 2:
        r = tree.get_midpoint_outgroup()
        tree.set_outgroup(r)
        x = tree.check_monophyly(values=members_found['Cyanobacteria'], target_attr="name")
        if x[0] == True:
            return 'Cyanobacteria', members_found['Cyanobacteria'], len(members_found['Cyanobacteria'])
        else:
            ## if cyanos are not monophyletic, consider other outgroups
            alternative_outgroup = consider_other_outgroups(members_found)
            return alternative_outgroup
    else:
        alternative_outgroup = consider_other_outgroups(members_found)
        return alternative_outgroup

def check_monophyly(tree, group_dict, mono):
    """
    Check monophyly of all known groups. Also check monophyly of all SAR202 members
    :param tree:
    :param group_dict:
    :param mono:
    :return: writes a text file of monophyly results
    """
    groups = list(set(group_dict.keys()))
    results = {}
    total_mono = 0
    sar202 = []
    dehalo_and_sar202 = []
    for g in groups:
        members = []
        for n in tree.traverse():
            if n.name in group_dict[g]:
                members.append(n.name)
                if g.startswith('Group'):
                    sar202.append(n.name)
                    dehalo_and_sar202.append(n.name)
                if g == "Dehalococcoidia":
                    dehalo_and_sar202.append(n.name)
        if len(members) > 1:
            x = tree.check_monophyly(values=members, target_attr="name")
            if x[0] == True:
                results[g] = "yes"
                total_mono += 1
            else:
                results[g] = "no"
        else:
            results[g] = "--"
    ## check SAR202 monophyly
    y = tree.check_monophyly(values=sar202, target_attr="name")
    if y[0] == True:
        results['all_sar202'] = "yes"
        total_mono += 1
    else:
        results['all_sar202'] = "no"

    ## check Dehalo + SAR202 monophyly
    y = tree.check_monophyly(values=dehalo_and_sar202, target_attr="name")
    if y[0] == True:
        results['dehalo_and_sar202'] = "yes"
        total_mono += 1
    else:
        results['dehalo_and_sar202'] = "no"

    sorted_mono_res = sorted(results.items(), key=operator.itemgetter(0))
    with open(mono, "w") as outfile:
        strn = ""
        for item in sorted_mono_res:
            strn += '{0}\t{1}\n'.format(item[0], item[1])
        strn += '{0}\t{1}\n'.format("Total", str(total_mono))
        outfile.write(strn)

def analyze_tree(tree, outfile, anno, monoout):
    """
    Main function to analyze the tree
    :param tree:
    :param outfile:
    :param anno:
    :param monoout:
    :return: tree file as pdf file
    """
    m = pd.read_csv(group_mapping, sep="\t", header=None, names=['taxon', 'group'])
    taxon_mapping = dict(zip(m['taxon'], m['group']))
    unique_groups = list(set(m['group']))

    group_dict = {}
    for u in unique_groups:
        group_dict[u] = []
    for k, v in taxon_mapping.iteritems():
        if taxon_mapping[k] in group_dict:
            group_dict[taxon_mapping[k]].append(k)

    adf = pd.read_csv(anno, sep="\t", header=None, names=['gene', 'anno'])
    anno_dict = dict(zip(adf['gene'], adf['anno']))
    title = "unknown gene"
    if tree.startswith("chloNOG"):
        gene = '.'.join(tree.split(".")[:3])
        if gene in anno_dict:
            title = gene + " (" + anno_dict[gene] + ")"
    elif tree.startswith("OC_"):
        gene = tree.split(".")[0]
        if gene in anno_dict:
            title = gene + " (" + anno_dict[gene] + ")"

    t = Tree(tree)
    outgroup_to_use = determine_outgroup(t, group_dict)

    print title

    if outgroup_to_use[2] == 0:
        print "None of the outgroups chosen are found in the tree. Rooted with mid-point rooting."
        r = t.get_midpoint_outgroup()
        t.set_outgroup(r)
    elif outgroup_to_use[2] == 1:
        t.set_outgroup(outgroup_to_use[1][0])
    elif outgroup_to_use[2] >= 2:
        ## do mid-point rooting first to get around the problem of some clades that can't be re-rooted right away
        r = t.get_midpoint_outgroup()
        t.set_outgroup(r)
        x = t.check_monophyly(values=outgroup_to_use[1], target_attr="name")
        if x[0] == True:
            print outgroup_to_use[0], "members are monophyletic"
            ## now, re-root with actual outgroup
            lca = t.get_common_ancestor(outgroup_to_use[1])
            t.set_outgroup(lca)
        else:
            ## use mid-point rooting if members are not mono (already mid-point rooted in the outer scope)
            print outgroup_to_use[0], "members are NOT monophyletic. Can't use as an outgroup. Rooted with mid-point rooting"

    check_monophyly(t, group_dict, monoout)

    t.ladderize(direction=1)
    ts = TreeStyle()
    ns = NodeStyle()
    ts.show_branch_support = True
    ts.extra_branch_line_color = "DarkGrey"
    ts.show_leaf_name = False
    ts.layout_fn = tree_layout
    #ts.branch_vertical_margin = 0
    ns['shape'] = "square"
    ns['size'] = 0
    ts.title.add_face(TextFace(title, fsize=8), column=0)
    for n in t.traverse():
        n.set_style(ns)
    t.render(outfile, w=1500, units="px", tree_style=ts)
    print " "

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script re-roots the tree with proper outgroup, checks for "
                                                 "monophyly of groups, and colors leaves based on known association")
    parser.add_argument("-t", "--tree", required=True, help="Tree file")
    parser.add_argument("-m", "--mapping", required=True, help="Mapping file (taxon to clade mapping file)")
    parser.add_argument("-c", "--colors", required=True, help="Color mapping file for groups")
    parser.add_argument("-a", "--anno", required=True, help="Annotations of genes")
    parser.add_argument("-n", "--mono", required=True, help="Output file to write monophyly check results")
    parser.add_argument("-o", "--outfile", required=True, help="Outfile name - PDF")
    args = parser.parse_args()
    global color_mapping
    color_mapping = args.colors
    global group_mapping
    group_mapping = args.mapping
    analyze_tree(args.tree, args.outfile, args.anno, args.mono)


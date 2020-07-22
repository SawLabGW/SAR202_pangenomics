#!/usr/bin/env python

__author__ = "Jimmy Saw"
"""
Usage example:

## compare my TARA bins amongst themselves
cd ~/cgrb_workspace/SAR202/final_set/unique_TARA_bins
osu_uniquefy_TARA_bins.py \
    -b bins.list -f ../all_TARA_bins \
    -t ../../../TARAOCEANS/OM.CompanionTables.txt \
    -m ../nucmer_output/coords.list -c among \
    -p ../all_checkm.txt -d ../nucmer_output > similar_bins_my_TARA.txt

## compare my TARA bins to Tully bins
cd ~/cgrb_workspace/SAR202/final_set/unique_TARA_bins
osu_uniquefy_TARA_bins.py \
    -b bins.list -f ../all_TARA_bins \
    -t ../../../TARAOCEANS/OM.CompanionTables.txt \
    -m ../nucmer_output/coords.list -c between \
    -p ../all_checkm.txt -d ../nucmer_output > similar_bins_my_TARA_vs_Tully.txt

## compare final bins to Delmont bins
cd ~/cgrb_workspace/TARAOCEANS/Delmont_bins/check_duplicates

osu_uniquefy_TARA_bins.py \
    -b bins.list -f all_bins \
    -t ../../OM.CompanionTables.txt \
    -m coords.list \
    -c final -p CheckM.tab \
    -d ANIm_all/ANIm_all/nucmer_output > similar_genomes.txt

"""

import os
import argparse
import pandas as pd
import itertools
import numpy as np
from Bio import SeqIO

def get_bins_list(accession, bins):
    """
    Returns list of bins belonging to a given accession number
    :param accession:
    :param bins:
    :return:
    """
    bins_list = []
    for bin in bins:
        parent_accession = bin.split(".")[0]
        if parent_accession == accession:
            bins_list.append(bin)
    return bins_list

def calculate_ANI(mummer_file_path):
    """
    Calculates Average Nucleotide Identity (ANI) between two bins compared using MUMMer
    :param mummer_file_path:
    :return: (ANI, matchlen1, matchlen2)
    """
    average = 0
    with open(mummer_file_path) as fh:
        lines = fh.readlines()
        total_id = 0
        matchlen1 = 0
        matchlen2 = 0
        counts = 0
        if len(lines) > 4:
            for line in lines[4:]:
                cols = line.split("\t")
                total_id += float(cols[6])
                matchlen1 += int(cols[4])
                matchlen2 += int(cols[5])
                counts += 1
            average = total_id / counts
    return (average, matchlen1, matchlen2)

def compare_permutations(permutations, dir_name, bin_folder, comp):
    print '{0:38}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}'.format('permutation', 'ANIm', 'mlen1',
                                                                                   'mlen2', 'cov1', 'cov2', 'b1size',
                                                                                   'b2size', 'b1comp', 'b2comp',
                                                                                   'cont1', 'cont2')
    completeness = {}
    with open(comp, "rU") as cf:
        lines = cf.readlines()
        for line in lines:
            x = line.split("\t")
            completeness[x[0]] = (float(x[12]), float(x[13].strip()))
    for perm in permutations:
        suffix = '.coords'
        m1_percent = 0
        m2_percent = 0
        cfile = perm + suffix
        fpath = os.path.join(dir_name, cfile)
        if os.path.isfile(fpath):
            ani, m1, m2 = calculate_ANI(fpath)
            bin1 = perm.split("_vs_")[0]
            bin2 = perm.split("_vs_")[1]
            bin1_path = os.path.join(bin_folder, bin1 + ".fa")
            bin2_path = os.path.join(bin_folder, bin2 + ".fa")
            bin1size = 0
            bin2size = 0
            comp1 = 0
            comp2 = 0
            cont1 = 0
            cont2 = 0
            if os.path.isfile(bin1_path):
                bin1size = np.sum([len(i.seq) for i in SeqIO.parse(bin1_path, "fasta")])
                m1_percent = (float(m1) / bin1size) * 100
            if os.path.isfile(bin2_path):
                bin2size = np.sum([len(i.seq) for i in SeqIO.parse(bin2_path, "fasta")])
                m2_percent = (float(m2) / bin2size) * 100
            if ani >= 99.0 and m1 >= 1000 and m2 >= 1000:
                if bin1 in completeness and bin2 in completeness:
                    comp1 = completeness[bin1][0]
                    comp2 = completeness[bin2][0]
                    cont1 = completeness[bin1][1]
                    cont2 = completeness[bin2][1]
                print '{0:38}\t{1:.2f}\t{2}\t{3}\t{4:.2f}\t{5:.2f}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}'.format(perm, ani, m1, m2,
                                                                                               m1_percent, m2_percent,
                                                                                               bin1size, bin2size,
                                                                                               comp1, comp2, cont1, cont2)

def compare_tully_and_my_bins(bins, dir_name, bin_folder, comp):
    bs = []
    with open(bins) as fh:
        lines = fh.readlines()
        for line in lines:
            bs.append(line.strip())
    mybins = []
    tbins = []
    for bin in bs:
        if bin.startswith("ERR"):
            mybins.append(bin)
        elif bin.startswith("TOBG"):
            tbins.append(bin)
    perms = map('_vs_'.join, itertools.chain(itertools.product(mybins, tbins), itertools.product(tbins, mybins)))
    compare_permutations(perms, dir_name, bin_folder, comp)

def compare_final_and_delmont(bins, dir_name, bin_folder, comp):
    bs = []
    with open(bins) as fh:
        lines = fh.readlines()
        for line in lines:
            bs.append(line.strip())
    final_bins = []
    delmont_bins = []
    for bin in bs:
        if bin.startswith("TARA_"):
            delmont_bins.append(bin)
        else:
            final_bins.append(bin)
    perms = map('_vs_'.join, itertools.chain(itertools.product(final_bins, delmont_bins),
                                             itertools.product(delmont_bins, final_bins)))
    compare_permutations(perms, dir_name, bin_folder, comp)

def compare_my_tara_bins(bins, meta, coords, dir_name, bin_folder, comp):
    """
    Compares bins from the same sample to see if they are almost identical.
    :param bins:
    :param meta:
    :param coords:
    :param dir_name:
    :param bin_folder:
    :return: Prints a list of bin permutations with their ANI, match lengths, match coverages, and parent sample name
    """
    completeness = {}
    with open(comp, "rU") as cf:
        lines = cf.readlines()
        for line in lines:
            x = line.split("\t")
            completeness[x[0]] = (float(x[12]), float(x[13].strip()))
    bs = []
    with open(bins) as fh:
        lines = fh.readlines()
        for line in lines:
            bs.append(line.strip())
    meta_df = pd.read_csv(meta, sep="\t", header=0)
    run_accs = meta_df['INSDC run accession number(s)']
    sample_labels = meta_df['Sample label [TARA_station#_environmental-feature_size-fraction]']
    print '{0:38}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:22}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}'.format('permutation', 'ANIm',
                                                     'mlen1', 'mlen2', 'cov1', 'cov2', 'sample', 'b1size', 'b2size',
                                                     'b1comp', 'b2comp', 'cont1', 'cont2')
    for sample, run in zip(sample_labels, run_accs):
        acc_list = [i for i in run.split('|')]
        bins_to_compare = []
        for acc in acc_list:
            y = get_bins_list(acc, bs)
            for x in y:
                bins_to_compare.append(x)
        permutations = itertools.permutations(bins_to_compare, 2)
        for perm in permutations:
            prefix = '_vs_'.join([perm[0], perm[1]])
            suffix = '.coords'
            m1_percent = 0
            m2_percent = 0
            cfile = prefix + suffix
            fpath = os.path.join(dir_name, cfile)
            if os.path.isfile(fpath):
                ani, m1, m2 = calculate_ANI(fpath)
                bin1 = perm[0]
                bin2 = perm[1]
                comp1 = 0
                comp2 = 0
                bin1size = 0
                bin2size = 0
                cont1 = 0
                cont2 = 0
                bin1_path = os.path.join(bin_folder, bin1 + ".fa")
                bin2_path = os.path.join(bin_folder, bin2 + ".fa")
                if os.path.isfile(bin1_path):
                    bin1size = np.sum([len(i.seq) for i in SeqIO.parse(bin1_path, "fasta")])
                    m1_percent = (float(m1) / bin1size) * 100
                if os.path.isfile(bin2_path):
                    bin2size = np.sum([len(i.seq) for i in SeqIO.parse(bin2_path, "fasta")])
                    m2_percent = (float(m2) / bin2size) * 100
                if ani >= 99.0 and m1 >= 1000 and m2 >= 1000:
                    if bin1 in completeness and bin2 in completeness:
                        comp1 = completeness[bin1][0]
                        comp2 = completeness[bin2][0]
                        cont1 = completeness[bin1][1]
                        cont2 = completeness[bin2][1]
                    print '{0:38}\t{1:.2f}\t{2}\t{3}\t{4:.2f}\t{5:.2f}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}'.format(
                                            prefix, ani, m1, m2, m1_percent, m2_percent, sample, bin1size, bin2size,
                                            comp1, comp2, cont1, cont2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script tries to identify redundant bins in TARA metagenomic samples")
    parser.add_argument("-b", "--bins", required=True, help="List of bins to check")
    parser.add_argument("-f", "--folder", required=True, help="")
    parser.add_argument("-t", "--tara_metadata", required=True, help="TARA metadata file")
    parser.add_argument("-m", "--mummer_list", required=True, help="List of MUMMer coord files")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing MUMMer coord files")
    parser.add_argument("-c", "--compare", required=True, help="Choose either: among or between")
    parser.add_argument("-p", "--completeness", required=True, help="Reformatted CheckM Completeness/Contamination file")
    args = parser.parse_args()
    if args.compare == 'among':
        compare_my_tara_bins(args.bins, args.tara_metadata, args.mummer_list, args.directory, args.folder,
                             args.completeness)
    elif args.compare == 'between':
        compare_tully_and_my_bins(args.bins, args.directory, args.folder, args.completeness)
    elif args.compare == 'final':
        compare_final_and_delmont(args.bins, args.directory, args.folder, args.completeness)


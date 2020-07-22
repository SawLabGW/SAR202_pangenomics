#!/usr/bin/env python
__author__ = "Jimmy Saw"

"""
usage example:
cd ~/cgrb_workspace/SAR202/refined_final_set/phylogenies/phylogenomics/markers139
osu_extract_single_copy_genes.py \
    -p pfam_ids.list -m tbls.list \
    -s ../seqs_links \
    -o single_copy_markers.txt \
    -d single_copy_markers

"""

import os
import argparse
import pandas as pd
import re
from Bio import SeqIO

def make_mapping_file(sdf, outdir):
    """
    Makes a mapping file for concatenation
    :param sdf:
    :param outdir:
    :return:
        makes and outdir if not present
        writes a mapping file (org_name, pfam_id, gene_name)
    """
    fname = "mapping.txt"
    outname = os.path.join(outdir, fname)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    with open(outname, "w") as outfile:
        for pf in sdf.columns:
            org_names = list(sdf[pf].index)
            for org, gene in zip(org_names, sdf[pf]):
                if len(gene) == 1:
                    str = '{0}\t{1}\t{2}\n'.format(org, pf, gene[0])
                    outfile.write(str)

def parse_hmmsearch_tbl(pfamlist, hmmlist, out):
    """
    Parses hmmsearch tbl out files and returns a dataframe of rows with counts less than 2
    :param pfamlist:
    :param hmmlist:
    :param out:
    :return: dataframe containing gene ids
    """
    pfams = [i.strip() for i in open(pfamlist).readlines()]
    pcounts = []
    genes_list = []
    bins_list = []

    hmm_results = [i.strip() for i in open(hmmlist).readlines()]
    for file in hmm_results:
        pdict = {}
        gdict = {}
        for p in pfams:
            pdict[p] = 0
            gdict[p] = []
        with open(file) as tbl:
            bin = re.sub(".tbl", "", file)
            bins_list.append(bin)
            lines = tbl.readlines()
            for line in lines:
                if not line.startswith("#"):
                    c = line.split()
                    gene = c[0]
                    pf = c[2]
                    ev = float(c[4])
                    if ev <= 1e-5:
                        pdict[pf] += 1
                        gdict[pf].append(gene)
        pcounts.append(pdict)
        genes_list.append(gdict)

    df = pd.DataFrame(pcounts, index=bins_list)
    genes_df = pd.DataFrame(genes_list, index=bins_list)

    # select only pfams with less than 2 (1 or 0)
    selected = df.ix[:, df.lt(2).all()]  #select all rows where all columns have less than 2
    gselected = genes_df.ix[:, df.lt(2).all()]  #positions of genes are same as counts, so you can use counts indices to extract
    print "original df:", df.shape
    print "selected df:", selected.shape
    #df.to_csv(out, sep="\t")
    selected.to_csv(out, sep="\t")
    #print gselected.head()
    #return selected
    return gselected

def extract_seqs(sdf, sfolder, outdir):
    """
    Extracts protein sequences of each organism grouped by pfam id (those with 0 counts are excluded)
    :param sdf:
    :param sfolder:
    :param outdir:
    :return:
    """
    for pf in sdf.columns:
        pflist = []
        org_names = list(sdf[pf].index)
        for org, gene in zip(org_names, sdf[pf]):
            if len(gene) == 1:
                gid = gene[0]
                seqfasta = os.path.join(sfolder, org + ".faa")
                seqs = [i for i in SeqIO.parse(seqfasta, "fasta")]
                for seq in seqs:
                    if gid == seq.id:
                        pflist.append(seq)
        print pf, len(pflist)
        fname = pf + ".faa"
        outname = os.path.join(outdir, fname)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        SeqIO.write(pflist, outname, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script extracts single-copy marker genes from a given set of genomes")
    parser.add_argument("-p", "--pfam_list", help="List of pfam ids from 139 bact marker genes set")
    parser.add_argument("-m", "--hmm_results", required=True, help="List of hmmsearch results to parse")
    parser.add_argument("-s", "--seqfolder", required=True, help="folder containing sequences to extract from")
    parser.add_argument("-o", "--outfile", required=True, help="Dataframe of pfams and counts saved")
    parser.add_argument("-d", "--dir", required=True, help="output directory to save fasta files")
    args = parser.parse_args()

    seqdf = parse_hmmsearch_tbl(args.pfam_list, args.hmm_results, args.outfile)
    extract_seqs(seqdf, args.seqfolder, args.dir)
    make_mapping_file(seqdf, args.dir)

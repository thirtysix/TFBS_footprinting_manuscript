#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 3.8.0 ###########################################################
# Libraries ####################################################################
from pathlib import Path
import pandas as pd
import fastparquet
import pyarrow
import os
import csv
import time
import json
import numpy as np
import seaborn as sns

from operator import itemgetter
from _collections import defaultdict
from collections import Counter

import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from math import log10,floor

from assign_dirs import *
################################################################################
# Description/Notes ############################################################
################################################################################
"""
OBJECTIVE:
Map all ChIP-Seq peaks cataloged in the GTRD database to the 5' and 3' ends of
the most highly expressed transcript for each gene cataloged in the GTEx
database.

OUTPUT:
Figure depicting summed counts of ChIP-Seq peaks overlapping regions centered on
5' and 3' ends of all most highly expressed transcripts.  Each subplot of the
figure is a composite of ~58,000 transcripts.

TO DO:
Add compression to table outfile/infile.
"""

################################################################################
# Base-level Functions #########################################################
################################################################################


def load_json(object_filename):
    """Load a json file to object."""

    if os.path.exists(object_filename):
        with open(object_filename, 'r') as object_file:
            return json.load(object_file)

def save_json(json_dict, json_filename):
    """Save json object to file."""
    
    with open(json_filename, 'w') as json_file:
        json.dump(json_dict, json_file)
################################################################################
# Task-specific Functions ######################################################
################################################################################
def overlap_range(x,y):
    return range(max(x[0], y[0]), min(x[-1], y[-1]))


def get_gtex_data():
    """Get GTEx transcript level TPM counts (v8; 2017-06-05) ~3.5GB"""

    print("Starting attempt to download {} data.".format('GTEx'))
    static_url = "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"

    # download archive file if it doesn't exist
    if not os.path.exists(gtex_transcripts_expr_fn):
        try:
            ssl._create_default_https_context = ssl._create_unverified_context # avoid 'SSL: CERTIFICATE_VERIFY_FAILED' error
            urllib.request.urlretrieve(static_url, gtex_transcripts_expr_fn)

        except:
            print("Download of {} data was unsuccessful".format('GTEx'))


def get_gtrd_data():
    """Get GTRD ChIP-Seq data.  This file appears to be updated semi-regularly so processing of it may need updating periodically."""

    print("Starting attempt to download {} data.".format('GTRD'))
    static_url = "http://gtrd.biouml.org:8888/downloads/20.06/intervals/chip-seq/Homo_sapiens_ChIP-seq_peaks_MACS2.zip"

    # download archive file if it doesn't exist
    if not os.path.exists(gtrd_peaks_fn):
        try:
            ssl._create_default_https_context = ssl._create_unverified_context # avoid 'SSL: CERTIFICATE_VERIFY_FAILED' error
            urllib.request.urlretrieve(static_url, gtrd_peaks_fn)
        except:
            print("Download of {} data was unsuccessful".format('GTRD'))


def parse_gtf():
    """
    Parse the human genome GTF file for transcript entries.
    """

    print ("parsing GTF")
    transcript_list = []

    # import as df
    ensgtf_df = pd.read_csv(human_gtf_data_filename,
                            compression='gzip',
                            sep="\t",
                            skiprows=5,
                            header=None,
                            names =
                            ['chrom','evidence','feature','start','end','dot1','strand','dot2','annotations'])

    # reduce to transcript rows
    ensgtf_transcript_df = ensgtf_df[ensgtf_df['feature']=='transcript']

    # parse annotations, row by row
    for i, row in ensgtf_transcript_df.iterrows():
        annotations_ls = row.annotations.split("; ")

        # only those rows with annotations are of interest, others lack gene names/ids
        if len(annotations_ls) >0:
            annotations_dict = {x.split(" ")[0]:x.split(" ")[1].replace('"', '') for x in annotations_ls}

            # add to list
            if "gene_id" in annotations_dict and "gene_name" in annotations_dict:
                transcript_data = ["chr"+str(row.chrom), row.strand, row.start, row.end, annotations_dict["gene_name"], annotations_dict["gene_id"], annotations_dict["transcript_id"], annotations_dict["transcript_biotype"].replace(";","")]
                transcript_list.append(transcript_data)

    # convert to dict for quick access
    human_gtf_dict = {x[6]:x for x in transcript_list}

    return human_gtf_dict


def gtrd_get_all_peaks():
    """
    Load GTRD data reduce to target columns (by use_cols),
    filter to entries with fold_enrichment above threshold.
    """

    # define in/out fn
    compiled_filtered_GTRD_peaks_fn = os.path.join(data_GTRD_dn, "compiled_filtered.GTRD_peaks.parquet")

    # compile and create df if parquet file doesn't exist
    if not os.path.exists(compiled_filtered_GTRD_peaks_fn):
        chr_list = ['X', '20', '1', '6', '3', '7', '12', '11', '4', '17', '2', '16', '8', '19', '9', '13', '5', '22', '14', '10', 'Y', '18', '15', '21', 'MT']

        # initiate empty df for concat
        gtrd_reduced_entries_df = pd.DataFrame({}, columns = ["#CHROM", "START", "END", "fold_enrichment"])

        # buffer list for holding parsed dfs
        concat_ls = []

        # iterate through all GTRD files
        for fn_i, gtrd_peaks_fn in enumerate(gtrd_peaks_fns):
            if fn_i%1000 == 0:
                print("GTRD file:", fn_i)    

            # load data, restrict cols used
            gtrd_peaks_df = pd.read_csv(gtrd_peaks_fn, usecols = ["#CHROM", "START", "END", "fold_enrichment"], sep="\t", dtype={"#CHROM":str, "START":int, "END":int, "fold_enrichment":float})

            # filter by fold enrichment
            reduced_entries_df = gtrd_peaks_df[gtrd_peaks_df['fold_enrichment'] >= gtrd_fold_enrichment_cutoff]

            # filter by chromosome
            reduced_entries_df = reduced_entries_df[reduced_entries_df['#CHROM'].isin(chr_list)]

            # add to concat list
            concat_ls.append(reduced_entries_df)

            # periodically concat
            if fn_i%50 == 0:
                concat_ls_df = pd.concat(concat_ls)
                gtrd_reduced_entries_df = pd.concat([gtrd_reduced_entries_df, concat_ls_df])
                concat_ls = []

        # final concat
        concat_ls_df = pd.concat(concat_ls)
        gtrd_reduced_entries_df = pd.concat([gtrd_reduced_entries_df, concat_ls_df])

        # sort entries
        gtrd_reduced_entries_df.sort_values(by=["#CHROM", "START"], inplace=True)

        # save to parquet
        gtrd_reduced_entries_df.to_parquet(compiled_filtered_GTRD_peaks_fn, index=False)

    else:
        gtrd_reduced_entries_df = pd.read_parquet(compiled_filtered_GTRD_peaks_fn)

    return gtrd_reduced_entries_df


def gtex_get_expressed_transcripts():
    """Load GTEX data to identify the single most-expressed transcript for each gene
    rather than testing all transcripts for each gene, many of which have same transcript start site."""
    
    gtex_most_expr_transcripts_fn = os.path.join(data_GTEx_dn, "gtex.most_expr_transcripts.tsv")

    if not os.path.exists(gtex_most_expr_transcripts_fn):

        gtex_iter = pd.read_csv(gtex_transcripts_expr_fn, sep="\t", header=2, iterator=True, chunksize=10**3)
        gtex_sum_df = pd.DataFrame(data={}, columns=["transcript_id", "gene_id", "transcript_expr_sum"])

        chunk_count = 0
        for chunk in gtex_iter:
            chunk_count += 1
            chunk['transcript_expr_sum'] = chunk.sum(axis=1)

            gtex_sum_df = pd.concat([gtex_sum_df, chunk[["transcript_id", "gene_id", "transcript_expr_sum"]]])

            print(chunk_count)
                                          
        # sort values gene_id and transcript_sum
        gtex_sum_sorted_df = gtex_sum_df.sort_values(['gene_id', 'transcript_expr_sum'], ascending=[True, False])
        gtex_most_expr_transcripts_df = gtex_sum_sorted_df.groupby(['gene_id']).head(1)
        gtex_most_expr_transcripts_df.to_csv(gtex_most_expr_transcripts_fn, sep="\t", index=False)    

    else:
        gtex_most_expr_transcripts_df = pd.read_csv(gtex_most_expr_transcripts_fn, sep="\t")

    # extract transcript ids of the GTEx most expressed transcripts
    gtex_most_expr_transcripts_ls = [x.split(".")[0] for x in gtex_most_expr_transcripts_df['transcript_id'] if len(x.split(".")[0]) == 15 and "ENST" in x]
    # filter down to those which are in the Ensembl meta data dict
    gtex_most_expr_transcripts_reduced_ls = [x for x in gtex_most_expr_transcripts_ls if x in human_gtf_dict] # reduce to those which have an entry in our Ensembl GTF data

    # get chrom/meta data for faster searching of chip-seq peaks
    gtex_most_expr_transcripts_meta_data_ls = []
    for tid in gtex_most_expr_transcripts_reduced_ls:
        # chrom, strand, start, end, gene_name, gene_id, transcript_id, transcript_biotype = human_gtf_dict[tid]
        gtex_most_expr_transcripts_meta_data_ls.append(human_gtf_dict[tid])

    # create df
    gtex_most_expr_transcripts_meta_data_df = pd.DataFrame(gtex_most_expr_transcripts_meta_data_ls, columns=['chrom', 'strand', 'start', 'end', 'gene_name', 'gene_id', 'transcript_id', 'transcript_biotype'])
    gtex_most_expr_transcripts_meta_data_df['chrom'] = gtex_most_expr_transcripts_meta_data_df['chrom'].str.replace("chr","")

    return gtex_most_expr_transcripts_meta_data_df


def calc_promoter_occupancy():
    promoter_occupancy_both_dict = {}
    
    for target_end, target_end_d in {"5prime": {"+":'start', "-":'end'}, "3prime":{"+":'end', "-":'start'}}.items():

        promoter_occupancy_dict_fn = os.path.join(data_GTRD_dn, ".".join(["promoter_occupancy_" + target_end + "_dict", "_".join(["FC", str(gtex_fold_enrichment_cutoff)]), str(window_size), "json"]))

        if not os.path.exists(promoter_occupancy_dict_fn):
            print("Occupancy dict does not exist, starting run.")

            # track progress
            start_time = time.time()
            transcript_count = 0

            # initiate a dict to track occupancy of the promoter
            promoter_occupancy_dict = Counter({x:0 for x in range(int(window_size/2)*-1, int(window_size/2)+1)})
            promoter_occupancy_ls = []

            # iterate through chromosomes for limiting/matching space of comparisons
            for target_chr in gtrd_reduced_entries_df['#CHROM'].unique():
                print("current chromosome:", target_chr)

                # GTEX top transcripts for this chrom
                gtex_reduced_entries_chr_df = gtex_most_expr_transcripts_meta_data_df[gtex_most_expr_transcripts_meta_data_df['chrom'] == target_chr]
                print("len(gtex_reduced_entries_chr_df)", len(gtex_reduced_entries_chr_df))

                # GTRD ChIP-Seq peaks for this chrom
                gtrd_reduced_entries_chr_df = gtrd_reduced_entries_df[gtrd_reduced_entries_df['#CHROM'] == target_chr]
                print("len(gtrd_reduced_entries_chr_df)", len(gtrd_reduced_entries_chr_df))
                

                if len(gtrd_reduced_entries_chr_df) >0:
                    for i, high_expr_transcript_row in gtex_reduced_entries_chr_df.iterrows():
                        high_expr_transcript_strand = high_expr_transcript_row['strand']

                        # define the transcription start based on the strand of the target high-expr transcript
                        high_expr_transcript_start = int(high_expr_transcript_row[target_end_d[high_expr_transcript_strand]])

                        # define window start and end locs
                        window_start = int(high_expr_transcript_start - window_size/2)
                        window_end = int(high_expr_transcript_start + window_size/2)

                        # reduce GTRD entries to those who have a start or end within the window
                        gtrd_entries_inside_window = gtrd_reduced_entries_chr_df.query('START >= @window_start & START <= @window_end | END >= @window_start & END<= @window_end')

                        # iterate through GTRD entries within window and determine overlap
                        for gtrd_entry_inside_window in gtrd_entries_inside_window.to_records(index=False):
                            intersection_locs = overlap_range([gtrd_entry_inside_window[1], gtrd_entry_inside_window[2]+1], [window_start, window_end + 1])

                            # handle directional switching by strand
                            if high_expr_transcript_strand == "+":
                                tss_relative_loc_ls = [x - high_expr_transcript_start for x in intersection_locs]
                            elif high_expr_transcript_strand == "-":
                                tss_relative_loc_ls = [high_expr_transcript_start - x for x in intersection_locs] # reverse the TSS relative loc, as negative strand heads reverse

                            # add TSS relative locs to list
                            promoter_occupancy_ls += tss_relative_loc_ls

                        # update count
                        transcript_count += 1

                        # update count dict periodically to reduce working memory requirements
                        if transcript_count%100 == 0:
                            print(transcript_count, "gtex completed")
                            promoter_occupancy_ls_count = Counter(promoter_occupancy_ls)
                            promoter_occupancy_dict += promoter_occupancy_ls_count
                            promoter_occupancy_ls = []

            # final update of count dict
            promoter_occupancy_ls_count = Counter(promoter_occupancy_ls)
            promoter_occupancy_dict += promoter_occupancy_ls_count

            # save count dict
            save_json(promoter_occupancy_dict, promoter_occupancy_dict_fn)

        else:
            promoter_occupancy_dict = load_json(promoter_occupancy_dict_fn)
            promoter_occupancy_dict = {int(k):v for k,v in promoter_occupancy_dict.items()}

        # add target end dict to combined dict for subplotting
        promoter_occupancy_both_dict[target_end] = promoter_occupancy_dict

        # optionally find peaks in distribution
        #find_target_end_peaks(promoter_occupancy_dict)

    return promoter_occupancy_both_dict


def find_target_end_peaks(promoter_occupancy_dict):
    # find peaks
    y = [promoter_occupancy_dict[loc] for loc in range(-1* int(window_size/2), int(window_size/2)+1)]    
    peaks = find_peaks(y, height = 1000000, threshold = 1, distance = 100, prominence=100)


def occupancy_plot():
    print("Making subplot distplot.")
    fig, axs = plt.subplots(2, 1, figsize=(5,5), sharex=True)
    ax_c = 0
    for target_end, promoter_occupancy_dict in promoter_occupancy_both_dict.items():
        ax=axs[ax_c]
        keys = [int(x) for x in list(promoter_occupancy_dict.keys())]
        values = list(promoter_occupancy_dict.values())
        sns.histplot(x=keys, weights=values, bins=4000, color=distplot_color_d[target_end], ax=ax)
        ax.set_title(target_end)
        ax_c+=1

    plt.xticks(range(-1 * int(window_size/2), int(window_size/2) +1, 5*10**(floor(log10(window_size))-1)))
    plt.xticks(rotation=60)
    plt.savefig(os.path.join(output_figures_dn, ".".join(["combined", "promoter_occupancy", "distplot", "_".join(["FC", str(gtex_fold_enrichment_cutoff)]), str(window_size), "svg"])))
    plt.clf()

################################################################################
# Initiating Variables #########################################################
################################################################################

# raw data filenames
human_gtf_data_filename = [os.path.join(data_raw_ensembl_dn, x) for x in os.listdir(data_raw_ensembl_dn) if "chr.gtf" in x and ".gz" in x][0]
gtrd_peaks_dn = os.path.join(data_raw_GTRD_dn, "Homo_sapiens_ChIP-seq_peaks_MACS2")
gtrd_peaks_fns = [os.path.join(gtrd_peaks_dn, x) for x in os.listdir(gtrd_peaks_dn) if os.path.isfile(os.path.join(gtrd_peaks_dn, x)) and os.path.splitext(x)[1]==".interval" and "~" not in x]
gtex_transcripts_expr_fn = os.path.join(data_raw_GTEx_dn, "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz")

# thresholds
gtex_fold_enrichment_cutoff = 10 # limit the number of GTEx peaks parsed by level of fold over-enrichment
gtrd_fold_enrichment_cutoff = 10 # limit the number of GTEx peaks parsed by level of fold over-enrichment
window_size = 40000 # size of window centered on TSS for collecting GTEx peaks

# var
distplot_color_d = {"5prime":"turquoise", "3prime":"grey"}
################################################################################
# Execution ####################################################################
################################################################################

confirm = input("This analysis requires downloading of >11GB of data files. To continue confirm by entering 'YES' (without quote marks):")
if confirm == 'YES':

    # download the relevant data sets
    get_gtex_data()
    get_gtrd_data()

# get Ensembl transcript data from GTF
human_gtf_dict = parse_gtf()

# identify most expressed transcripts for each gene, create df with metadata
gtex_most_expr_transcripts_meta_data_df = gtex_get_expressed_transcripts()

# concat or load filtered GTRD peaks
gtrd_reduced_entries_df = gtrd_get_all_peaks()

# Promoter occupancy of ChIP-Seq
promoter_occupancy_both_dict = calc_promoter_occupancy()

# subplots of 5prime and 3prime ChIP-Seq occupancy
occupancy_plot()




    






                             
                         










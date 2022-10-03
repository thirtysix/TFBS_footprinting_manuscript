#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 3.8.0 ###########################################################
# Libraries ####################################################################
##from pathlib import Path
from assign_dirs import *
import os
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np

import csv
import json
from collections import Counter
from operator import itemgetter
import pprint
from scipy.stats import ttest_ind
from scipy.stats import ttest_rel
from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon

from statsmodels.stats.multitest import multipletests

import time
################################################################################
# Description/Notes ############################################################
################################################################################
"""

"""

################################################################################
# Base-level Functions #########################################################
################################################################################

def load_json(object_filename):
    """Load a json file to object."""

    if os.path.exists(object_filename):
        with open(object_filename, 'r') as object_file:
            return json.load(object_file)

def save_json(json_d, json_filename):
    """Save json object to file."""
    
    with open(json_filename, 'w') as json_file:
        json.dump(json_d, json_file)

def overlap_range(x,y):
    x.sort()
    y.sort()

    return range(max(x[0], y[0]), min(x[-1], y[-1])+1)

################################################################################
# Task-specific Functions ######################################################
################################################################################

def get_all_results(target_species, results_dns):
    """
    For target species: Read all results files and save, or load existing file.
    """

    target_all_results_fn = os.path.join(data_TFBS_analyses_dn, target_species, ".".join([target_species, "all_results", "parquet"]))

    if not os.path.exists(target_all_results_fn):
        print("run get_all_results", target_species, "file doesn't exist")
        print("...begin parsing")

        # get all result file names as list
        results_fns = [os.path.join(results_dn, "TFBSs_found.sortedclusters.csv") for results_dn in results_dns if os.path.exists(os.path.join(results_dn, "TFBSs_found.sortedclusters.csv"))]

        # concatenate results from all files into a single dataframe
        target_all_results_df = pd.concat((pd.read_csv(results_fn, usecols=['binding prot.', 'strand', 'start', 'end', 'TSS-relative start', 'TSS-relative end', 'PWM score', 'combined\naffinity\nscore', 'combined\naffinity\nscore\np-value']).assign(dirname = os.path.basename(os.path.dirname(results_fn))) for results_fn in results_fns))

        # add result key column for finding matched result in other_species
        print("Add result key column")
        target_all_results_df['result_key'] = target_all_results_df['dirname'] + "_" + target_all_results_df['strand'].astype(str) + "_" + target_all_results_df['start'].astype(str) + "_" + target_all_results_df['end'].astype(str)

        # create a column containing the range from start-end, to filter for those that overlap the variant position
        print("Add column for hits range")
        target_all_results_df['range'] = [list(range(i, j+1)) for i, j in target_all_results_df[['start', 'end']].values]

        # sort all results by combined affinity score (this order is used in the first step of determining critical pvals).
        target_all_results_df.sort_values(by='combined\naffinity\nscore\np-value', inplace=True)

        # save to parquet file
        target_all_results_df.to_parquet(target_all_results_fn)
        
    else:
        print("run get_all_results", target_species)
        print("...load existing file")
        target_all_results_df = pd.read_parquet(target_all_results_fn)

    return target_all_results_df


def get_worst_scores():
    """
    Calculate the worst combined affinity score after excluding all of the
    '>0.1' CAS pvalue entries; those with <=0.1 CAS pvalue to be substituted
    as the paired score when a putative positive hit is found in a target
    species but the other species lacks any significant hit.  This
    conservatively estimates differences in scores between species.
    """

    print("Calculating worst scores")
    human_worst_score_d = {tf_name:min(human_all_results_df[(human_all_results_df['binding prot.'].isin([tf_name])) & (human_all_results_df['combined\naffinity\nscore\np-value']!='>0.1')]['PWM score'].tolist()) for tf_name in human_all_results_df['binding prot.'].unique()}
    neand_worst_score_d = {tf_name:min(neand_all_results_df[(neand_all_results_df['binding prot.'].isin([tf_name])) & (neand_all_results_df['combined\naffinity\nscore\np-value']!='>0.1')]['PWM score'].tolist()) for tf_name in neand_all_results_df['binding prot.'].unique()}

    return human_worst_score_d, neand_worst_score_d


def determine_adjusted_pval(target_pvalues_ls, mtc_method, critical_fdr):
    """
    Determine a critical p-value for binding of each TF.
    """

    try:
        target_pvalues_ls.sort()
    except:
        print(len(target_pvalues_ls))
        print([x for x in target_pvalues_ls if isinstance(x, str)])
            
    # calculate bh critical pval for this list of pvals        
    bool_array, corrected_pvals, corrected_alpha_sidak, corrected_alpha_bonf = multipletests(target_pvalues_ls, method=mtc_method, alpha=critical_fdr, is_sorted=True)
    critical_pval = 0
    if any(bool_array):
        critical_pval_i = max(list(np.where(bool_array)[0]))
        critical_pval = target_pvalues_ls[critical_pval_i]
    
    return critical_pval


def parse_tfbs_results(target_species, target_all_results_df):
    """
    Get results satisfying pval threshold for the target organism (dict).
    """

    print("parse_tfbs_results", target_species)

    best_results_d = {tf_name:{} for tf_name in all_tfs}

    # edit the dataframe based on the desired slicing
    target_all_results_df_tmp = target_all_results_df.copy()

    all_accepted_pvals = []
    for result_description_i, result_description in enumerate(target_all_results_df_tmp['dirname'].unique()):        
        
        # generate a slice of all results df specific for this directory
        dirname_target_all_results_df_tmp = target_all_results_df_tmp.loc[np.in1d(target_all_results_df_tmp['dirname'], [result_description])]

        # keep track of progress
        if result_description_i>0:
            if result_description_i%100  == 0:
                print(result_description_i)             
        
        # ensure results exist for this directory
        if len(dirname_target_all_results_df_tmp) > 0:

            # iterate through each of the TFs which have results in the filtered df
            for tf_name in dirname_target_all_results_df_tmp['binding prot.'].unique():

                # new df is subset of all results for this dir which are for this TF
                tf_dirname_target_all_results_df_tmp = dirname_target_all_results_df_tmp.loc[np.in1d(dirname_target_all_results_df_tmp['binding prot.'], [tf_name])]

                # identify the best score for this TF in this results file
                tf_results_best_result = tf_dirname_target_all_results_df_tmp.iloc[0]
                tf_results_best_result_strand = tf_results_best_result['strand']
                tf_results_best_result_start = tf_results_best_result['start']
                tf_results_best_result_end = tf_results_best_result['end']
                tf_results_best_result_len = tf_results_best_result['end'] - tf_results_best_result['start']
                tf_results_best_result_pval = tf_results_best_result['combined\naffinity\nscore\np-value']
                result_key = tf_results_best_result['result_key']

                # get the calculated benjamini-hochberg critical p-value for this TF
                target_pvalues_ls = tf_dirname_target_all_results_df_tmp['combined\naffinity\nscore\np-value'].tolist()
                # account for entries which have been pre-filtered by PWM p-value during TFBSfootprinter run, on both strands
                # this increases the size of the pvalue list and results in a more stringent BH adjusted critical p-value
                target_pvalues_ls = target_pvalues_ls + [0.1]*(((window_size-tf_results_best_result_len) * 2) - len(target_pvalues_ls))
                # calculate critical pvalue for this TF in this analysis of this promoter
                critical_pval = determine_adjusted_pval(target_pvalues_ls, mtc_method, critical_fdr_tfbs)

                # keep results which satisfy critical pval
                if tf_results_best_result_pval <= critical_pval:
                    all_accepted_pvals.append(tf_results_best_result_pval)
                    # add score for this best hit to dict, indexed by TF and transcript result
                    best_results_d[tf_name][result_key] = tf_results_best_result['PWM score']

    all_accepted_pvals.sort()
    print("highest and lowest accepted pvals:", all_accepted_pvals[0], all_accepted_pvals[-1])

    return best_results_d


def compare_promoters(target_species, other_species):
    print("run compare_promoters", target_species, "vs.", other_species)

    # establish which of the species is being compared against the other
    if target_species == "human":
        target_all_results_df = human_all_results_df
        other_all_results_df = neand_all_results_df
        other_worst_score_d = neand_worst_score_d
        
    elif target_species == "neand":
        target_all_results_df = neand_all_results_df
        other_all_results_df = human_all_results_df
        other_worst_score_d = human_worst_score_d

    # get results satisfying pval threshold for the target organism (dict)
    target_best_results_d = parse_tfbs_results(target_species, target_all_results_df)
    
    # compare target vs. other ###########################################################################################
    paired_results_d = {}

    # iterate through tfs in target organism (dict) and corresponding best results in each directory (sub-dict)
    for tf_name, target_tf_best_results_d in target_best_results_d.items():

        # initiate entry in paired dict for this TF
        if tf_name not in paired_results_d:
            paired_results_d[tf_name] = {}

        # get the results df this TF in the other species
        other_tf_results_df = other_all_results_df.loc[np.in1d(other_all_results_df['binding prot.'], [tf_name])] # get all the results for this tf in the other_species results
        
        # iterate through each of the top hits (satisfying threshold) in the target species
        for target_tf_best_results_dn, target_tf_best_results_score in target_tf_best_results_d.items():

            # from the matching TF results, get those which correspond to the target's result key (key: directory, strand, TSS relative start, TSS relative end)
            other_tf_result_key_results_df = other_tf_results_df.loc[np.in1d(other_tf_results_df['result_key'], [target_tf_best_results_dn])]

            # if there are scores for this TF at this location (directory) in the other species, get the max score
            if len(other_tf_result_key_results_df) > 0:
                other_tf_best_results_score = max(other_tf_result_key_results_df['PWM score'].tolist())
            # otherwise get the worst score among the ALL the results for this TF which were not '>0.1'
            else:
                other_tf_best_results_score = other_worst_score_d[tf_name]

            # add the paired results [target species score, other species score] for this tf at this location
            paired_results_d[tf_name][target_tf_best_results_dn] = [target_tf_best_results_score, other_tf_best_results_score]

    # stat test ##########################################################################################################
    diff_ls_stattest = []
    for tf_name, target_other_scores_paired_d in paired_results_d.items():
        paired_d_values = target_other_scores_paired_d.values()
        print(tf_name, len(paired_d_values))
        if len(paired_d_values) >= n_threshold:
            
            target_scores_ls = [x[0] for x in paired_d_values]
            other_scores_ls = [x[1] for x in paired_d_values]
            
            # calculate wilcoxon stat
            stat_wil, pval_wil = wilcoxon(target_scores_ls, other_scores_ls, alternative=stat_alternative)

            # add entry
            diff_ls_stattest.append([tf_name, stat_wil, pval_wil, len(paired_d_values)])

    # output results for this analysis of thresholds
    stattest_df = pd.DataFrame(diff_ls_stattest, columns=["tf_name","stat","pval","n"])
    stattest_df.sort_values(by="pval", inplace=True)

    return paired_results_d, stattest_df


def do_multitests(target_species, critical_fdr, stattest_df):
    """Multiple testing correction"""
    
    wilcoxon_pvals = stattest_df['pval'].tolist()

    # target method compare target and other
    bool_array, corrected_pvals, corrected_alpha_sidak, corrected_alpha_bonf = multipletests(wilcoxon_pvals, method=mtc_method, alpha=critical_fdr, is_sorted=True)
    passing_count = len([x for x in bool_array if x])
    print(passing_count)
    passing_tfs = stattest_df.iloc[:passing_count]['tf_name'].tolist()

    # output results for this analysis of thresholds
    corrected_pval_colname = 'qval'
    stattest_df[corrected_pval_colname] = corrected_pvals
    diff_ls_sorted_stattest_fn = os.path.join(data_DBTFs_dn, ".".join([str(x) for x in [target_species, "wilcoxon_pwm.exact_loc", TSS_radius_th, critical_fdr_tfbs, critical_fdr_dbtb, n_threshold, stat_alternative, mtc_method, "tsv"]]))
    stattest_df.to_csv(diff_ls_sorted_stattest_fn, index=False, sep="\t")

    return stattest_df, passing_tfs

################################################################################
# Initiating Variables #########################################################
################################################################################
# dirnames
# TFBS_footprinter results dirnames
human_results_dns = [os.path.join(data_TFBS_analyses_human_results_dn, x) for x in os.listdir(data_TFBS_analyses_human_results_dn)]
neand_results_dns = [os.path.join(data_TFBS_analyses_neand_results_dn, x) for x in os.listdir(data_TFBS_analyses_neand_results_dn)]

# ensembl transcript id, gene id, gene name table
ens_transcript_fn = os.path.join(data_raw_ensembl_dn, "mart_export.transcript_id.gene_id.gene_name.txt")
ens_transcript_df = pd.read_csv(ens_transcript_fn, sep="\t")
ens_transcript_d = dict(zip(ens_transcript_df['Transcript stable ID'], ens_transcript_df['Gene name']))

# thresholds
##tfbs_result_pvalue_th = 0.01
critical_fdr_tfbs = 0.01
critical_fdr_dbtb = 0.10
mtc_method = 'fdr_bh'
stat_alternative = "greater"
TSS_radius_th = 2500 
window_size = 50 # size of window used in "neanderthal_variants" script
n_threshold = 5


################################################################################
# Execution ####################################################################
################################################################################
# Get results for both human and neanderthal
human_all_results_df = get_all_results("human", human_results_dns)
neand_all_results_df = get_all_results("neanderthal", neand_results_dns)
all_tfs = human_all_results_df['binding prot.'].unique()

# get the worst scores for each species - which are for a <0.1 CAS pvalue
human_worst_score_d, neand_worst_score_d = get_worst_scores()

# limit to those results which occur within TSS_radius_th bp (up/downstream) of TSS
print("Filter to hits inside target window")
human_all_results_df = human_all_results_df[abs((human_all_results_df['TSS-relative start'] + human_all_results_df['TSS-relative end'])/2)<=TSS_radius_th]
neand_all_results_df = neand_all_results_df[abs((neand_all_results_df['TSS-relative start'] + neand_all_results_df['TSS-relative end'])/2)<=TSS_radius_th]

# limit to those results which overlap the var position
print("Filter to hits which overlap variant")
human_var_mask = human_all_results_df['range'].apply(lambda x: int(window_size/2) in x)
human_all_results_df = human_all_results_df[human_var_mask]
neand_var_mask = neand_all_results_df['range'].apply(lambda x: int(window_size/2) in x)
neand_all_results_df = neand_all_results_df[neand_var_mask]

# sort result dfs by descending PWM score for easy comparison 
human_all_results_df.sort_values(by='PWM score', inplace=True, ascending=False)
neand_all_results_df.sort_values(by='PWM score', inplace=True, ascending=False)

# run the comparisons in both directions
passing_tfs_ls = []
comparison_pairs = [["human", "neand"],["neand", "human"]]
comparison_pairs = [["neand", "human"], ["human", "neand"]]
for comparison_pair in comparison_pairs:
    print("Starting comparison", comparison_pair)
    target_species, other_species = comparison_pair
    paired_results_d, stattest_df = compare_promoters(target_species, other_species)
    stattest_df, passing_tfs = do_multitests(target_species, critical_fdr_dbtb, stattest_df)
    passing_tfs_ls.append(passing_tfs)

    # for all tfs, output the paired scores and locations
    # for those results passing the p-value threshold in the target_species and for those tfs which pass wilcoxon testing for differential binding (stronger)
    tfs_results_out_ls = [] 
    for tf, tf_paired_results in paired_results_d.items():

        for tf_paired_result_k, tf_paired_result_v in tf_paired_results.items():
            tf_paired_result_k = tf_paired_result_k.replace("(", "_").replace(")_", "_")
            row = [tf] + tf_paired_result_k.split("_") + tf_paired_result_v
            tfs_results_out_ls.append(row)

    # output table of paired results
    tfs_results_out_df = pd.DataFrame(tfs_results_out_ls, columns = ["tf", "transcript", "region_start", "region_end", "frame_pval_cutoff", "strand", "hit_start", "hit_end", " ".join([target_species, "score"]), " ".join([other_species, "score"])])
    tfs_results_out_df['promoter_gene'] = [ens_transcript_d[x] if x in ens_transcript_d else "" for x in tfs_results_out_df['transcript'].tolist()]
    tfs_results_out_df['target_larger'] = tfs_results_out_df[" ".join([target_species, "score"])] > tfs_results_out_df[" ".join([other_species, "score"])]
    tfs_results_out_fn = os.path.join(data_DBTFs_dn, ".".join([str(x) for x in [target_species, other_species, "all_paired_results", TSS_radius_th, critical_fdr_tfbs, critical_fdr_dbtb, n_threshold, stat_alternative, mtc_method, "tsv"]]))
    tfs_results_out_df.to_csv(tfs_results_out_fn, sep="\t", index=False)

    # for all passing tfs, output the paired scores and locations
    # for those results passing the p-value threshold in the target_species and for those tfs which pass wilcoxon testing for differential binding (stronger)
    passing_tfs_results_out_ls = [] 
    for passing_tf in passing_tfs:
        passing_tf_paired_results = paired_results_d[passing_tf]

        for passing_tf_paired_result_k, passing_tf_paired_result_v in passing_tf_paired_results.items():
            passing_tf_paired_result_k = passing_tf_paired_result_k.replace("(", "_").replace(")_", "_")
            row = [passing_tf] + passing_tf_paired_result_k.split("_") + passing_tf_paired_result_v
            passing_tfs_results_out_ls.append(row)

    # output table of paired results
    passing_tfs_results_out_df = pd.DataFrame(passing_tfs_results_out_ls, columns = ["tf", "transcript", "region_start", "region_end", "frame_pval_cutoff", "strand", "hit_start", "hit_end", " ".join([target_species, "score"]), " ".join([other_species, "score"])])
    passing_tfs_results_out_df['promoter_gene'] = [ens_transcript_d[x] if x in ens_transcript_d else "" for x in passing_tfs_results_out_df['transcript'].tolist()]
    passing_tfs_results_out_df['target_larger'] = passing_tfs_results_out_df[" ".join([target_species, "score"])] > passing_tfs_results_out_df[" ".join([other_species, "score"])]
    passing_tfs_results_out_fn = os.path.join(data_DBTFs_dn, ".".join([str(x) for x in [target_species, other_species, "paired_results", TSS_radius_th, critical_fdr_tfbs, critical_fdr_dbtb, n_threshold, stat_alternative, mtc_method, "tsv"]]))
    passing_tfs_results_out_df.to_csv(passing_tfs_results_out_fn, sep="\t", index=False)



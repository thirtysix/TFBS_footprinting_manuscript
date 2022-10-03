#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 2.7.0 ###########################################################
# Libraries ####################################################################
from assign_dirs import *
import os
import csv
import pprint
import json
from operator import itemgetter

import math
import numpy as np
import pandas as pd
from numpy import median
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from matplotlib.cbook import boxplot_stats

import scipy.cluster.hierarchy as sch
################################################################################
# Description/Notes ############################################################
################################################################################
"""

"""

################################################################################
# Base-level Functions #########################################################
################################################################################


################################################################################
# Task-specific Functions ######################################################
################################################################################
def clean_jaspar_names(tf_names_ls):
    """
    Clean names of jaspar transcription factor names.
    MSX3 <- lost in humans.
    RHOX11 <- only present in 3 species.
    DUX <- mouse only gene.
    EWSR1 <- didn't end up in the Ensembl BioMart export.
    MIX-A <- jaspar says present in xenopus laevis.
    NR1A4 <- not in Ensembl for human.
    """

    special_dict = {"EWSR1-FLI1" : ["EWSR1","FLI1"]}
    names_list = []

    # split the combined names
    for uncleaned_jaspar_id in tf_names_ls:
        uncleaned_jaspar_id = uncleaned_jaspar_id.upper()
            
        split_names = uncleaned_jaspar_id.split("::")
        for name in split_names:
            names_list.append(name)

    # replace variants
    for i, name in enumerate(names_list):
        names_list[i] = name.replace("(VAR.2)","").replace("(VAR.3)","")

    tmp_list = []
    for i, name in enumerate(names_list):
        if name in special_dict:
            tmp_list += special_dict[name]
        else:
            tmp_list.append(name)

    names_list = list(set(tmp_list))
    names_list.sort()

    return names_list


def get_cleaned_target_tfs(target_tfs_filename):
    """
    Get DB TFs from results files.  Reduce/clean from JASPAR format to single gene names.
    """

    # capture expression across tissues in FANTOM dataset, for each of the target TFs
    target_tfs_cleaned = []
    
    # load DBTF stat test results
    target_tfs_df = pd.read_csv(target_tfs_filename, sep="\t")

    # reduce to those entries satisfying p-value/n thresholds
    target_tfs_df = target_tfs_df[(target_tfs_df['qval']<=critical_fdr) & (target_tfs_df['n']>=n_threshold)]

    # extract the JASPAR differentially binding TF names
    target_tfs = target_tfs_df['tf_name'].tolist()

    # clean the JASPAR-format target TF to its component TFs (e.g., "NR1H3::RXRA")
    file_target_tfs_cleaned = clean_jaspar_names(target_tfs)
    target_tfs_cleaned += file_target_tfs_cleaned

    # reduce to unique cleaned tf names
    target_tfs_cleaned = list(set(target_tfs_cleaned))

    return target_tfs_cleaned


def update_column_data(columns_metadata_df):

    # reformat the ages with leading zero
    raw_ages = columns_metadata_df['age'].tolist()
    formatted_ages = ["0"+x if len(x.split(" ")[0]) == 1 else x for x in raw_ages]
    columns_metadata_df['age'] = formatted_ages

    # new column for age category + structure
    age_categories = []
    for age in columns_metadata_df['age'].tolist():
        # categorize ages to 3 groups
        if age in pcw_levels:
            age_category = "pcw"
        elif age in mos_levels:
            age_category = "early"
        else:
            age_category = "late"

        age_categories.append(age_category)
        
    columns_metadata_df['age_category'] = age_categories
    columns_metadata_df['age_category+structure'] = columns_metadata_df['age_category']+"."+columns_metadata_df['structure_name']

    return columns_metadata_df


def heatmap():
    """
    Generate a heatmap.
    """

    tissue_names = columns_metadata_df['age_category+structure'].tolist()

    # slice tpm expr df to those genes which are in our DBTFs list
    gene_tissue_expr_df = tpm_expression_matrix_df.copy()

    # get row ids of DBTFs genes from rows metadata
    DBTF_rows_metadata_df = rows_metadata_df[rows_metadata_df['gene_symbol'].isin(human_target_tfs_cleaned)]
    DBTF_row_ids = DBTF_rows_metadata_df['row_num'].tolist()
    DBTF_row_gene_names = DBTF_rows_metadata_df['gene_symbol']

    # slice expr df to those genes which are in our DBTFs list
    DBTF_gene_tissue_expr_df = gene_tissue_expr_df[gene_tissue_expr_df[0].isin(DBTF_row_ids)]
    print(DBTF_gene_tissue_expr_df)
    # remove row id col (at pos 0) from expr df
    DBTF_gene_tissue_expr_df.drop(DBTF_gene_tissue_expr_df.columns[0], axis=1, inplace=True)

    # add tissue names as columns
    DBTF_gene_tissue_expr_df.columns = tissue_names

    # tissue names are age category + structure, which means duplicates.  Average them.
    # appropriated from: https://stackoverflow.com/questions/40311987/pandas-mean-of-columns-with-the-same-names/40312254
    import numbers
    DBTF_gene_tissue_expr_df = DBTF_gene_tissue_expr_df.groupby(by=DBTF_gene_tissue_expr_df.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])

    # add a pseudocount prior to log transformation
    DBTF_gene_tissue_expr_df = DBTF_gene_tissue_expr_df+0.1

    # log transform df
    DBTF_gene_tissue_expr_log_df = np.log2(DBTF_gene_tissue_expr_df)

    # add column of gene symbols
    DBTF_gene_tissue_expr_log_df.set_index(DBTF_row_gene_names, inplace=True)

    # life stages, corresponding to age category
    life_stages = [x.split(".")[0] for x in DBTF_gene_tissue_expr_log_df.columns]
    
    # transform df
    DBTF_gene_tissue_expr_log_df_T = DBTF_gene_tissue_expr_log_df.T
    
    # set row colors
    lut_d = {"pcw":sns.color_palette("PiYG", 10)[0], "early":sns.color_palette("PiYG", 10)[2], "late":sns.color_palette("PiYG", 10)[9]}
    lut_d = {"pcw":sns.color_palette("Reds_r", 10)[0], "early":sns.color_palette("Reds_r", 10)[5], "late":sns.color_palette("Reds_r", 10)[8]}
    row_colors = [lut_d[life_stage] for life_stage in life_stages]

    # make clustermap
    sns.set(font_scale=0.25)
    row_dist = sch.distance.pdist(DBTF_gene_tissue_expr_log_df_T)
    row_linkage = sch.linkage(row_dist, method='ward', metric="euclidean")
    col_dist = sch.distance.pdist(DBTF_gene_tissue_expr_log_df_T.T)
    col_linkage = sch.linkage(col_dist, method='ward', metric="euclidean")
    g = sns.clustermap(DBTF_gene_tissue_expr_log_df_T, row_linkage=row_linkage, col_linkage=col_linkage, standard_scale=0, cmap="mako", row_colors=row_colors)
    col_clusters = sch.fcluster(col_linkage, 0.5*col_dist.max(), 'distance')
##    plt.tick_params(axis='both', which='major', labelsize=14)

    # save clustermap
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    fig_filename = os.path.join(output_figures_ABA_brainspan_dn, ".".join(["ABA_brainspan", "cluster", "004","svg"]))
    figsize=(50,25)
    plt.savefig(fig_filename, bbox_inches = "tight")

    # get and save labels for other plots
    ylabels = [y.get_text() for y in list(g.ax_heatmap.yaxis.get_majorticklabels())]
    xlabels = [x.get_text() for x in list(g.ax_heatmap.xaxis.get_majorticklabels())]

    # export xlabels
    xlabels_ofn = os.path.join(output_figures_ABA_brainspan_dn, "xlabels_ABA.csv")
    xlabels_df = pd.DataFrame(xlabels, columns = ["tf_names"])
    xlabels_df.to_csv(xlabels_ofn, sep="\t", index=False)

    # export clusters
    genes_clusters_ofn = os.path.join(output_figures_ABA_brainspan_dn, "clusters_ABA.tsv")
    genes_clusters_df = pd.DataFrame(zip(list(DBTF_gene_tissue_expr_log_df_T.columns), col_clusters), columns = ["tf_names", "cluster"])
    genes_clusters_df.to_csv(genes_clusters_ofn, sep="\t", index=False)

    # clear the plot
    plt.clf()

    return DBTF_gene_tissue_expr_log_df_T, ylabels, xlabels, genes_clusters_df


def plot_cluster_expression():
    """
    Creat a plot of expression for each gene cluster, over timepoints.
    """

    # unique clusters
    clusters = genes_clusters_df['cluster'].unique()
    clusters.sort()
    # iterate through clusters
    for cluster in clusters:
        print(cluster)

        # initiate ls for custom data related to this cluster
        cluster_df_ls = []

        # DBTF genes in this cluster
        cluster_genes = genes_clusters_df[genes_clusters_df['cluster']==cluster]['tf_names'].tolist()

        # get data for each sample, for each gene (old code; but works)
        for i, row in columns_metadata_df.iterrows():
            age_weeks_log = row['age_weeks_log']
            age_weeks = row['age_weeks']
            structure_name = row['structure_name']

            # get data for each gene in this cluster in this sample
            for cluster_gene in cluster_genes:
                cluster_gene_expr_log = DBTF_gene_tissue_expr_log_T_df.iloc[i][cluster_gene]
                if isinstance(cluster_gene_expr_log, pd.Series):
                    cluster_gene_expr_log = sum(cluster_gene_expr_log)/len(cluster_gene_expr_log)
                cluster_df_ls.append([age_weeks, age_weeks_log, structure_name, cluster_gene, cluster_gene_expr_log])

        # create df of targeted data for this cluster
        cluster_df = pd.DataFrame(cluster_df_ls, columns=['age_weeks', 'age_weeks_log', 'structure_name', 'cluster_gene', 'expr_log'])
        # limit to those structure which were targeted to have several timepoints, and which occur at <=92 weeks.
        cluster_df = cluster_df[cluster_df['structure_name'].isin(target_structure_names)]
        cluster_df = cluster_df[cluster_df['age_weeks']<=144]

        # build plot
        fig, axes = plt.subplots(2, 1, sharex=True, figsize=(20,12), gridspec_kw={'height_ratios': [4, 1]})

        # lineplot
        ax1=axes[0]
        age_descriptions_d = {8:"early fetal", 9:"early fetal", 12:"early fetal", 13:"early mid-fetal", 16:"early mid-fetal", 17:"early mid-fetal", 19:"early mid-fetal", 21:"early mid-fetal", 22:"early mid-fetal", 24:"late fetal", 25:"late fetal", 26:"late fetal", 35:"late fetal", 37:"late fetal", 56:"neonatal and early infancy", 80:"late infancy", 92:"early childhood", 144:"early childhood"}
        sns.lineplot(x='age_weeks', y='expr_log', data=cluster_df, ax=ax1)
        xticks = sorted(list(set(cluster_df['age_weeks'].tolist())))
        xticks_labels = ["-".join([str(x)+"w", age_descriptions_d[x]]) for x in xticks]
        ax1.set_xticks(xticks)
        ax1.set_xticklabels(xticks_labels, rotation=30, ha='right', size=10)
        ax1.set_xlabel("Age (weeks)", fontsize=30)
        ax1.set_ylabel("TPM (log)", fontsize=30)
        ax1.set_title(" ".join(["Cluster", str(cluster), "genes expression during development"]), fontsize=20)
        ax1.tick_params(axis='both', which='major', labelsize=14)

        # add bars for age categories
        ax2=axes[1]
        bars = [(8, 4), (12, 12), (24, 13), (37,43), (80,64)]
        bar_colors = sns.color_palette("Reds_r", len(bars))
        ax2.broken_barh(bars, (0, 2), facecolors=bar_colors)
        ax2.grid(False)
        ax2.set_yticks([])

        # add labels for the barplots
        age_descriptions_texts = ["Early fetal (8-12w)", "Mid fetal (12-24w)", "Late fetal (24-37w)", "Infancy", "Childhood"]
        for age_descriptions_text, age_range, bar_color in zip(age_descriptions_texts, bars, bar_colors):
            text_color = 'w'
            bar_color_lumi = math.sqrt(0.299*(bar_color[0]**2) + 0.587*(bar_color[1]**2) + 0.114*(bar_color[2]**2))
            if bar_color_lumi > 0.65:
                text_color='dimgray'

            # add the text
            ax2.text(age_range[0]+(age_range[1]/2), 0.2, age_descriptions_text, rotation=90, color=text_color, fontsize=14)

        # aesthetics
        plt.tight_layout()

        # save figure
        fig_filename = os.path.join(output_figures_ABA_brainspan_dn, ".".join(["ABA_brainspan", "cluster", str(cluster), "lineplot", "svg"]))
        plt.savefig(fig_filename)

################################################################################
# Initiating Variables #########################################################
################################################################################
# variables from DBTF analysis
tfbs_result_pvalue_th = 0.01
critical_fdr = 0.10
mtc_method = 'fdr_bh'
stat_alternative = "greater"
TSS_radius_th = 2500 
window_size = 50 # size of window used in "neanderthal_variants" script
n_threshold = 5

# ABA brainspan data
columns_metadata_filename = os.path.join(data_raw_ABA_brainspan_dn, "columns_metadata.csv")
columns_metadata_df = pd.read_csv(columns_metadata_filename, sep=",")
rows_metadata_filename = os.path.join(data_raw_ABA_brainspan_dn, "rows_metadata.csv")
rows_metadata_df = pd.read_csv(rows_metadata_filename, sep=",")
expression_matrix_filename = os.path.join(data_raw_ABA_brainspan_dn, "expression_matrix.jaspar.tpm.parquet")
tpm_expression_matrix_df = pd.read_parquet(expression_matrix_filename)

# ABA brainspan categories
pcw_levels = ['08 pcw', '09 pcw', '12 pcw', '13 pcw', '16 pcw', '17 pcw', '19 pcw', '21 pcw', '24 pcw', '25 pcw', '26 pcw', '35 pcw', '37 pcw', '04 mos']
mos_levels = ['10 mos','01 yrs','02 yrs','03 yrs','04 yrs']
years_levels = ['08 yrs', '11 yrs', '13 yrs', '15 yrs', '18 yrs', '19 yrs', '21 yrs', '23 yrs', '30 yrs', '36 yrs', '37 yrs', '40 yrs']

# DBTFs
human_target_tfs_filename = os.path.join(data_DBTFs_dn, ".".join([str(x) for x in ["human", "wilcoxon_pwm.exact_loc", str(TSS_radius_th), tfbs_result_pvalue_th, critical_fdr, n_threshold, stat_alternative, mtc_method, "tsv"]]))
human_target_tfs_cleaned = get_cleaned_target_tfs(human_target_tfs_filename) # base names of DB TFs

################################################################################
# Execution ####################################################################
################################################################################
# update the column metadata
columns_metadata_df = update_column_data(columns_metadata_df)
# create heatmap
DBTF_gene_tissue_expr_log_df_T, ylabels, xlabels, genes_clusters_df = heatmap()

# clusters expression
columns_metadata_df.drop('column_num', inplace=True, axis=1)
week_convert_dict = {"pcw":[1, 0], "mos":[4, 40], "yrs":[52, 40]}
life_stage_order_weeks = [int(x.split(" ")[0]) *  week_convert_dict[x.split(" ")[1]][0] + week_convert_dict[x.split(" ")[1]][1] for x in columns_metadata_df['age'].tolist()]
life_stage_order_weeks_log = [round(math.log(x, 2), 3) for x in life_stage_order_weeks]
columns_metadata_df['age_weeks'] = life_stage_order_weeks
columns_metadata_df['age_weeks_log'] = life_stage_order_weeks_log

# tpm expression data
# get row ids of DBTFs genes from rows metadata
DBTF_rows_metadata_df = rows_metadata_df[rows_metadata_df['gene_symbol'].isin(human_target_tfs_cleaned)]
DBTF_row_ids = DBTF_rows_metadata_df['row_num'].tolist()
DBTF_row_gene_names = DBTF_rows_metadata_df['gene_symbol']
DBTF_gene_tissue_expr_df = tpm_expression_matrix_df[tpm_expression_matrix_df[0].isin(DBTF_row_ids)]
DBTF_gene_tissue_expr_df.drop(DBTF_gene_tissue_expr_df.columns[0], axis=1, inplace=True)
DBTF_gene_tissue_expr_df.columns = columns_metadata_df['age_category+structure'].tolist()
DBTF_gene_tissue_expr_df = DBTF_gene_tissue_expr_df + 0.1
DBTF_gene_tissue_expr_log_df = np.log2(DBTF_gene_tissue_expr_df)
DBTF_gene_tissue_expr_log_T_df = DBTF_gene_tissue_expr_log_df.T
DBTF_gene_tissue_expr_log_T_df.columns = list(DBTF_row_gene_names)

# identify structures with several timepoints
structure_names = columns_metadata_df[columns_metadata_df['age_weeks']<=92]['structure_name'].unique()
target_structure_names = [structure_name for structure_name in structure_names if len(columns_metadata_df[(columns_metadata_df['age_weeks']<=144) & (columns_metadata_df['structure_name']==structure_name)]['age_weeks'].unique()) >= 12 ]

# plot expression over time for each gene cluster
plot_cluster_expression()




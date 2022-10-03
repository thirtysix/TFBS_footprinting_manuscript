#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 3.7.5 ###########################################################
# Libraries ####################################################################
import os
import time
from datetime import datetime
from shutil import copyfile
from assign_dirs import *

from timeit import default_timer as timer

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
import math

import itertools

################################################################################
# Description/Notes ############################################################
################################################################################
"""

"""

################################################################################
# Base-level Functions #########################################################
################################################################################
def copy_script():
    # make a copy of the current script file and save it to the output dir
    cur_script_filename = __file__
    now = datetime.now()
    datetime_str = now.strftime("%Y-%m-%d--%H.%M.%S")
    cur_script_datetime_filename = cur_script_filename+"."+datetime_str
    output_cur_script_datetime_filename = os.path.join(figure_dirname, os.path.basename(cur_script_datetime_filename))
    print("copied scriptfilename", output_cur_script_datetime_filename)
    copyfile(cur_script_filename, output_cur_script_datetime_filename)

def export_to_parquet(adata):
    import pyarrow.parquet as pq
    import pyarrow as pa

    #pd.DataFrame(data=adata.raw.X, index=adata.obs_names, columns=adata.raw.var_names).to_csv('adata_raw_x.csv')
    df = pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names)
    df_pa_table = pa.Table.from_pandas(df)
    parquet_filename = "./adata_x.parquet"
    pq.write_table(df_pa_table, parquet_filename, use_dictionary=True, compression='snappy')


##def export_raw_w():
##    
##    # export the germline data for coexpression analysis
##    sc.pp.normalize_total(adata, target_sum=1e6)
##    adata_germ = adata[adata.obs['TissueName']=='Germline']
##    adata_germ.raw = adata_germ
##    t = adata_germ.raw.X.toarray()
##    outdf = pd.DataFrame(data=t, index=adata_germ.obs_names, columns=adata_germ.raw.var_names)
##    outdf = outdf.transpose()
##    outdf.to_csv('adata_raw_x.germline.csv')


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
################################################################################
# Task-specific Functions ######################################################
################################################################################



################################################################################
# Figure Variables #############################################################
################################################################################

def set_dirs():
    # name dirs & set outdir
    figure_dirname = os.path.join(output_figures_ABA_cortical_dn, ".".join(["figures", analysis_description]))
    qc_figure_dirname = os.path.join(figure_dirname, "QC")
    hvg_dirname = os.path.join(figure_dirname, "HVG")
    tsne_figure_dirname = os.path.join(figure_dirname, "tsne")
    umap_figure_dirname = os.path.join(figure_dirname, "umap")
    fa_figure_dirname = os.path.join(figure_dirname, "fa")
    dotplot_figure_dirname = os.path.join(figure_dirname, "dotplot")
    paga_figure_dirname = os.path.join(figure_dirname, "paga")
    barplot_figure_dirname = os.path.join(figure_dirname, "celltype_barplot")
    violinplot_figure_dirname = os.path.join(figure_dirname, "violinplot")
    matrixplot_figure_dirname = os.path.join(figure_dirname, "matrixplot")
    heatmapplot_figure_dirname = os.path.join(figure_dirname, "heatmapplot")
    tracksplot_figure_dirname = os.path.join(figure_dirname, "tracksplot")
    rank_genes_tsv_dirname = os.path.join(figure_dirname, "rank_genes_tsv")
    other_figure_dirname = os.path.join(figure_dirname, "other")

    # make the dirs
    for mk_dirname in [figure_dirname, qc_figure_dirname, hvg_dirname, umap_figure_dirname, fa_figure_dirname, dotplot_figure_dirname, paga_figure_dirname, barplot_figure_dirname, violinplot_figure_dirname, matrixplot_figure_dirname, heatmapplot_figure_dirname, tracksplot_figure_dirname, other_figure_dirname, rank_genes_tsv_dirname]:
        if not os.path.exists(mk_dirname):
            os.mkdir(mk_dirname)

    return figure_dirname, qc_figure_dirname, hvg_dirname, tsne_figure_dirname, umap_figure_dirname, fa_figure_dirname, dotplot_figure_dirname, paga_figure_dirname, barplot_figure_dirname, violinplot_figure_dirname, matrixplot_figure_dirname, heatmapplot_figure_dirname, tracksplot_figure_dirname, other_figure_dirname, rank_genes_tsv_dirname


def figure_params():
    #Define a nice colour map for gene expression
    colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
    colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
    colorsComb = np.vstack([colors3, colors2])
    mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

    # figures settings
    sc.settings.figdir = figure_dirname # set output dir
    sc.set_figure_params(dpi_save=1000)
    sns_color_pal = "RdYlBu"

    return mymap, sns_color_pal


def pretty_colors(obs_name):
    
    unique_obs_count = len(list(set(adata.obs[obs_name].to_list())))
    obs_colors = sns.color_palette(sns_color_pal, unique_obs_count)

    return obs_colors


def target_obs_marker_genes_heatmap(obs):
    """
    Make heatmap of expression in current obs categories, using the target genes
    which are marker genes in categories of this obs.
    Return list of clustered ordered target marker genes.
    """
    # heatmap of target_obs_marker_genes

    # get the marker genes associated with this obs
    target_obs_marker_genes = obs_marker_genes_d[obs]

    # initiate dict for holding target marker genes
    combined_d = {}

    # iterate through target marker genes and get expr vals
    for target_obs_marker_gene in target_obs_marker_genes:
        expr_ls = [x[0] for x in adata[:,target_obs_marker_gene].X.todense().tolist()]
        combined_d[target_obs_marker_gene] = expr_ls

    # save the obs values
    combined_d[obs] = adata.obs[obs].tolist()

    # create df
    combined_df = pd.DataFrame(combined_d)

    # use means for each obs level across all genes
    fig_ls = []
    for obs_n in combined_df[obs].unique():
        row_vals = []
        row_vals.append(obs_n)
        
        obs_combined_df = combined_df[combined_df[obs]==obs_n]
  
        for target_obs_marker_gene in target_obs_marker_genes:
            target_obs_marker_gene_mean = obs_combined_df[target_obs_marker_gene].mean()
            row_vals.append(target_obs_marker_gene_mean)

        fig_ls.append(row_vals)

    # make figure
    fig_df = pd.DataFrame(fig_ls, columns = ["obs_n"] + target_obs_marker_genes)
    fig_df.pop("obs_n")
    fig = plt.figure(figsize=(5, 20))
    cm = sns.clustermap(fig_df, method='ward', metric='euclidean')
    plt.tight_layout()
    plt.savefig(os.path.join(other_figure_dirname, ".".join([obs, "marker_genes_heatmap.svg"])))

    # get coexpression ordered gene list back
    xlabels = [x.get_text() for x in list(cm.ax_heatmap.xaxis.get_majorticklabels())]

    return xlabels


def cortical_layer_heatmap(obs):
    # cortical layer heatmap
    # get the ordered cluster labels - naming depends on scanpy version
    ordered_cluster_label_ls = adata.uns['dendrogram_'+obs]['categories_ordered']
##    ordered_cluster_label_ls = adata.uns['dendrogram_['+obs+"]"]['categories_ordered']
    # get the cortical layers
    cortical_layers_ls = sorted(list(adata.obs['cortical_layer_label'].cat.categories))

    cluster_label_layer_ls_ls = []
    for cluster_label in ordered_cluster_label_ls:

        cluster_label_layer_ls = []
        cluster_label_df = adata.obs[adata.obs[obs]==cluster_label]

        for cortical_layer in cortical_layers_ls:
            cluster_label_layer_ls.append(len(cluster_label_df[cluster_label_df['cortical_layer_label'] == cortical_layer])+1)


        cluster_label_layer_ls_sum = sum(cluster_label_layer_ls)
        cluster_label_layer_ls_max = max(cluster_label_layer_ls)
        
        cluster_label_layer_ls_ratio = [x/cluster_label_layer_ls_max for x in cluster_label_layer_ls]
        cluster_label_layer_ls_ls.append(cluster_label_layer_ls_ratio)


    cluster_label_layer_df = pd.DataFrame(cluster_label_layer_ls_ls, columns = cortical_layers_ls, index=ordered_cluster_label_ls)

    fig = plt.figure(figsize=(5, 20))
    ax = sns.heatmap(cluster_label_layer_df, cmap="YlGnBu", square=True, linewidths=.5)
    plt.tight_layout()
    plt.savefig(os.path.join(other_figure_dirname, obs+".cortical_layer.max.svg"))


def cortical_region_heatmap(obs):
    # cortical region heatmap
    # get the ordered cluster labels - naming depends on scanpy version    
    ordered_cluster_label_ls = adata.uns['dendrogram_'+obs]['categories_ordered']
##    ordered_cluster_label_ls = adata.uns['dendrogram_['+obs+"]"]['categories_ordered']
    # get the cortical layers
    regions_ls = sorted(list(adata.obs['region_label'].cat.categories))

    cluster_label_layer_ls_ls = []
    for cluster_label in ordered_cluster_label_ls:

        cluster_label_layer_ls = []
        cluster_label_df = adata.obs[adata.obs[obs]==cluster_label]

        for region in regions_ls:
            cluster_label_layer_ls.append(len(cluster_label_df[cluster_label_df['region_label'] == region])+1)

        cluster_label_layer_ls_sum = sum(cluster_label_layer_ls)
        cluster_label_layer_ls_max = max(cluster_label_layer_ls)
        
        cluster_label_layer_ls_ratio = [x/cluster_label_layer_ls_max for x in cluster_label_layer_ls]
        cluster_label_layer_ls_ls.append(cluster_label_layer_ls_ratio)


    cluster_label_layer_df = pd.DataFrame(cluster_label_layer_ls_ls, columns = regions_ls, index=ordered_cluster_label_ls)

    fig = plt.figure(figsize=(5, 20))
    ax = sns.heatmap(cluster_label_layer_df, cmap=sns.color_palette("ch:s=-.2,r=.6", as_cmap=True), square=True, linewidths=.5)
    plt.tight_layout()
    plt.savefig(os.path.join(other_figure_dirname, obs+".region.max.svg"))


def nuclei_barplot(obs):
    # nuclei barplot
    ordered_cluster_label_ls = adata.uns['dendrogram_'+obs]['categories_ordered']
##    ordered_cluster_label_ls = adata.uns['dendrogram_['+obs+"]"]['categories_ordered']
    nuclei_ls = []

    for cluster_label in ordered_cluster_label_ls:

        cluster_label_df = adata.obs[adata.obs[obs]==cluster_label]
        nuclei_ls.append([cluster_label, len(cluster_label_df)])

    nuclei_df = pd.DataFrame(nuclei_ls, columns=[obs, 'N Nuclei'])
        
    fig = plt.figure(figsize=(5, 20))
    bp = sns.barplot(data=nuclei_df, y=obs, x='N Nuclei', color='grey', saturation=0.5, orientation='horizontal', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(other_figure_dirname, obs+".N_nuclei.svg"))
    plt.clf()
    
    
################################################################################
# Initiating Variables #########################################################
################################################################################
# multiple data files
dataset_dir = data_raw_ABA_cortical_dn

# single data file
dataset_fn = os.path.join(dataset_dir, "matrix.csv")
dataset_h5ad_fn = os.path.join(dataset_dir, "matrix.h5ad")
features_fn = os.path.join(dataset_dir, "metadata.csv")
features_df = pd.read_csv(features_fn)
tsne_fn = os.path.join(dataset_dir, "tsne.csv")
tsne_df = pd.read_csv(tsne_fn)

# name your targets
# DBTFs
# variables from DBTF analysis
tfbs_result_pvalue_th = 0.01
critical_fdr = 0.10
mtc_method = 'fdr_bh'
stat_alternative = "greater"
TSS_radius_th = 2500 
window_size = 50 # size of window used in "neanderthal_variants" script
n_threshold = 5
human_target_tfs_filename = os.path.join(data_DBTFs_dn, ".".join([str(x) for x in ["human", "wilcoxon_pwm.exact_loc", str(TSS_radius_th), tfbs_result_pvalue_th, critical_fdr, n_threshold, stat_alternative, mtc_method, "tsv"]]))
target_genes = get_cleaned_target_tfs(human_target_tfs_filename) # base names of DB TFs

# dataset characteristics
species = 'human'
# dataset observations of interest
obss = ["cluster_label", "class_label", "subclass_label", "region_label", "cortical_layer_label"] 

# User-defined variables
primary_outformat = "svg"
dotsize=1.25
n_top_genes = 5000 # top genes to analyze if relevant
n_comps=15 # PCA components
n_neighbors = 15 # neighbors
hvg_flavor = "seurat_v3" # highly variable genes algo
batch_correction_key = "cortical_layer_label"
experiment_description = "ABA_cortical_areas_precommit"

# descriptions used to make and set dirnames
accessory_description_ls = ["DEG_"+str(n_top_genes),
                            "Ncomps_"+str(n_comps),
                            "NN_" + str(n_neighbors),
                            hvg_flavor]
accessory_description = ".".join(accessory_description_ls)
analysis_description = "_".join([experiment_description, accessory_description])


figure_dirname,qc_figure_dirname,hvg_dirname,tsne_figure_dirname,umap_figure_dirname,fa_figure_dirname,dotplot_figure_dirname,paga_figure_dirname,barplot_figure_dirname,violinplot_figure_dirname,matrixplot_figure_dirname,heatmapplot_figure_dirname,tracksplot_figure_dirname,other_figure_dirname,rank_genes_tsv_dirname = set_dirs()
# save a copy of the script to figure output dir
copy_script()
# set figure params
mymap, sns_color_pal = figure_params()

################################################################################
# Execution ####################################################################
################################################################################
# LOAD DATA ####################################################################
# TSV/CSV
##if not os.path.exists(dataset_h5ad_fn):
##    adata = sc.read_csv(dataset_fn, delimiter=",") # dense TSV/CSV file
##    adata.X = csr_matrix(adata.X)
##    adata.write_h5ad(dataset_h5ad_fn)
##
##else:
##    adata = sc.read_h5ad(dataset_h5ad_fn) # load dataset

adata = sc.read_csv(dataset_fn, delimiter=",") # dense TSV/CSV file
adata.X = csr_matrix(adata.X)

print(adata)

# make var_names (cell ids) unique
adata.var_names_make_unique()


# CHARACTERISTICS ##############################################################
# get sample characteristics if needed
# add target sample feature metadata
for obs in obss:
    adata.obs[obs] = features_df[obs].tolist()

# original experiment clusters
adata.obs['cluster_order'] = features_df['cluster_order'].tolist()
adata.obs['cluster_order'] = adata.obs['cluster_order'].astype('category')

# get gene names and reduce target genes to those which actually occur in this dataset
gene_names = list(adata.var_names)
target_genes = [x for x in target_genes if x in gene_names]
print(target_genes)

# Quality control - calculate QC covariates
# remove original experiment identified outliers
adata.obs["outlier_call"] = features_df["outlier_call"].tolist()
adata = adata[adata.obs['outlier_call']!=True]

# add cluster label top level
adata.obs['cluster_label'] = adata.obs['cluster_label'].fillna("")
adata.obs['cluster_label_top'] = [" ".join(x.split(" ")[:3]) for x in adata.obs['cluster_label'].tolist()]
adata.obs['cluster_label_top'] = adata.obs['cluster_label_top'].astype('category')
obss.append('cluster_label_top')

# add cluster label top agg
adata.obs['cluster_label_top_agg'] = [" ".join(x.split(" ")[:2]) for x in adata.obs['cluster_label'].tolist()]
adata.obs['cluster_label_top_agg'] = adata.obs['cluster_label_top_agg'].astype('category')
obss.append('cluster_label_top_agg')

# add region layer class
adata.obs['region_layer_class'] = [" ".join(x) for x in zip(adata.obs['region_label'].tolist(), adata.obs['cortical_layer_label'].tolist(), adata.obs['class_label'].tolist())]
adata.obs['region_layer_class'] = adata.obs['region_layer_class'].astype('category')
obss.append('region_layer_class')


# import original experiment tsne, rather than calculate
adata.obsm['X_tsne'] = tsne_df[['tsne_1', 'tsne_2']].values

# counts
adata.obs['n_counts'] = adata.X.sum(1)
adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
adata.obs['n_genes'] = (adata.X > 0).sum(1)

## Quality control - plot QC metrics
##Sample quality plots
sc.settings.figdir = qc_figure_dirname # set output dir
groupby=""
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
##sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'mt_frac'], jitter=0.4, multi_panel=True, save=".".join(["counts", "panel", primary_outformat]), show=False)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True, save=".".join(["counts", "panel", primary_outformat]), show=False)
t1 = sc.pl.violin(adata, 'n_counts', size=2, log=True, cut=0, save=".".join(["counts", groupby, primary_outformat]), show=False)
##t2 = sc.pl.violin(adata, 'mt_frac', save=".".join(["mtfrac", groupby, primary_outformat]), show=False)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True, save=".".join(["counts", "panel", primary_outformat]), show=False)
t1 = sc.pl.violin(adata, 'n_counts', size=2, log=True, cut=0, save=".".join(["counts", groupby, primary_outformat]), show=False)


###Data quality summary plots
####p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac', save=os.path.join(qc_figure_dirname, ".".join(["n_counts_n_genes", primary_outformat])), show=False)
####p2 = sc.pl.scatter(adata[adata.obs['n_counts']<15000], 'n_counts', 'n_genes', save=".".join(["n_counts_n_genes", "lower", primary_outformat]), color='mt_frac', show=False)
##p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', save=".".join(["n_counts_n_genes", primary_outformat]), show=False)
##plt.show()
##
###Thresholding decision: counts
##p3 = sns.distplot(adata.obs['n_counts'], kde=False)
##plt.show()
##p4 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<500000], kde=False, bins=40)
##plt.show()
##p5 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']>1000000], kde=False, bins=40)
##plt.show()
##
###Thresholding decision: genes
##p6 = sns.distplot(adata.obs['n_genes'], kde=False, bins=40)
##plt.show()
##p7 = sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']<1000], kde=False, bins=40)
##plt.show()


# FILTERING ####################################################################
min_counts = 75000
min_counts = 300000
max_counts = 3500000
mt_frac_upper_threshold = 0.2
min_genes = 2000
max_genes = 100000
min_cells = adata.n_obs * 0.005
min_cells = 50

# starting counts
print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))

# Filter cells according to identified QC thresholds:
sc.pp.filter_cells(adata, min_counts = min_counts)
print('Number of cells after min count filter: {:d}'.format(adata.n_obs))
sc.pp.filter_cells(adata, max_counts = max_counts)
print('Number of cells after max count filter: {:d}'.format(adata.n_obs))
##adata = adata[adata.obs['mt_frac'] < mt_frac_upper_threshold]
##print('Number of cells after MT filter: {:d}'.format(adata.n_obs))

# Filter cells by gene counts
sc.pp.filter_cells(adata, min_genes = min_genes)
print('Number of cells after gene filter: {:d}'.format(adata.n_obs))
##sc.pp.filter_cells(adata, max_genes = max_genes)
##print('Number of cells after gene filter: {:d}'.format(adata.n_obs))

# Filter genes by cells expressed
sc.pp.filter_genes(adata, min_cells=min_cells)
print('Number of cells after all filter: {:d}'.format(adata.n_obs))
print('Number of genes after all filter: {:d}'.format(adata.n_vars))

### save adata to parquet
##export_to_parquet(adata)

# Refresh target genes based on newly filtered gene list #######################
gene_names = list(adata.var_names)
target_genes = [x for x in target_genes if x in gene_names]


# NORMALIZE & LOG-TRANSFORM ####################################################
adata.layers['counts'] = adata.X.copy()
print("NORMALIZE")
start=timer()
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
adata.raw = adata
end=timer()
print(f'{(end-start):.4}', "seconds")


# VARIABLE GENES ###############################################################
print("HIGHLY VARIABLE GENES")
sc.settings.figdir = hvg_dirname
start=timer()
sc.pp.highly_variable_genes(adata, flavor=hvg_flavor, n_top_genes=n_top_genes, batch_key=batch_correction_key, layer='counts') # for use with seurat_v3
end=timer()
print(f'{(end-start):.4}', "seconds")
print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
sc.pl.highly_variable_genes(adata, save=".".join(['highly_variable', primary_outformat]), show=False)


# PCA & NEAREST NEIGHBORS ######################################################
print("PCA")
start=timer()
sc.pp.pca(adata, n_comps=n_comps, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=n_neighbors)
end=timer()
print(f'{(end-start):.4}', "seconds")


# CALCULATE VISUALIZATIONS  ####################################################
# TSNE positions have been calculated in original experiment, and are imported
print("UMAP")
start=timer()
sc.tl.umap(adata, min_dist=0.1)
end=timer()
print(f'{(end-start):.4}', "seconds")

print("FORCE DIRECTED MODEL")
start=timer()
sc.tl.draw_graph(adata)
end=timer()
print(f'{(end-start):.4}', "seconds")


# RANK GENES ###################################################################
obs_marker_genes_d = {}

for obs in obss:
    print(obs)

    # RANK GENES
    rank_genes_key = 'rank_genes.'+str(obs)
    sc.tl.rank_genes_groups(adata, groupby=obs, key_added=rank_genes_key)
    rank_genes_result = adata.uns[rank_genes_key]

    # can plot ranked genes if interested
    #sc.pl.rank_genes_groups(adata, key=rank_genes_key, n_genes=25, sharey=False, save=".".join(["rank_gene_groups", str(obs), primary_outformat]), show=False, use_raw=False)

    # save ranked genes to table
    rank_genes_ofn = os.path.join(rank_genes_tsv_dirname, rank_genes_key+".tsv")
    rank_genes_df = pd.DataFrame(rank_genes_result['names']).head(200)
    rank_genes_df.to_csv(rank_genes_ofn, sep="\t", index=False) # export the x top ranked genes per cluster

    # save ranked genes full to table
    rank_genes_full_ofn = os.path.join(rank_genes_tsv_dirname, rank_genes_key+".full.tsv")
    groups = rank_genes_result['names'].dtype.names
    print(groups)
    rank_genes_full_df = pd.DataFrame({group + '_' + key[:1]: rank_genes_result[key][group] for group in groups for key in ['names','logfoldchanges','scores','pvals_adj']})
    rank_genes_full_df.to_csv(rank_genes_full_ofn, sep="\t", index=False) # all top ranked genes and enrichment/pvalue per cluster

    # check target genes in differentially expressed genes 
    target_marker_genes = []
    for cluster_id in rank_genes_df.columns:
        cluster_top_degs = rank_genes_df[cluster_id].tolist()[:100]
        cluster_top_degs_target_genes = [x for x in target_genes if x in cluster_top_degs]
        print(cluster_id, len(cluster_top_degs_target_genes), cluster_top_degs_target_genes)
        target_marker_genes+=cluster_top_degs_target_genes

    # save to dict
    obs_marker_genes_d[obs] = list(set(target_marker_genes))
    print(obs, list(set(target_marker_genes)))


### PLOT TARGET GENES | FORCE-DIRECTED & UMAP #####################################
## There are many target genes (100+) only output FA and UMAP figures for each if desired, uncomment
##sc.settings.figdir = fa_figure_dirname # set output dir
##for target_gene in target_genes:
##    sc.pl.draw_graph(adata, color=target_gene, legend_loc='on data', title=target_gene, frameon=True, save="."+ ".".join([target_gene, primary_outformat]), show=False, use_raw=False, color_map=mymap, size=dotsize)
##
##sc.settings.figdir = umap_figure_dirname
##for target_gene in target_genes:
##    sc.pl.umap(adata, color=target_gene, legend_loc='on data', title=target_gene, frameon=True, save="."+ ".".join([target_gene, primary_outformat]), show=False, use_raw=False, color_map=mymap, size=dotsize)


# OBS PLOTS ####################################################################
sc.set_figure_params(dpi_save=600) # for these larger figures, lower the DPI

for obs in obss:

    # overlap in target genes and marker genes for this obs
    target_obs_marker_genes = obs_marker_genes_d[obs]

    # creation of the heatmap produces an cluster-ordered list of target genes
    target_obs_marker_genes = target_obs_marker_genes_heatmap(obs)

    # colors based on number of categories in this obs
    obs_colors = pretty_colors(obs)

    # calculate dendrogram for this obs
    sc.tl.dendrogram(adata, groupby=obs, cor_method='spearman', linkage_method='ward', var_names=['mean_counts'], key_added='dendrogram_'+obs)

    # FA
    sc.settings.figdir = fa_figure_dirname # set output dir
    sc.pl.draw_graph(adata, color=obs, legend_loc='on data', title=obs, frameon=True, save="."+ ".".join([obs, primary_outformat]), show=False, use_raw=False, palette=obs_colors, legend_fontsize = 4, size=dotsize, add_outline=True, outline_width=(0.15, 0.025), edges=True, edges_width=0.01)
    sc.pl.draw_graph(adata, color=obs, title=obs, frameon=True, save="."+ ".".join([obs, "legend_off", primary_outformat]), show=False, use_raw=False, palette=obs_colors, legend_fontsize = 4, size=dotsize, add_outline=True, outline_width=(0.15, 0.025), edges=True, edges_width=0.01)

    # dotplot    
    sc.settings.figdir = dotplot_figure_dirname # set output dir
    sc.pl.dotplot(adata, target_obs_marker_genes, groupby=obs, save=".".join([obs, 'dendrogram', primary_outformat]), show=False, use_raw=False, dendrogram=True)

    # matrixplot
    sc.settings.figdir = matrixplot_figure_dirname
    sc.pl.matrixplot(adata, target_obs_marker_genes, groupby=obs, save=".".join([obs, 'dendrogram_scaled', 'marker_genes', primary_outformat]), show=False, use_raw=False, dendrogram=True, standard_scale='var', cmap='Reds')
    sc.pl.matrixplot(adata, target_obs_marker_genes, groupby=obs, save=".".join([obs, 'dendrogram', 'marker_genes', primary_outformat]), show=False, use_raw=False, dendrogram=True, cmap='Reds')

    # custom heatmaps using Seaborn; will be output to other_figure_dirname
    cortical_layer_heatmap(obs)
    cortical_region_heatmap(obs)
    nuclei_barplot(obs)


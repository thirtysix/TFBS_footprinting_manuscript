#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 3.8.0 ###########################################################
# Libraries ####################################################################
from assign_dirs import *
import os
import pprint
import json
import csv
from operator import itemgetter
##from pathlib import Path

import pandas as pd
##import seaborn as sns
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from matplotlib.cbook import boxplot_stats

import math
import numpy as np
from numpy import median

from collections import Counter
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


def get_protein_atlas_enriched_tissues_genes():
    """
    ProteinAtlas lists which genes have enriched expression in various tissues.
    """

    # protein atlas enriched tissues
    all_elevated_tissues_raw = protein_atlas_data['RNA tissue specific nTPM'].tolist()
    all_elevated_tissues_pairs_lsls = [x.split(";") for x in all_elevated_tissues_raw]
    all_elevated_tissues_pairs_ls = [tissue_pair for all_elevated_tissues_ls in all_elevated_tissues_pairs_lsls for tissue_pair in all_elevated_tissues_ls]
    all_elevated_tissues_ls = list(set([x.split(": ")[0] for x in all_elevated_tissues_pairs_ls]))
    print(all_elevated_tissues_ls)
    
    # assign ProteinAtlas bulk RNA-Seq tissues to categories
    tissue_categories_dict = {
        "blood & immune":['blood', 'bone marrow', 'lymphoid tissue', ],
        "brain, eye & retina":["brain", 'retina'],
        "endocrine & exocrine":['adrenal gland', 'parathyroid gland', 'pancreas', 'pituitary gland', 'thyroid gland',],
        "female/male":['breast', 'cervix, uterine', 'endometrium 1', 'fallopian tube', 'ovary', 'placenta', 'vagina', 'ductus deferens', 'epididymis', 'prostate', 'seminal vesicle', 'testis', 'cervix'],
        "gastrointestinal":['esophagus', 'gallbladder', 'intestine', 'pancreas', 'salivary gland', 'stomach 1', 'tongue', 'liver'],
        "urinary":['kidney', 'urinary bladder']}

    # identify any tissues which are in the proteinatlas 'RNA tissue specific nTPM' column which are not accounted for in our dict
    tissues_lsls =  list(tissue_categories_dict.values())
    tissues_ls = [x for subls in tissues_lsls for x in subls]
    print("Not categorized proteinatlas tissues", [x for x in all_elevated_tissues_ls if x not in tissues_ls])
    
    for tissue_group, tissue_terms in tissue_categories_dict.items():
        tissue_group_protein_atlas_data_df = protein_atlas_data[protein_atlas_data['RNA tissue specific nTPM'].str.contains("|".join(tissue_terms))]

        # get the primary gene names
        protein_atlas_tissue_genes_dict[tissue_group] = tissue_group_protein_atlas_data_df['Gene'].tolist()

        # add the gene synonyms - which exist as a list
        gene_synonyms_raw = tissue_group_protein_atlas_data_df['Gene synonym'].tolist()
        gene_synonyms_raw_lsls = []
        for x in gene_synonyms_raw:
            gene_synonyms_raw_lsls.append(x.split(", "))
            
        gene_synonyms = [synonym for gene_synonyms_raw_ls in gene_synonyms_raw_lsls for synonym in gene_synonyms_raw_ls]
        protein_atlas_tissue_genes_dict[tissue_group] += gene_synonyms

    return tissue_categories_dict


def get_protein_atlas_enriched_celltypes_genes():
    """
    ProteinAtlas lists which genes have enriched expression in various cell types.
    """

    # file associating cell types to top level tissue
    celltype_categories_fn = os.path.join(data_raw_ProteinAtlas_dn, "proteinatlas_celltypes_tissue.tsv")
    celltype_categories_df = pd.read_csv(celltype_categories_fn, sep="\t")
    celltype_categories_dict = {tissue_category:celltype_categories_df[celltype_categories_df['tissue group']==tissue_category]['name'].to_list() for tissue_category in tissue_categories_dict.keys()}

    # for each tissue group and associated cell type terms, add gene names and gene name synonyms to dic
    for tissue_group, celltype_terms in celltype_categories_dict.items():
        if len(celltype_terms) > 0:
            tissue_group_protein_atlas_data_df = protein_atlas_data[protein_atlas_data['RNA single cell type specific nTPM'].str.contains("|".join(celltype_terms))]

            # get the primary gene names
            protein_atlas_tissue_genes_dict[tissue_group]+=tissue_group_protein_atlas_data_df['Gene'].tolist()

            # add the gene synonyms - which exist as a list
            gene_synonyms_raw = tissue_group_protein_atlas_data_df['Gene synonym'].tolist()
            gene_synonyms_raw_lsls = []
            for x in gene_synonyms_raw:
                gene_synonyms_raw_lsls.append(x.split(", "))
                
            gene_synonyms = [synonym for gene_synonyms_raw_ls in gene_synonyms_raw_lsls for synonym in gene_synonyms_raw_ls]
            protein_atlas_tissue_genes_dict[tissue_group] += gene_synonyms


def get_protein_atlas_enriched_brain_celltypes_genes():
    """
    ProteinAtlas lists which genes have enriched expression in brain cell types.
    """
    
    # protein atlas enriched brain region
    protein_atlas_data_brain = protein_atlas_data[protein_atlas_data['RNA brain regional specific nTPM']!=""]
    # get the primary gene names
    protein_atlas_tissue_genes_dict['brain, eye & retina']+=protein_atlas_data_brain['Gene'].tolist()

    # add the gene synonyms - which exist as a list
    gene_synonyms_raw = protein_atlas_data_brain['Gene synonym'].tolist()
    gene_synonyms_raw_lsls = []
    for x in gene_synonyms_raw:
        gene_synonyms_raw_lsls.append(x.split(", "))
        
    gene_synonyms = [synonym for gene_synonyms_raw_ls in gene_synonyms_raw_lsls for synonym in gene_synonyms_raw_ls]
    protein_atlas_tissue_genes_dict['brain, eye & retina'] += gene_synonyms


def get_protein_atlas_enriched_blood_celltypes_genes():
    """
    ProteinAtlas lists which genes have enriched expression in blood cell types.
    """

    # protein atlas enriched immune & blood
    protein_atlas_data_blood = protein_atlas_data[protein_atlas_data['RNA blood cell specific nTPM']!=""]
    # get the primary gene names
    ##protein_atlas_tissue_genes_dict['immune']+=protein_atlas_data_brain['Gene'].tolist()
    protein_atlas_tissue_genes_dict['blood & immune']=protein_atlas_data_blood['Gene'].tolist()

    # add the gene synonyms - which exist as a list
    gene_synonyms_raw = protein_atlas_data_blood['Gene synonym'].tolist()
    gene_synonyms_raw_lsls = []
    for x in gene_synonyms_raw:
        gene_synonyms_raw_lsls.append(x.split(", "))
        
    gene_synonyms = [synonym for gene_synonyms_raw_ls in gene_synonyms_raw_lsls for synonym in gene_synonyms_raw_ls]
    protein_atlas_tissue_genes_dict['blood & immune'] += gene_synonyms


def get_protein_atlas_enriched_cellline_genes():
    """
    ProteinAtlas lists which genes have enriched expression in cell lines.
    """
    # protein atlas enriched cell line
    cellline_categories_fn = os.path.join(data_raw_ProteinAtlas_dn, "proteinatlas_celllines_tissue.tsv")
    cellline_categories_df = pd.read_csv(cellline_categories_fn, sep="\t")
    cellline_categories_dict = {tissue_category:cellline_categories_df[cellline_categories_df['tissue group']==tissue_category]['name'].to_list() for tissue_category in tissue_categories_dict.keys()}

    for tissue_group, cellline_terms in cellline_categories_dict.items():
        if len(cellline_terms) > 0:
            tissue_group_protein_atlas_data_df = protein_atlas_data[protein_atlas_data['RNA cell line specific nTPM'].str.contains("|".join(cellline_terms))]

            # get the primary gene names
            protein_atlas_tissue_genes_dict[tissue_group]+=tissue_group_protein_atlas_data_df['Gene'].tolist()

            # add the gene synonyms - which exist as a list
            gene_synonyms_raw = tissue_group_protein_atlas_data_df['Gene synonym'].tolist()
            gene_synonyms_raw_lsls = []
            for x in gene_synonyms_raw:
                gene_synonyms_raw_lsls.append(x.split(", "))
                
            gene_synonyms = [synonym for gene_synonyms_raw_ls in gene_synonyms_raw_lsls for synonym in gene_synonyms_raw_ls]
            protein_atlas_tissue_genes_dict[tissue_group] += gene_synonyms

    
def get_DICE_immune_genes():

    """
    Parse DICE immune project expression data, and identify genes with above average expression in any of the immune cell types.
    """

    # ensembl gene ids - for matching with DICE entries
    # can be downloaded directly: http://www.ensembl.org/biomart/martview/4d9ab7b71bd18892033d0ebbff1d8d70?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.external_gene_name&FILTERS=&VISIBLEPANEL=filterpanel
    ensembl_fn = os.path.join(data_raw_ensembl_dn, "mart_export.gene_id.gene_name.tsv")
    ensembl_df = pd.read_csv(ensembl_fn, sep="\t")
    ensembl_df.columns=["ens_id", "gene_name"]

    # parse the immune gene names
    # can be downloaded directly: https://dice-database.org/download/mean_tpm_merged.csv
    immunity_genes_filename = os.path.join(data_raw_DICE_dn, "mean_tpm_merged.csv")
    immunity_genes_df = pd.read_csv(immunity_genes_filename, sep=",")
    immunity_genes_df_columns = list(immunity_genes_df.columns)
    immunity_genes_df_columns[0] = "ens_id"
    immunity_genes_df.columns = immunity_genes_df_columns

    # merge the immunity genes and the ensembl data - so that there are gene names associated
    immunity_genes_df = pd.merge(immunity_genes_df, ensembl_df, on="ens_id")

    # get 90th quan of all columns - immune cell type
    immunity_genes_quan_df = immunity_genes_df.quantile(0.95)

    # identify genes with expression above cell type quantile, across any cell type
    above_quan_ls = []
    cell_type_ls = []
    for i, row in immunity_genes_df.iterrows():
        above_average = False
        cell_types = []

        # check each column (cell type)
        for colname, value in row.items():
            if colname in immunity_genes_quan_df:
                col_thresh = immunity_genes_quan_df[colname]
                if value >=col_thresh:
                    above_average = True
                    cell_types.append(colname)

        above_quan_ls.append(above_average)
        cell_type_ls.append("|".join(cell_types))
        
    # above quan qualifier, and associated cell types, to df
    immunity_genes_df["above_quan"] = above_quan_ls
    immunity_genes_df["cell_type"] = cell_type_ls

    DICE_immunity_genes = immunity_genes_df[immunity_genes_df["above_quan"]==True]["gene_name"].tolist()

    return DICE_immunity_genes


def get_box_genes():
    """
    Genes in homeobox, FOX, and SOX gene groups.  As defined by: https://www.genenames.org/data/genegroup/#!/
    """
    
    homeobox_genes_ls = ['CERS1','CERS2','CERS3','CERS4','CERS5','CERS6','CUX1','CUX2','CUX2P1','CUX2P2','ONECUT1','ONECUT2','ONECUT3','SATB1','SATB2','HMBOX1','HNF1A','HNF1B','CDX1','CDX2','CDX4','EVX1','EVX2','GBX1','GBX2','GSX1','GSX2','HOXA1','HOXA10','HOXA11','HOXA13','HOXA2','HOXA3','HOXA4','HOXA5','HOXA6','HOXA7','HOXA9','HOXB1','HOXB13','HOXB2','HOXB3','HOXB4','HOXB5','HOXB6','HOXB7','HOXB8','HOXB9','HOXC10','HOXC11','HOXC12','HOXC13','HOXC4','HOXC5','HOXC6','HOXC8','HOXC9','HOXD1','HOXD10','HOXD11','HOXD12','HOXD13','HOXD3','HOXD4','HOXD8','HOXD9','MEOX1','MEOX2','MNX1','PDX1','ISL1','ISL2','LHX1','LHX2','LHX3','LHX4','LHX5','LHX6','LHX8','LHX9','LMX1A','LMX1B','BARHL1','BARHL2','BARX1','BARX2','BSX','DBX1','DBX2','DLX1','DLX2','DLX3','DLX4','DLX5','DLX6','EMX1','EMX2','EN1','EN2','HHEX','HLX','HMX1','HMX2','HMX3','LBX1','LBX2','MSX1','MSX2','MSX2P1','NANOG','NANOGP1','NANOGP10','NANOGP11','NANOGP2','NANOGP3','NANOGP4','NANOGP5','NANOGP6','NANOGP7','NANOGP8','NANOGP9','NKX1-1','NKX1-2','NKX2-1','NKX2-2','NKX2-3','NKX2-4','NKX2-5','NKX2-6','NKX2-8','NKX3-1','NKX3-2','NKX6-1','NKX6-2','NKX6-3','NOTO','TLX1','TLX2','TLX3','VAX1','VAX2','VENTX','VENTXP1','VENTXP2','VENTXP3','VENTXP4','VENTXP5','VENTXP6','VENTXP7','HDX','POU1F1','POU2F1','POU2F2','POU2F3','POU3F1','POU3F2','POU3F3','POU3F4','POU4F1','POU4F2','POU4F3','POU5F1','POU5F1B','POU5F1P2','POU5F1P3','POU5F1P4','POU5F1P5','POU5F1P6','POU5F1P7','POU5F2','POU6F1','POU6F2','ALX1','ALX3','ALX4','ARGFX','ARGFXP1','ARGFXP2','ARX','CRX','DMBX1','DPRX','DPRXP1','DPRXP2','DPRXP3','DPRXP4','DPRXP5','DPRXP6','DPRXP7','DRGX','DUX4L1','DUX4L10','DUX4L11','DUX4L12','DUX4L13','DUX4L14','DUX4L15','DUX4L16','DUX4L17','DUX4L18','DUX4L19','DUX4L2','DUX4L20','DUX4L21','DUX4L22','DUX4L23','DUX4L24','DUX4L25','DUX4L26','DUX4L27','DUX4L28','DUX4L29','DUX4L3','DUX4L31','DUX4L32','DUX4L33','DUX4L34','DUX4L35','DUX4L36','DUX4L37','DUX4L4','DUX4L40','DUX4L41','DUX4L42','DUX4L43','DUX4L44','DUX4L45','DUX4L46','DUX4L47','DUX4L48','DUX4L49','DUX4L5','DUX4L50','DUX4L51','DUX4L52','DUX4L6','DUX4L7','DUX4L8','DUX4L9','DUXA','DUXAP1','DUXAP10','DUXAP2','DUXAP3','DUXAP4','DUXAP5','DUXAP6','DUXAP7','DUXAP8','DUXAP9','DUXB','ESX1','GSC','GSC2','HESX1','HOPX','ISX','LEUTX','MIXL1','NOBOX','OTP','OTX1','OTX2','OTX2P1','PAX2','PAX3','PAX4','PAX5','PAX6','PAX7','PAX8','PHOX2A','PHOX2B','PITX1','PITX2','PITX3','PROP1','PRRX1','PRRX2','RAX','RAX2','RHOXF1','RHOXF2','RHOXF2B','SEBOX','SHOX','SHOX2','TPRX1','TPRX1P1','TPRX2P','TPRXL','UNCX','VSX1','VSX2','PROX1','PROX2','SIX1','SIX2','SIX3','SIX4','SIX5','SIX6','IRX1','IRX1P1','IRX2','IRX3','IRX4','IRX4P1','IRX5','IRX6','MEIS1','MEIS2','MEIS3','MEIS3P1','MEIS3P2','MKX','PBX1','PBX2','PBX2P1','PBX3','PBX4','PKNOX1','PKNOX2','TGIF1','TGIF1P1','TGIF2','TGIF2LX','TGIF2LY','TGIF2P1','ADNP','ADNP2','HOMEZ','TSHZ1','TSHZ2','TSHZ3','ZEB1','ZEB2','ZEB2P1','ZFHX2','ZFHX3','ZFHX4','ZHX1','ZHX2','ZHX3']
    sox_genes_ls = ['SOX1','SOX2','SOX3','SOX4','SOX5','SOX6','SOX7','SOX8','SOX9','SOX10','SOX11','SOX12','SOX13','SOX14','SOX15','SOX17','SOX18','SOX21','SOX30','SRY']
    fox_genes_ls = ['FOXA1','FOXA2','FOXA3','FOXB1','FOXB2','FOXC1','FOXC2','FOXD1','FOXD2','FOXD3','FOXD4','FOXD4L1','FOXD4L3','FOXD4L4','FOXD4L5','FOXD4L6','FOXE1','FOXE3','FOXF1','FOXF2','FOXG1','FOXH1','FOXI1','FOXI2','FOXI3','FOXJ1','FOXJ2','FOXJ3','FOXK1','FOXK2','FOXL1','FOXL2','FOXM1','FOXN1','FOXN2','FOXN3','FOXN4','FOXO1','FOXO3','FOXO4','FOXO6','FOXP1','FOXP2','FOXP3','FOXP4','FOXQ1','FOXR1','FOXR2','FOXS1']
    box_genes_ls = homeobox_genes_ls+sox_genes_ls+fox_genes_ls

    return box_genes_ls


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
    

def heatmap(target_tfs_fantom_cage_exp_data_df):
    """
    Generate a heatmap.
    """

    # log-normalize TPM expression vals
    expr_df = target_tfs_fantom_cage_exp_data_df.copy()
    expr_df.set_index('gene_name', inplace=True)
    expr_df = expr_df+0.1
    expr_df = expr_df.apply(np.log2)
    tissue_names = expr_df.columns

    # identify life stage of each tissue
    life_stages = []
    for tissue_name in tissue_names:
        if "fetal" in tissue_name.lower() or "newborn" in tissue_name.lower():
            life_stages.append("fetal/newborn")
        else:
            life_stages.append("adult")
        
    # Re-orient and add tissue names and lifestages data
    gene_tissue_expr_dict = expr_df.T.to_dict(orient='list')
    gene_tissue_expr_dict["tissue"] = tissue_names
    gene_tissue_expr_dict["life_stage"] = life_stages

    # build new dataframe for heatmap
    heatmap_df = pd.DataFrame(data=gene_tissue_expr_dict)
    print(heatmap_df)
    heatmap_df.index = heatmap_df['tissue']
    heatmap_df.pop("tissue")
    life_stages = heatmap_df.pop("life_stage")

    # set row colors
    lut = dict(zip(life_stages.unique(), [sns.color_palette("RdYlBu", 10)[3], sns.color_palette("PiYG", 10)[2]]))
    row_colors = life_stages.map(lut)
    
    # make clustermap
    sns.set(font_scale=0.25)

    
    g = sns.clustermap(heatmap_df, method=method, metric=metric, standard_scale=0, cmap="mako", row_colors=row_colors)

    # save clustermap
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    fig_filename = os.path.join(output_figures_FANTOM_clustermap_dn, ".".join([species, metric, method, "cluster", "svg"]))
    figsize=(50,25)
    plt.savefig(fig_filename, bbox_inches = "tight")

    # get and save clustered labels for other plots
    ylabels = [y.get_text() for y in list(g.ax_heatmap.yaxis.get_majorticklabels())]
    xlabels = [x.get_text() for x in list(g.ax_heatmap.xaxis.get_majorticklabels())]
    
    # output the clustered gene order to file
    xlabels_ofn = os.path.join(output_figures_FANTOM_clustermap_dn, ".".join([species, "xlabels.csv"]))
    xlabels_df = pd.DataFrame(xlabels, columns=["tf_names"])
    xlabels_df.to_csv(xlabels_ofn, index=False)

    # clear the plot
    plt.clf()

    return ylabels, xlabels


def tf_genes_class_colors(xlabels):
    """
    Generate a multi-plot showing presence of target genes in those tissue groups
    for there is evidence of enrichment.
    """

    # set colors
    colors_palette_dict = tissue_category_color_dict

    # tissue groups of interest
    tissue_groups = ["developmental", "brain, eye & retina", "blood & immune", "female/male", "gastrointestinal", "urinary", "endocrine & exocrine"]

    # initiate figure and subplots
    fig, axes = plt.subplots(nrows=1, ncols=len(tissue_groups), sharey=True, figsize=(1,10))
    
    # iterate through tissues
    for i, tissue_group in enumerate(tissue_groups):

        # set subplot axis
        ax = axes[i]
        ax.set_title(tissue_group, rotation=90)

        # get genes enriched in this tissue group
        genes_list = protein_atlas_tissue_genes_dict[tissue_group]

        # palette of colors depending on presence of gene in list of tissue group genes
        gene_class_ls = [colors_palette_dict[tissue_group] if xlabel in genes_list else sns.color_palette("RdBu_r", 7)[3] for xlabel in xlabels]

        # make a df for easier plotting
        gene_class_df = pd.DataFrame(data={"gene":xlabels, "expression":[1 for x in xlabels]})

        # make ax barplot
        bp = sns.barplot(x="expression", y="gene", data=gene_class_df, palette=gene_class_ls, edgecolor="none", ax=ax, linewidth=0, ci=None)
        bp.set(xticklabels=[])

        # labels
        ax.set_ylabel("")
        ax.set_xlabel("")
        
    # adjust plot aesthetics
    plt.subplots_adjust(wspace=0)
    plt.xticks([])

    # save plot
    fig_filename = os.path.join(output_figures_FANTOM_clustermap_dn, ".".join([species, metric, method, "gene_tissues_groups", "combined", "svg"]))
    fig.savefig(fig_filename)
    plt.clf()


# unique DB TFs by species
def get_unique_db_TFs():
    unique_db_TFs_d = {}
    human_target_tfs_filename = os.path.join(data_DBTFs_dn, ".".join([str(x) for x in ["human", "wilcoxon_pwm.exact_loc", str(TSS_radius_th), tfbs_result_pvalue_th, critical_fdr, str(n_threshold), stat_alternative, mtc_method, "tsv"]]))
    neand_target_tfs_filename = os.path.join(data_DBTFs_dn, ".".join([str(x) for x in ["neand", "wilcoxon_pwm.exact_loc", str(TSS_radius_th), tfbs_result_pvalue_th, critical_fdr, str(n_threshold), stat_alternative, mtc_method, "tsv"]]))
    human_target_tfs_cleaned = get_cleaned_target_tfs(human_target_tfs_filename) # base names of DB TFs
    neand_target_tfs_cleaned = get_cleaned_target_tfs(neand_target_tfs_filename) # base names of DB TFs
    unique_human_db_tfs = [x for x in human_target_tfs_cleaned if x not in neand_target_tfs_cleaned]
    unique_neand_db_tfs = [x for x in neand_target_tfs_cleaned if x not in human_target_tfs_cleaned]
    unique_db_TFs_d = {"human":unique_human_db_tfs, "neand":unique_neand_db_tfs}

    return unique_db_TFs_d


def agg_expression_barplot():
    """
    Aggregate expression of target TFs in top x FANTOM tissues.
    """

    # get FANTOM expression data and reduce to DB TFs
    expr_df = target_tfs_fantom_cage_exp_data_df.copy()
    expr_df.pop("gene_name")

    # sum expression for each tissue, across the DB TFs
    barplot_sum_df = pd.DataFrame(expr_df.sum(axis=0), columns=["expression"])
    print(barplot_sum_df)

    # assign colors to each of the tissue groups
    clrs = []
    for ylabel in ylabels:
        tissue_groups_color = "white"
        
        if ylabel in FANTOM_tissue_categories_d:
            tissue_group = FANTOM_tissue_categories_d[ylabel]

            if tissue_group in tissue_category_color_dict:
                tissue_groups_color = tissue_category_color_dict[tissue_group]
            else:
                print("Undefined tissue group color:", ylabel)
        else:
            print("Unknown tissue group:", ylabel)
            
            
        clrs.append(tissue_groups_color)
        
    # build figure
    fig = plt.figure()
    bp = sns.barplot(x=barplot_sum_df.index, y="expression", data = barplot_sum_df, order = ylabels, palette=clrs, edgecolor = "none")
    fig = bp.get_figure()
    plt.xticks(rotation=90)
    plt.tight_layout()

    # save figure
    fig_filename = os.path.join(output_figures_FANTOM_clustermap_dn, ".".join([species, "agg_barplot", "svg"]))
    fig.savefig(fig_filename)
    plt.clf()


def agg_expression_barplot_all_tissues():
    """
    Aggregate expression of target TFs in top x FANTOM tissues.
    """

    # get FANTOM expression data and reduce to DB TFs
    expr_df = fantom_cage_exp_data_df[fantom_cage_exp_data_df['gene_name'].isin(target_tfs_cleaned)]
    expr_df.pop("gene_name")

    # sum expression for each tissue, across the DB TFs
    barplot_sum_df = pd.DataFrame(expr_df.sum(axis=0), columns=["expression"])
    barplot_sum_df.sort_values("expression", ascending=False, inplace=True)
    print(barplot_sum_df)

    # assign colors to each of the tissue groups
    clrs = []
    for ylabel in barplot_sum_df.index:
        tissue_groups_color = "white"
        
        if ylabel in FANTOM_tissue_categories_d:
            tissue_group = FANTOM_tissue_categories_d[ylabel]

            if tissue_group in tissue_category_color_dict:
                tissue_groups_color = tissue_category_color_dict[tissue_group]

        else:
            print(ylabel)
             
        clrs.append(tissue_groups_color)
        
    # build figure
    fig = plt.figure()
    bp = sns.barplot(x=barplot_sum_df.index, y="expression", data = barplot_sum_df, palette=clrs, edgecolor = "none")
    fig = bp.get_figure()
    plt.xticks(rotation=90)
    plt.tight_layout()

    # save figure
    fig_filename = os.path.join(output_figures_FANTOM_clustermap_dn, ".".join([species, "agg_barplot", "all_tissues", "svg"]))
    fig.savefig(fig_filename)
    plt.clf()


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

# FANTOM expression data 
fantom_cage_exp_data_fn = os.path.join(data_raw_FANTOM_dn, "hg38_fair+new_genes_phase1and2_tpm.osc.txt.healthy_samples_jaspar.top_level_tissues.parquet")
fantom_cage_exp_data_df = pd.read_parquet(fantom_cage_exp_data_fn)

# number of FANTOM tissues to include in heatmap; after summing expression of
# target TF genes in all FANTOM tissues, the top number to include
tissues_to_include = 100

# load table of FANTOM tissues assigned to general categories
FANTOM_tissue_categories_fn = os.path.join(data_raw_FANTOM_dn, "FANTOM.grouped_tissue_categories.tsv")
FANTOM_tissue_categories_df = pd.read_csv(FANTOM_tissue_categories_fn, sep="\t")
FANTOM_tissue_categories_df.fillna("", inplace=True)
FANTOM_tissue_categories_d = dict(zip(FANTOM_tissue_categories_df.grouped_tissue, FANTOM_tissue_categories_df.category))

# color palette for FANTOM tissue categories
tissue_category_color_dict = {"brain, eye & retina":sns.color_palette("Paired")[3],
                              "gastrointestinal":sns.color_palette("Paired")[8],
                              "urinary":sns.color_palette("Paired")[6],                              
                              "muscle":sns.color_palette("Paired")[11],
                              "female/male":sns.color_palette("Paired")[5],
                              "fetal growth":sns.color_palette("Paired")[4],
                              "skin & connective":sns.color_palette("Paired")[7],
                              "epithelial":sns.color_palette("Paired")[2],
                              "adipose & soft":sns.color_palette("Paired")[10],                              
                              "respiratory":sns.color_palette("Paired")[0],
                              "developmental":sns.color_palette("GnBu_d")[0],
                              "blood & immune":sns.color_palette("Paired")[9],
                              "endocrine & exocrine":sns.color_palette("GnBu_d")[4],
                              "circulatory":sns.color_palette("Paired")[1]}

################################################################################
# Execution ####################################################################
################################################################################

# Genes associated to various tissues ##########################################
# ProteinAtlas gene tissue enrichment
protein_atlas_data_filename = os.path.join(data_raw_ProteinAtlas_dn, "proteinatlas.21.0.parquet")
protein_atlas_data = pd.read_parquet(protein_atlas_data_filename)
protein_atlas_data.fillna("", inplace=True)
protein_atlas_tissue_genes_dict = {}

# ProteinAtlas enriched tissue genes
tissue_categories_dict = get_protein_atlas_enriched_tissues_genes()

# ProteinAtlas enriched cell type genes             
get_protein_atlas_enriched_celltypes_genes()

# ProteinAtlas enriched brain cell type genes
get_protein_atlas_enriched_brain_celltypes_genes()

# ProteinAtlas enriched blood cell type genes
get_protein_atlas_enriched_blood_celltypes_genes()

# ProteinAtlas enriched cell line genes
get_protein_atlas_enriched_cellline_genes()

# DICE immune genes
DICE_immunity_genes = get_DICE_immune_genes()
protein_atlas_tissue_genes_dict["blood & immune"]+= DICE_immunity_genes

# developmental (BOX) genes
box_genes_ls = get_box_genes()
protein_atlas_tissue_genes_dict["developmental"] = box_genes_ls

unique_db_TFs_d = get_unique_db_TFs()

for species in ["human", "neand"]:

    # Expression of target genes ###################################################
    target_tfs_cleaned = unique_db_TFs_d[species]
    target_tfs_fantom_cage_exp_data_df = fantom_cage_exp_data_df[fantom_cage_exp_data_df['gene_name'].isin(target_tfs_cleaned)]

    if len(target_tfs_fantom_cage_exp_data_df):

        # restrict to top X tissues with highest aggregate expression among DB TFs
        expr_df = fantom_cage_exp_data_df[fantom_cage_exp_data_df['gene_name'].isin(target_tfs_cleaned)]
        expr_df.pop("gene_name")
        expr_sum_df = pd.DataFrame(expr_df.sum(axis=0), columns=["expression"])
        expr_sum_df.sort_values("expression", ascending=False, inplace=True)
        top_tissues_ls = list(expr_sum_df.index[:tissues_to_include-1])
        drop_tissues_ls = list(expr_sum_df.index[tissues_to_include-1:])
        target_tfs_fantom_cage_exp_data_df.drop(drop_tissues_ls, axis=1, inplace=True)

        # Figures ######################################################################
        # calculate heatmap
        method, metric ="ward", "euclidean"
    
        ylabels, xlabels = heatmap(target_tfs_fantom_cage_exp_data_df)

        # plot genes vs tissue categories
        tf_genes_class_colors(xlabels)

        # plot FANTOM tissues vs tissue categories
        agg_expression_barplot()
        agg_expression_barplot_all_tissues()




        
    

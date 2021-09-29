#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 3.7.0 ###########################################################
# Libraries ####################################################################
import os
from pathlib import Path
import pandas as pd

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
def assign_dirs():
    
    """Assign dirnames to global variables, and build dirs if they don't already
    exist.  Original list follows, most up to date list will be found in
    "../dirnames.tsv"
    
    data_raw_dn				../data_raw
    data_raw_ensembl_dn			../data_raw/Ensembl
    data_raw_ABA_brainspan_dn		../data_raw/Allen_Brain_Atlas_brainspan
    data_raw_ABA_cortical_dn		../data_raw/Allen_Brain_Atlas_cortical
    data_raw_catalog_of_changes_dn	../data_raw/CatalogOfChanges
    data_raw_DICE_dn			../data_raw/DICE
    data_raw_FANTOM_dn			../data_raw/FANTOM
    data_raw_GTEx_dn			../data_raw/GTEx
    data_raw_GTRD_dn			../data_raw/GTRD
    data_raw_ProteinAtlas_dn		../data_raw/ProteinAtlas
    data_dn				../data
    data_DBTFs_dn			../data/DBTFs
    data_GTRD_dn			../data/GTRD
    data_TFBS_analyses_dn		../data/TFBS_footprinter_analysis
    data_TFBS_analyses_human_dn		../data/TFBS_footprinter_analysis/human
    data_TFBS_analyses_neand_dn		../data/TFBS_footprinter_analysis/neanderthal
    data_promoter_variants_dn		../data/promoter_variants
    output_dn				../output
    output_figures_dn			../output/figures
    """

    print("assigning global dirname variables")

    dirnames_fn = "../dirnames.tsv"
    dirnames_df = pd.read_csv(dirnames_fn, sep="\t")
    for i, row in dirnames_df.iterrows():
        dirname_var = row['dirname_var']
        rel_path = row['rel_path']
        Path(rel_path).mkdir(exist_ok=True, parents=True)
        globals()[dirname_var] = rel_path


################################################################################
# Initiating Variables #########################################################
################################################################################


################################################################################
# Execution ####################################################################
################################################################################
assign_dirs()

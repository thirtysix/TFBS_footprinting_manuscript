#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 3.8.0 ###########################################################
# Libraries ####################################################################
import os
import csv
import shutil
from assign_dirs import *

################################################################################
# Description/Notes ############################################################
################################################################################
"""
Downloaded files from the human analysis are copied to a neanderthal dir,
and edited to match the variant cataloged in the CatalogOfChanges data.
TFBS_footprinter is then run on these created Neanderthal sequence versions.

Analysis utilizes TFBS_footprinter friendly tables built by previous steps.
TFBS_footprinter downloads data for each analysis via Ensembl REST.
If internet or server connection errors occur, this script can be re-run and
analyses will restart from where they left off.

The complete TFBS_footprinter table can be manually split and TFBS_footprinter
run via command line if desired.
"""

################################################################################
# Base-level Functions #########################################################
################################################################################
def file_to_datalist(data_filename):
    """
    Starting with a filename, import and convert data to a list.
    """

    with open(data_filename, 'r') as data_file:
        csv_reader = csv.reader(data_file, delimiter = ",")
        all_data = list(csv_reader)

    return all_data

################################################################################
# Task-specific Functions ######################################################
################################################################################
def replace_hominin_base(hsa_base, other_hominin_base, other_hominin_results_dir, hsap_sequence_filename):
    """
    If the SNP catalogued on this line is a variant in the other species,
    copy the human dir and results to the other hominin parent folder and change
    the base to the other hominin and then save the sequence file so that it can
    be analyzed for TFBSs.
    """

    conjugate_dict = {"A":"T","T":"A","G":"C","C":"G"}

    # track discrepancies in CatalogOfChanges bases vs downloaded seqs
    non_match_count = 0
    matching_missing_other_base_count = 0

    if hsa_base != other_hominin_base and other_hominin_base!="":
        out_dir = os.path.join(other_hominin_results_dir, result_dir)
        
        # copy the full human result directory and contents to the other hominin results dir
        if os.path.exists(full_hsap_result_dirpath):
            if not os.path.exists(out_dir):
                shutil.copytree(full_hsap_result_dirpath, out_dir)
        else:
            print("Expected human result dir does not exist:", full_hsap_result_dirpath)

        # delete the hsap results files, these are the only files which are not needed
        copied_sorted_clusters_filename = os.path.join(out_dir, "TFBSs_found.sortedclusters.csv")
        copied_cluster_dict_filename = os.path.join(out_dir, "cluster_dict.json")
        if os.path.exists(copied_sorted_clusters_filename):
            os.remove(copied_sorted_clusters_filename)
        if os.path.exists(copied_cluster_dict_filename):
            os.remove(copied_cluster_dict_filename)
        
        # open the copied hsap sequence in the new other hominin dir, change the snp base to the ancestral/other hominin variant
        other_hominin_sequence_filename = os.path.join(out_dir, "alignment_cleaned.fasta")
        if os.path.exists(other_hominin_sequence_filename):
            
            # open the sequence file
            with open(hsap_sequence_filename, 'r') as hsap_sequence_file:
                read_lines = hsap_sequence_file.readlines()
                header = read_lines[0]
                human_sequence = read_lines[1]

                # convert the human sequence to the other hominin_sequence
                other_hominin_seq = human_sequence

                if strand == "+":            
                    if not hsa_base==human_sequence[25]:
                        non_match_count+=1 # track sequences which are not as expected by CatalogOfChanges
                        print( "no base match:", result_dir)
                    other_hominin_seq = list(other_hominin_seq)
                    other_hominin_seq[25] = other_hominin_base
                    other_hominin_seq = "".join(other_hominin_seq)
                    
                # the sequences were retrieved on just the positive strand because both sides are scanned in TFBS_footprinter
                # however, the other hominin variant is given in the negative strand (when strand is negative)
                # therefore we will just replace the target base on the positive strand with the conjugate of the neanderthal negative strand base
                if strand == "-":
                    if not hsa_base==conjugate_dict[human_sequence[25]]:
                        non_match_count+=1 # track sequences which are not as expected by CatalogOfChanges
                        print("no base match:", result_dir)
                    replace_chars = ""
                    other_hominin_seq = list(other_hominin_seq)

                    for char in other_hominin_base:
                        replace_chars+=conjugate_dict[char]
                    other_hominin_seq[25] = replace_chars
                    other_hominin_seq = "".join(other_hominin_seq)

            # write the new other hominin sequence to the other hominin seq file in the other hominin results dir
            # this folder will then be able to be scanned as a other hominin sequence for TFBSs
            with open(other_hominin_sequence_filename, 'w') as other_hominin_sequence_outfilename:
                other_hominin_sequence_outfilename.write(header)
                other_hominin_sequence_outfilename.write(other_hominin_seq)
            
        
        else:
            print (other_hominin_sequence_filename, "doesn't exist!")

    else:
        # track the entries which are not a human vs. other variant
        matching_missing_other_base_count+=1


################################################################################
# Initiating Variables #########################################################
################################################################################
# target comparison and transcript type
primary_comparison = 'hsap2nea_snp'
primary_transcript_type = "protein_coding"
var_type = 'snps'
promoter_upstream_dist = 2500
promoter_downstream_dist = 2500
promoter_overlap_window_size = promoter_upstream_dist + promoter_downstream_dist
window_size = 50

# filenames
# human_results_table_fn contains the dirnames and primary filename of the TFBS_footprinter analysis of the human promoters
# and human/neanderthal/chimp bases at this location
human_results_table_fn = os.path.join(data_promoter_variants_dn, ".".join([primary_comparison, var_type, primary_transcript_type, "results_mut_type_filenames", str(promoter_overlap_window_size), str(window_size), "csv"]))
human_results = file_to_datalist(human_results_table_fn)


################################################################################
# Execution ####################################################################
################################################################################

# iterate through each of the expected human result dirs
for row in human_results[1:]:
    result_dir, result_filename, mut_type, mut_len, strand, Hsap2Nea_snp_status, Hsap2Pan_snp_status, nea2pan_snp_status, hsap2nea2pan_snp_status, hsa_base, neaH_base, neaC_base, pan_base, transcript_biotype = row

    # current human result
    full_hsap_result_dirpath = os.path.join(data_TFBS_analyses_human_results_dn, result_dir)
    hsap_sequence_filename = os.path.join(data_TFBS_analyses_human_results_dn, result_dir, "alignment_cleaned.fasta")
    
    # replace the human base with the ancestral neandertal base
    replace_hominin_base(hsa_base, neaH_base, data_TFBS_analyses_neand_results_dn, hsap_sequence_filename)











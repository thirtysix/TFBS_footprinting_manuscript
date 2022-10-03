#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 3.8.0 ###########################################################
# Libraries ####################################################################
import os
import pprint
from assign_dirs import *
from pyliftover import LiftOver
from operator import itemgetter
import bisect

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import seaborn as sns

import itertools

################################################################################
# Description/Notes ############################################################
################################################################################
"""
Check upstream of the transcript start site ~500/1000/2500/5000 bp.
Check the 5' UTR regions and 3' UTR regions.

Error-check that pyliftover is giving the appropriate new positions (hg18>38).
Consider larger analysis using Denisovan genome as well.

Ensemlb GTF entry Biotypes:
['3prime_overlapping_ncRNA',
 'IG_C_gene',
 'IG_C_pseudogene',
 'IG_D_gene',
 'IG_J_gene',
 'IG_J_pseudogene',
 'IG_V_gene',
 'IG_V_pseudogene',
 'IG_pseudogene',
 'TEC',
 'TR_C_gene',
 'TR_D_gene',
 'TR_J_gene',
 'TR_J_pseudogene',
 'TR_V_gene',
 'TR_V_pseudogene',
 'antisense',
 'bidirectional_promoter_lncRNA',
 'lincRNA',
 'macro_lncRNA',
 'miRNA',
 'misc_RNA',
 'non_coding',
 'non_stop_decay',
 'nonsense_mediated_decay',
 'polymorphic_pseudogene',
 'processed_pseudogene',
 'processed_transcript',
 'protein_coding',
 'pseudogene',
 'rRNA',
 'rRNA_pseudogene',
 'retained_intron',
 'ribozyme',
 'sRNA',
 'scRNA',
 'scaRNA',
 'sense_intronic',
 'sense_overlapping',
 'snRNA',
 'snoRNA',
 'transcribed_processed_pseudogene',
 'transcribed_unitary_pseudogene',
 'transcribed_unprocessed_pseudogene',
 'translated_processed_pseudogene',
 'unitary_pseudogene',
 'unprocessed_pseudogene',
 'vaultRNA']

"""

################################################################################
# Base-level Functions #########################################################
################################################################################


################################################################################
# Task-specific Functions ######################################################
################################################################################
def convert_genomic_coords(chrom, coord):
    """
    Use PyLiftover to convert genomic coordinates from one genome version to
    another.
    """

    converted_tup_list = lo.convert_coordinate(chrom, coord)
    
    if len(converted_tup_list) >0:
        if len(converted_tup_list[0]) == 4:
            converted_coord = converted_tup_list[0][1]

    else:
        converted_coord = None        

    return converted_coord


def fill_promoter_start(row):
    """
    Based on genomic strand, define the promoter start position.
    """

    if row['strand'] == "+":
        promoter_start_int = row['start'] - int(promoter_overlap_window_size/2)
    if row['strand'] == "-":
        promoter_start_int = row['end'] - int(promoter_overlap_window_size/2)

    return promoter_start_int


def fill_promoter_end(row):
    """
    Based on genomic strand, define the promoter end position.
    """

    if row['strand'] == "+":
        promoter_end_int = row['start'] + int(promoter_overlap_window_size/2)
    if row['strand'] == "-":
        promoter_end_int = row['end'] + int(promoter_overlap_window_size/2)

    return promoter_end_int


def fill_TSS(row):
    """
    Based on genomic strand, define the TSS position.
    """
    if row['strand'] == "+":
        TSS=row['start']
    if row['strand'] == "-":
        TSS=row['end']

    return TSS


def parse_gtf():
    """
    Parse the human genome GTF file.
    Extract annotation data for each entry,
    create DF.
    """

    print ("parsing GTF")
    transcript_list = []

    # import as df
    ensgtf_df = pd.read_csv(human_gtf_data_fn,
                            compression='gzip',
                            sep="\t",
                            skiprows=5,
                            header=None,
                            names =
                            ['chrom','evidence','feature','start','end','dot1','strand','dot2','annotations'],
                            usecols=['chrom','feature','start','end','strand','annotations'])

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

    # create dict for quick referencing
    human_gtf_dict = {x[6]:x for x in transcript_list}
    
    # create new df with transcript annotations parsed to columns
    transcript_data_df = pd.DataFrame(transcript_list, columns=['chrom', 'strand', 'start', 'end', 'gene_name', 'gene_id', 'transcript_id', 'transcript_biotype'])

    # create cols for promoter start and end based on strand
    transcript_data_df['promoter_start'] = transcript_data_df.apply(lambda row : fill_promoter_start(row), axis=1)
    transcript_data_df['promoter_end'] = transcript_data_df.apply(lambda row : fill_promoter_end(row), axis=1)
    transcript_data_df['TSS'] = transcript_data_df.apply(lambda row : fill_TSS(row), axis=1)

    return transcript_list, human_gtf_dict, transcript_data_df


def convert_variant_data():
    """
    Convert the hg18 variant data to hg38.
    """

    print("Convert hg18 variant locs to hg38.")
    # copy the hg18 data as template for hg38 data
    hg38_snp_data_df = hg18_snp_data_df.copy()
    hg38_snp_locs = []

    # get all hg18 chromosome names
    hg18_snp_chrs = list(hg18_snp_data_df['HSA_Chr'])
    
    # convert hg18 variant locs to 1-based from 0-based
    hg18_snp_locs = [int(x)+1 for x in hg18_snp_data_df['HSA_Pos']]

    # zip together
    hg18_snp_chrs_locs = zip(hg18_snp_chrs, hg18_snp_locs)

    # iterate through each chrom/loc and convert using liftover
    for hg18_snp_chr_loc in hg18_snp_chrs_locs:
        hg18_snp_chr, hg18_snp_loc = hg18_snp_chr_loc
        hg38_snp_loc = convert_genomic_coords(hg18_snp_chr, hg18_snp_loc)
        hg38_snp_locs.append(hg38_snp_loc)

    # replace list of hg18 locs in hg38 df with hg38 locs
    hg38_snp_data_df['HSA_Pos'] = hg38_snp_locs

    # remove locations which did not have an hg18>hg38 convertible location
    hg38_snp_data_df = hg38_snp_data_df[hg38_snp_data_df['HSA_Pos']!=None]

    return hg38_snp_data_df


def variants_in_promoters():
    """
    For each snp variant identify a single transcript (at most).
    Variant should overlap the 'promoter' region as defined by
    promoter_overlap_window_size variable.
    """
    
    for target_transcript_type in target_transcript_types:
        for var_type in  var_types: # option to address indels later

            # output filenames
            tfbs_footprinter_table_ofn = os.path.join(data_TFBS_analyses_human_dn, ".".join([comparison, var_type, target_transcript_type, "tfbs_footprinter_table", str(promoter_overlap_window_size), str(window_size), "csv"]))
            tfbs_footprinter_table_neand_ofn = os.path.join(data_TFBS_analyses_neand_dn, ".".join([comparison, var_type, target_transcript_type, "tfbs_footprinter_table", str(promoter_overlap_window_size), str(window_size), "csv"]))
            results_mut_type_filenames_ofn = os.path.join(data_promoter_variants_dn, ".".join([comparison, var_type, target_transcript_type, "results_mut_type_filenames", str(promoter_overlap_window_size), str(window_size), "csv"]))

            # do not proceed if output table exists
            if not os.path.exists(results_mut_type_filenames_ofn):
            
                # length of target mutation, limited to snps now
                mut_len = {'snps':1}[var_type] # for indels mutation len would need to be parsed

                # filter Ensembl transcripts to those of the target biotype (e.g., 'protein_coding')
                transcript_data_biotype_df = transcript_data_df[transcript_data_df['transcript_biotype']==target_transcript_type]

                # initiate output file lists
                tfbs_footprinter_ls = []
                results_mut_type_filenames_ls = []

                # have fun, keep track, win
                variants_done_count = 0

                # iterate through each of the chromosomes in the ensembl transcript data
                for ens_chr in transcript_data_biotype_df['chrom'].unique():

                    # ens transcripts ['chrom', 'strand', 'start', 'end', 'gene_name', 'gene_id', 'transcript_id', 'transcript_biotype']
                    transcript_data_biotype_chr_df = transcript_data_biotype_df[transcript_data_biotype_df['chrom']==ens_chr]

                    # variants with this ens chr
                    hg38_snp_chr_data_df = hg38_snp_data_df[(hg38_snp_data_df['HSA_Chr']==ens_chr) & hg38_snp_data_df[comparison]==True]

                    # only proceed if there are variants (snps) on this chromosome
                    if len(hg38_snp_chr_data_df) > 0:

                        # build interval index
                        idx = pd.IntervalIndex.from_arrays(transcript_data_biotype_chr_df['promoter_start'], transcript_data_biotype_chr_df['promoter_end'])

                        # iterate through all of the hg38 snps on this chrom
                        for i, hg38_snp in hg38_snp_chr_data_df.iterrows():
                            hg38_snp_loc = hg38_snp['HSA_Pos']
                            hsap2nea_snp_status = hg38_snp['hsap2nea_snp']
                            hsap2pan_snp_status = hg38_snp['hsap2pan_snp']
                            nea2pan_snp_status = hg38_snp['nea2pan_snp']
                            hsap2nea2pan_snp_status = hg38_snp['hsap2nea2pan_snp']
                            
                            hsa_base = hg38_snp['#HSA_Base']
                            neaH_base = hg38_snp['NEA_BaseH']
                            neaC_base = hg38_snp['NEA_BaseC']
                            pan_base = hg38_snp['PAN_Base']
                            
                            # identify transcripts which occur within the promoter_start and promoter_end
                            overlapping_ens_transcripts_df = transcript_data_biotype_chr_df[idx.contains(hg38_snp_loc)]

                            # find transcript which the variant is closest to the TSS of
                            if len(overlapping_ens_transcripts_df) > 0:
                                overlapping_ens_transcripts_df['abs_dist_to_var'] = abs(overlapping_ens_transcripts_df['TSS'] - hg38_snp_loc)
                                variant_closest_transcript = overlapping_ens_transcripts_df.sort_values(by=['abs_dist_to_var']).iloc[0]

                                # create TFBS footprinter table entry
                                transcript_tss = variant_closest_transcript['TSS']                        
                                var_center = int((hg38_snp_loc + hg38_snp_loc + 1)/2)
                                transcript_id = variant_closest_transcript['transcript_id']
                                transcript_strand = variant_closest_transcript['strand']

                                # for TFBS_footprinter define start and end locations relative to variant
                                if transcript_strand == "+":
                                    tss_relative_start = int((var_center - transcript_tss ) - window_size/2)
                                    tss_relative_end = int((var_center - transcript_tss) + window_size/2)
                                if transcript_strand == "-":
                                    tss_relative_start = int((transcript_tss - var_center) - window_size/2)
                                    tss_relative_end = int((transcript_tss - var_center) + window_size/2)

                                # convert to format friendly for TFBS_footprinter runs
                                tfbs_footprinter_start = tss_relative_start*-1
                                tfbs_footprinter_end = tss_relative_end      

                                # create TFBS_footprinter analysis table entry 
                                tfbs_footprinter_ls.append([transcript_id, "", tfbs_footprinter_start, tfbs_footprinter_end, tfs_to_plot, pval, pval_c])
                                result_dir = "".join([transcript_id+"("+str(tfbs_footprinter_start)+"_"+str(tfbs_footprinter_end)+")"+"_"+pval])
                                result_filename = os.path.join(result_dir, "TFBSs_found.sortedclusters.csv")
                                results_mut_type_filenames_ls.append([result_dir, result_filename, var_type, mut_len, transcript_strand, str(hsap2nea_snp_status), str(hsap2pan_snp_status), str(nea2pan_snp_status), str(hsap2nea2pan_snp_status), hsa_base, neaH_base, neaC_base, pan_base, target_transcript_type])

                            # have fun, keep track, win
                            variants_done_count+=1
                            if variants_done_count%100000 == 0:
                                print(comparison, "variants assessed:", variants_done_count)

                # write results dirs table
                results_mut_type_filenames_df = pd.DataFrame(results_mut_type_filenames_ls, columns=["result_dir", "result_filename", "mut_type", "mut_len", "strand", "hsap2nea snp status","hsap2pan snp status", "nea2pan snp status", "hsap2nea2pan snp status", 'hsa_base', 'neaH_base', 'neaC_base', 'pan_base', 'transcript_biotype'])
                results_mut_type_filenames_df.to_csv(results_mut_type_filenames_ofn, sep=",", index=False)

                # write TFBS_footprinter table
                if comparison == primary_comparison and target_transcript_type == primary_transcript_type:
                    tfbs_footprinter_table_df = pd.DataFrame(tfbs_footprinter_ls, columns=["ens_tid", "target_tfs_filename","before TSS","after TSS","tfs to plot", "pval", "pval_c"])
                    tfbs_footprinter_table_df.to_csv(tfbs_footprinter_table_ofn, sep=",", index=False)
                    tfbs_footprinter_table_df.to_csv(tfbs_footprinter_table_neand_ofn, sep=",", index=False)

            
def venn_plot():
    """
    Create a venn diagram plot from the number of variants between species.
    Requires the matplotlib_venn library.
    Target comparison is for snps in protein_coding transcripts promoter region.
    """

    from matplotlib import pyplot as plt
    from matplotlib_venn import venn3, venn3_circles
    import numpy as np

    # specify target variation type and transcript type
    var_type = "snps"
    target_transcript_type = "protein_coding"

    # specify format of outfigure
    fig_outfmt = "svg"

    # figure aesthetics
    color_pal = ['red', 'green', 'blue']
    plt.figure(figsize=(10,10))

    # labels correspond to default comparison order
    labels = ["Homo sapiens", "Neanderthal", "Chimpanzee"]

    # get counts of variants of type 'var_type' and in transcript type of 'target_transcript_type', for each of the comparisons made
    comparison_variants_count_d = {}

    for comparison in comparisons:
        # load variants for the comparison of these two species
        results_mut_type_filenames_ofn = os.path.join(data_promoter_variants_dn, ".".join([comparison, var_type, target_transcript_type, "results_mut_type_filenames", str(promoter_overlap_window_size), str(window_size), "csv"]))
        results_mut_type_filenames_df = pd.read_csv(results_mut_type_filenames_ofn, sep=",")
        comparison_variants_count_d[comparison] = len(results_mut_type_filenames_df)

    # format of subset sizes for 'venn3': (Abc, aBc, ABc, abC, AbC, aBC, ABC)
    set_default_size = 200000
    subset_sizes = (set_default_size, set_default_size, comparison_variants_count_d[comparisons[0]], set_default_size, comparison_variants_count_d[comparisons[1]], comparison_variants_count_d[comparisons[2]], comparison_variants_count_d[comparisons[3]])
    v = venn3(subsets = subset_sizes, set_labels = labels, set_colors=color_pal, alpha=0.5)

    # save the figure to the figure dir
    plt.savefig(os.path.join(output_figures_dn, ".".join([target_transcript_type, "venn_plot", fig_outfmt])))

    
################################################################################
# Initiating Variables #########################################################
################################################################################
# pyliftover querty/target genomes
lo = LiftOver('hg18', 'hg38')

# desired promoter boundaries for analysis of all human promoters
promoter_upstream_dist = 2500
promoter_downstream_dist = 2500
pval = str(1)
pval_c = str(0.1)
tfs_to_plot = str(10)

# target data filenames
##indel_fn = os.path.join(data_raw_catalog_of_changes_dn, "combined_indel_anno.har.filter.tsv")
snp_fn = os.path.join(data_raw_catalog_of_changes_dn, "combined_SNP_anno.ns.har.filter.tsv")
human_gtf_data_fn = [os.path.join(data_raw_ensembl_dn, x) for x in os.listdir(data_raw_ensembl_dn) if "chr.gtf" in x and ".gz" in x][0]

# vars
primary_comparison = 'hsap2nea_snp' # the target species comparison
var_types = ['snps'] # variant types to use in comparisons
primary_transcript_type = "protein_coding" # the transcript types to use in analysis

# thresholds/windows
window_size = 50
promoter_overlap_window_size = promoter_upstream_dist + promoter_downstream_dist


################################################################################
# locate variants occurring in promoters #######################################
################################################################################
print("Analysis takes ~{} minutes to run to completion.".format(75))

# load hg18 snp data
hg18_snp_data_df = pd.read_csv(snp_fn, sep="\t")
print("A total of {} comparison variants will be assessed for overlap with promoter regions.".format(len(hg18_snp_data_df)))

# bools for presence of snps
hg18_snp_data_df['mut_len'] = [1]* len(hg18_snp_data_df)
hg18_snp_data_df['hsap2nea_snp'] = (hg18_snp_data_df['#HSA_Base']!=hg18_snp_data_df['NEA_BaseH']) | (hg18_snp_data_df['#HSA_Base']!=hg18_snp_data_df['NEA_BaseC'])
hg18_snp_data_df['hsap2pan_snp'] = hg18_snp_data_df['#HSA_Base']!=hg18_snp_data_df['PAN_Base']
hg18_snp_data_df['nea2pan_snp'] = (hg18_snp_data_df['PAN_Base']!=hg18_snp_data_df['NEA_BaseH']) | (hg18_snp_data_df['PAN_Base']!=hg18_snp_data_df['NEA_BaseC'])
hg18_snp_data_df['hsap2nea2pan_snp'] = ((hg18_snp_data_df['PAN_Base']!=hg18_snp_data_df['NEA_BaseH']) & (hg18_snp_data_df['#HSA_Base']!=hg18_snp_data_df['PAN_Base']) &
                                        (hg18_snp_data_df['#HSA_Base']!=hg18_snp_data_df['NEA_BaseH'])) | ((hg18_snp_data_df['PAN_Base']!=hg18_snp_data_df['NEA_BaseC']) &
                                         (hg18_snp_data_df['#HSA_Base']!=hg18_snp_data_df['PAN_Base']) & (hg18_snp_data_df['#HSA_Base']!=hg18_snp_data_df['NEA_BaseC']))

# convert the all variant locations from hg18 to hg38
hg38_snp_data_df_orig = convert_variant_data()

# parse Ensembl transcripts and create dictionary
transcript_list, transcript_data_dict, transcript_data_df = parse_gtf()

# Analysis of protein coding snps in all comparisons
# Focus on protein coding transcripts
comparisons = ['hsap2nea_snp', 'nea2pan_snp', 'hsap2pan_snp', 'hsap2nea2pan_snp']
target_transcript_types = ["protein_coding"]
for comparison in comparisons:
    print("Comparing:", comparison)
    
    # limit to target comparison
    # generate a copy of the original data and slice to this comparison
    hg38_snp_data_df = hg38_snp_data_df_orig.copy()
    hg38_snp_data_df = hg38_snp_data_df[hg38_snp_data_df[comparison] == True]

    # identify transcripts whose promoters contain a variant, and create TFBS_footprinter table for analysis
    variants_in_promoters()


################################################################################
# promoter occupancy by transcript type ########################################
################################################################################
"""
For each of the desired transcript biotypes, identify the frequency that
genetic variants occur at each location in the defined promoter region --
e.g., +/- 2,500 bp.
To-do: plot the results as histogram of generic promoter.
"""
##
### Analysis of protein coding snps in comparisons
### Analysis of all transcript type snps in primary_comparison
##comparisons = [primary_comparison]
##target_transcript_types = transcript_data_df['transcript_biotype'].unique()
##for comparison in comparisons:
##    print("Comparing:", comparison)
##
##    # limit to target comparison
##    # generate a copy of the original data and slice to this comparison
##    hg38_snp_data_df = hg38_snp_data_df_orig.copy()
##    hg38_snp_data_df = hg38_snp_data_df[hg38_snp_data_df[comparison] == True]
##
##    # identify transcripts whose promoters contain a variant, and create TFBS_footprinter table for analysis
##    variants_in_promoters()
##
##
##var_type = 'snps'
##comparison = 'hsap2nea_snp'
##target_transcript_types = transcript_data_df['transcript_biotype'].unique()
##
##occupancy_ls = []
##transcript_data_df['range'] = [range(x, x+promoter_upstream_dist+promoter_downstream_dist) for x in transcript_data_df['promoter_start'].tolist()]
##for target_transcript_type in target_transcript_types:
##    print(target_transcript_type)
##    results_fn = os.path.join(data_promoter_variants_dn, ".".join([comparison, var_type, target_transcript_type, "results_mut_type_filenames", str(promoter_overlap_window_size), str(window_size), "csv"]))    
##    results_df = pd.read_csv(results_fn)
##
##    target_transcript_type_promoter_breadth = 0
##
##    transcript_data_biotype_df = transcript_data_df[transcript_data_df['transcript_biotype']==target_transcript_type]
##    for ens_chrom in transcript_data_df['chrom'].unique():
##        chr_complete_ranges = []
##        transcript_data_biotype_chrom_df = transcript_data_biotype_df[transcript_data_biotype_df['chrom']==ens_chrom]
##
##        chr_unique_nts = set(itertools.chain.from_iterable(transcript_data_biotype_chrom_df['range'].tolist()))
##
##        target_transcript_type_promoter_breadth += len(chr_unique_nts)
##
##    occupancy = len(results_df)/target_transcript_type_promoter_breadth
##    occupancy_ls.append([target_transcript_type, len(results_df), len(transcript_data_biotype_df), target_transcript_type_promoter_breadth, occupancy])
##    occupancy_df = pd.DataFrame(occupancy_ls, columns=['target_transcript_type', 'variant count', 'transcript type count', 'promoter nts', 'occupancy'])
##    print(target_transcript_type, len(results_df), target_transcript_type_promoter_breadth, occupancy)
##        

################################################################################
# Venn diagram of snps in protein coding transcripts, between species ##########
################################################################################
##venn_plot()








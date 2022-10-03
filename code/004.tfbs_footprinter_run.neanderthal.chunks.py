#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 3.8.0 ###########################################################
# Libraries ####################################################################
import os
import subprocess
import pandas as pd

from assign_dirs import *

################################################################################
# Description/Notes ############################################################
################################################################################
"""
A copy of the output of the TFBS_footprinter analysis of human promoters has
been made and the human sequences modified to ancestral Neanderthal seqs.
This script will now run TFBS_footprinter on the modified sequences and will not
need to download seq data from Ensembl REST.

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


################################################################################
# Task-specific Functions ######################################################
################################################################################


################################################################################
# Initiating Variables #########################################################
################################################################################
# TFBS_footprinter outputs results in dir it is run from, so change there
os.chdir(data_TFBS_analyses_neand_dn) 

# var
thread_c = 2 # number of instances to run at once, one per core (est. 1-5GB memory requirement each)
primary_comparison = 'hsap2nea_snp'
primary_transcript_type = "protein_coding"
var_type = 'snps'
promoter_upstream_dist = 2500
promoter_downstream_dist = 2500
promoter_overlap_window_size = promoter_upstream_dist + promoter_downstream_dist
window_size = 50

################################################################################
# Execution ####################################################################
################################################################################
# Neanderthal run
print("Analysis will run with {} instances of TFBS_footprinter.".format(thread_c))
print("Each instance will use ~1-5GB of memory.  With 5 instances running, the total analysis can take 5+ hours.")
try:
    tfbs_footprinter_table_fn = ".".join([primary_comparison, var_type, primary_transcript_type, "tfbs_footprinter_table", str(promoter_overlap_window_size), str(window_size), "csv"])
    tfbs_footprinter_table_size = len(pd.read_csv(tfbs_footprinter_table_fn))

    chunksize=int(tfbs_footprinter_table_size/thread_c)+1

    procs = []

    try:
        for i, chunk in enumerate(pd.read_csv(tfbs_footprinter_table_fn, chunksize=chunksize)):
            tfbs_footprinter_table_bn, ext = os.path.splitext(tfbs_footprinter_table_fn)
            tfbs_footprinter_table_chunk_fn = ".".join([tfbs_footprinter_table_bn, 'chunk{}.csv'.format(i)])
            chunk.to_csv(tfbs_footprinter_table_chunk_fn, index=False)

            args = ["tfbs_footprinter3", "-t", tfbs_footprinter_table_chunk_fn, "-no"]
            tfbs_run_str = " ".join(args)
            print(tfbs_run_str)

            procs.append(subprocess.Popen(args))

        for p in procs:
            p.wait()

    except:
        print("Error running tfbs_footprinter.  Check log file")
        
except:
    print("No tfbs_footprinter table file found")








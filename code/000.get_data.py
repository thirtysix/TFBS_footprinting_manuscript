#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 3.7.0 ###########################################################
# Libraries ####################################################################
import os
from pathlib import Path
from assign_dirs import *

import pandas as pd
import gzip
import shutil
import tarfile

import ftplib
import urllib.request
import ssl

################################################################################
# Description/Notes ############################################################
################################################################################
"""
Get variants.
Get Ensembl GTF.
"""

################################################################################
# Base-level Functions #########################################################
################################################################################


################################################################################
# Task-specific Functions ######################################################
################################################################################
def get_ens_gtf():
    """Ensembl GTF updates frequently.  Get the most current version."""

    # ftp var
    ens_ftp = "ftp.ensembl.org"
    ens_curr_release_dn = "pub/current_gtf/homo_sapiens/"

    # connect to Ensembl FTP
    ftp = ftplib.FTP(ens_ftp)
    ftp.login() # anonymous
    ftp.cwd(ens_curr_release_dn)
    ftp_fns = ftp.nlst()

    # identify the filename of the current release
    curr_release_chr_gtf_fn = [x for x in ftp_fns if 'chr.gtf.gz' in x][0] # assumes only one '.chr'
    curr_release_chr_gtf_url = "https://"+os.path.join(ens_ftp, ens_curr_release_dn, curr_release_chr_gtf_fn)
    curr_release_chr_gtf_ofn = os.path.join(data_raw_ensembl_dn, curr_release_chr_gtf_fn)

    # download archive file if it doesn't exist
    if not os.path.exists(curr_release_chr_gtf_ofn):
        try:
            ssl._create_default_https_context = ssl._create_unverified_context # avoid 'SSL: CERTIFICATE_VERIFY_FAILED' error
            urllib.request.urlretrieve(curr_release_chr_gtf_url, curr_release_chr_gtf_ofn)
        except:
            print("Download of {} data was unsuccessful".format('Ensembl GTF'))

    else:
        print("{} data already exists".format('Ensembl GTF'))

def get_variants():
    """Get variants: Modern human vs. Neanderthal vs. Chimp"""

    # ftp var
    variants_static_url = "https://"+"ftp.ebi.ac.uk/pub/databases/ensembl/neandertal/CatalogOfChanges.tgz"
    variants_ofn = os.path.join(data_raw_catalog_of_changes_dn, os.path.basename(variants_static_url))
    print(variants_ofn)

    # download archive file if it doesn't exist
    if not os.path.exists(variants_ofn):
        try:
            ssl._create_default_https_context = ssl._create_unverified_context # avoid 'SSL: CERTIFICATE_VERIFY_FAILED' error
            urllib.request.urlretrieve(variants_static_url, variants_ofn)
        except:
            print("Download of {} data was unsuccessful".format('variant'))

    else:
        print("{} data already exists".format('variant'))
    # extract file
    tar = tarfile.open(variants_ofn)
    for item in tar:
        tar.extract(item, os.path.dirname(variants_ofn))


################################################################################
# Initiating Variables #########################################################
################################################################################


################################################################################
# Execution ####################################################################
################################################################################
print("Will download Ensembl GTF data and Neanderthal/Human/Chimp variant data files, each ~<50MB. Variant data will be unzipped, and its contents are ~250MB.")
get_ens_gtf()
get_variants()











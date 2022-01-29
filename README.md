

## TFBS_footprinter manuscript analyses

# 1. Background
This repository provides Python scripts used to perform analyses and generate images for the pre-print article found at [biorxiv](https://www.biorxiv.org/content/10.1101/2020.09.04.282749v2).
Implements the transcription factor binding site prediction tool [TFBS_footprinter3](https://github.com/thirtysix/TFBS_footprinting3) with additional documentation available at [ReadTheDocs](https://tfbs-footprinting.readthedocs.io/en/latest/index.html)

---
---
# 2. Instructions
We will give instructions to set up a miniconda environment to contain the dependencies needed and to run the manuscript's analysis scripts.

## 2.1 Install miniconda
Miniconda is the barebones version of the larger Conda package. We will use this so that we choose only the dependencies that are needed and therefore reduce the installation size and time. The Miniconda installation instructions are here: [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

## 2.2 Create environment and activate
```
    $ conda create -n tfbs_footprinter_analyses python=3.7
    $ conda activate tfbs_footprinter_analyses
```

## 2.3 Install tfbs_footprinter3
```
    $ conda install -c thirtysix tfbs-footprinting3 -c conda-forge
```

## 2.4 Install analyses dependencies
```
    $ conda install -c anaconda seaborn
    $ conda install -c bioconda pyliftover
```

---
---
# 3. Script files
### __000.get_data.py__
**DESCRIPTION:** 
Will download Ensembl GTF data and Neanderthal/Human/Chimp variant data files, each ~<50MB. Variant data will be unzipped, and its contents are ~250MB.
__INPUTS:__ N/A
__OUTPUTS:__
* (/data_raw/Ensembl) Ensembl GTF data.
* (/data_raw/CatalogOfChanges) Neanderthal/Human/Chimp variant data files.

---
### __001.neanderthal_variants.py__
**DESCRIPTION:** 
Identify modern human vs. Neanderthal single-nucleotide polymorphisms (SNPs) which occur within 2,500 bp (upstream/downstream) of a human protein-coding transcript transcription start site.
__INPUTS:__
- (/data_raw/Ensembl) Ensembl GTF data.
- (/data_raw/CatalogOfChanges) Neanderthal/Human/Chimp variant data files.

__OUTPUTS:__
* (/data/TFBS_footprinter_analysis/human) TFBS_footprinter run friendly table 
* (/data/promoter_variants) Summary table of MH vs. Neanderthal variants (e.g., SNPs) inside defined promoter region (e.g., +/- 2,500 bp of TSSs) of target transcript type  (e.g., protein-coding).

---
### __002.tfbs_footprinter_run.human.chunks.py__
**DESCRIPTION:** 
Run a transcription factor binding site (TFBS) prediction using tfbs_footprinter on the human version of the human vs. Neanderthal SNPs, which have been identified within 2,500 bp of a human protein-coding transcript transcription start site (TSS).
__INPUTS:__
* (/data/TFBS_footprinter_analysis/human) TFBS_footprinter run friendly table.

__OUTPUTS:__
* (/data/TFBS_footprinter_analysis/human/tfbs_results) TFBS_footprinter run human results.

---
### __003.human2neandertal_chimp.py__
**DESCRIPTION:** 
Downloaded files from the human analysis are copied to a neanderthal dir, and edited to match the variant cataloged in the CatalogOfChanges data. TFBS_footprinter will in the next step be run on these created Neanderthal sequence versions.
__INPUTS:__
* (/data/TFBS_footprinter_analysis/human/tfbs_results) TFBS_footprinter run results.

__OUTPUTS:__
* (/data/TFBS_footprinter_analysis/neanderthal/tfbs_results) TFBS_footprinter downloaded data copied to Neanderthal directory, and modified based on SNP data.

---
### __004.tfbs_footprinter_run.neanderthal.chunks.py__
**DESCRIPTION:** 
Run a TFBS prediction using tfbs_footprinter on the Neanderthal version of all human vs. Neanderthal SNPs, which have been identified within 2,500 bp of a human protein-coding transcript transcription start site.
__INPUTS:__
* (/data/TFBS_footprinter_analysis/neanderthal) TFBS_footprinter run friendly table.

__OUTPUTS:__
* (/data/TFBS_footprinter_analysis/neanderthal/tfbs_results) TFBS_footprinter run Neanderthal results.

---
### __005.tf_frame_score_changes.py__
**DESCRIPTION:** 
Match and analyze differences in predicted TFBSs at modern human vs. Neanderthal SNPs. Produce a list of those TFs which have a statistically significant difference in binding between species; differentially binding TFs (DBTFs).
__INPUTS:__
* (/data/TFBS_footprinter_analysis/human/tfbs_results) TFBS_footprinter run human results.
* (/data/TFBS_footprinter_analysis/neanderthal/tfbs_results) TFBS_footprinter run Neanderthal results.

__OUTPUTS:__
* (/data/DBTFs) Tables of statistical significant differentially binding TFs, for each species.  Tables cataloging all significant TFBS_footprinter putative TFBSs overlapping each SNP, for each species.

---
### __006.clustermap.py__
**DESCRIPTION:** 
Using the list DBTFs, extract expression data from the FANTOM dataset and perform cluster analysis and produce clustermap figure. Using data from The Protein Atlas, identify which DBTFs have tissue-specific expression and produce matrix figure.
__INPUTS:__
* (/data_raw/FANTOM) Bulk RNA-Seq gene expression data from the FANTOM 5 project.  Custom parquet file has been created to limit the FANTOM 5 data to expression of only those TF genes described in the JASPAR database.  Likewise, disease tissues have been removed and remaining healthy tissues grouped to top-level umbrella tissues (\~220MB \*.txt.gz to \~1.1MB \*.parquet)
* (/data_raw/ProteinAtlas) ProteinAtlas annotations of gene to tissue types.
__OUTPUTS:__
* (/output/figures/FANTOM_clustermap) Clustermap depicting correlations of expression for 97 DB TFs in 100 tissues. 

---
### __999.scanpy_run.py__
**DESCRIPTION:** 
Analyze single cell RNA-Seq data from human cortex by [Hodge et al., Conserved cell types with divergent features in human versus mouse cortex](https://www.nature.com/articles/s41586-019-1506-7)
__INPUTS:__
* (/data/TFBS_footprinter_analysis/human/tfbs_results) TFBS_footprinter run human results.
* (/data_raw/Allen_Brain_Atlas_cortical) Single-cell RNA-Seq data from human cortex.

__OUTPUTS:__
* (/output/figures/cortical_single_cell) Force directed graph plots of cell type clusters and matrix plot figures and tables.
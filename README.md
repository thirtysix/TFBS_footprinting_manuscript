

## TFBS_footprinter manuscript analyses

# 1. Background
This repository provides Python scripts used to perform analyses and generate images for the pre-print article found at [biorxiv](https://www.biorxiv.org/content/10.1101/2020.09.04.282749v2).
Implements the transcription factor binding site prediction tool [TFBS_footprinter3](https://github.com/thirtysix/TFBS_footprinting3) with additional documentation available at [ReadTheDocs](https://tfbs-footprinting.readthedocs.io/en/latest/index.html)

---
---
# 2. Instructions - Setup
We will give instructions to set up a miniconda environment to contain the dependencies needed and to run the manuscript's analysis scripts.


## 2.1 Install miniconda
Miniconda is the barebones version of the larger Conda package. We will use this so that we choose only the dependencies that are needed and therefore reduce the installation size and time. The Miniconda installation instructions are here: [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

## 2.2 Create environment and activate
```
    $ conda create -n tfbs_footprinter_analyses python=3.8
    $ conda activate tfbs_footprinter_analyses
```

## 2.3 Install tfbs_footprinter3
```
    $ conda install -c thirtysix tfbs-footprinting3 -c conda-forge
```

## 2.4 Install analyses dependencies
```
    $ conda install -c anaconda seaborn Pillow
    $ conda install -c bioconda pyliftover
    $ conda install -c conda-forge pyarrow=7.0.0 fastparquet
    $ pip install wget
```

## 2.5 Download data

Move to our code directory.
```
    $ cd ./code
```

Download Ensembl and Neanderthal variant data.
```
    $ python3 000.get_data.py
```



# 3. Instructions - Promoterome analysis
## 3.0 - Bypass 
The tasks done in Promoterome analysis section have significant computational, network (80,000 Ensembl REST calls), disk space (300,000 files; ~2GB), and time (10-30 hours) requirements.
These can be bypassed by directly downloading complete results files (two ~60MB files) from the OSF.io repository and placing them in the relevant subfolder within (/data/TFBS_footprinter_analysis)
 - [modern human results parquet file](https://osf.io/download/p8umw/) to (/data/TFBS_footprinter_analysis/human)
 - [Neanderthal results parquet file](https://osf.io/download/w56gr/) to (/data/TFBS_footprinter_analysis/neanderthal)



## 3.1 Identify variants 
Those that occur in promoters of protein coding transcripts (+/- 2,500 bp from TSS).
```
    $ python3 001.neanderthal_variants.py
```

## 3.2 Run TFBS_footprinter3 - Human
Targets are the human verion of variant locations (\~22,000 positions/transcript promoters).
 - Generates ~154,000 files; Time and RAM requirement depends on number of threads run.
 - Number of simultaneous threads can be set in script file (\~1 to 5GB required for each thread), default is 2 threads.
 - Will download sequences (\~22,000) centered on each variant (e.g., 50bp) using Ensembl REST
 - Score with all JASPAR TF models.

```
    $ python3 002.tfbs_footprinter_run.human.chunks.py
```

## 3.3 Copy downloaded human promoter sequences
Create the Neanderthal sequence versions by substituting in the Neanderthal variant to the copied human sequence.
 
```
    $ python3 003.human2neandertal.py
```



## 3.4 Run TFBS_footprinter3 - Neanderthal
Targets are the Neanderthal verion of variant locations (\~22,000 positions/transcript promoters).
 - Generates ~154,000 files; Time and RAM requirement depends on number of threads run.
 - Number of simultaneous threads can be set in script file (\~1 to 5GB required for each thread), default is 2 threads.
 - Sequences do not need to be redownloaded, as they are copied from modern human version and modified in last step.
 - Score with all JASPAR TF models.

```
    $ python3 004.tfbs_footprinter_run.neanderthal.chunks.py
```



# 4. Instructions - Analysis of the differentially binding TFs (DB TFs)
Identify the DB TFs and perform various analyses of DB TF expression in secondary datasets.

## 4.1 Compare TF binding affinities
Identify dfferences in TF binding scores between modern human vs. Neanderthal variants.
 - Identify significant TF binding sites, with Benjimini-Hochbergy correction filtered at 0.01 FDR.
 - Compare significant results between modern human and Neanderthal with Wilcoxon rank statistical test with Benjimini-Hochbergy correction filtered at 0.1 FDR.

```
    $ python3 005.tf_frame_score_changes.py
```


## 4.2 Generate cluster/heatmap - FANTOM
Of differentially binding TFs gene expression, using all healthy whole body tissues in FANTOM dataset. A truncated (to only TF genes) and compressed (to parquet format) FANTOM dataset file is provided as part of this repository (~400 MB to ~1 MB).
```
    $ python3 006.clustermap.fantom_expression_tf_genes.py
```


## 4.3 Generate cluster/heatmap and expression module plots - Allen Brain Atlas
Of differentially binding TFs gene expression, using all healthy brain tissues in Allen Brain Atlas dataset. A truncated (to only TF genes) and compressed (to parquet format) Allen Brain Atlas dataset file is provided as part of this repository (~180 MB to ~3 MB).
```
    $ python3 007.allen_brain_atlas.py
```


## 4.4 Analysis of single cell data from cortical regions
See script description below regarding \~5GB of data that needs to be manually downloaded from [https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq](https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq)
Additional Python libraries need to be installed to process the single cell data.
```
    $ conda install -c conda-forge python-igraph
    $ pip install scanpy==1.8.2
    $ conda install -c conda-forge louvain umap-learn fa2
```

Run the analysis script:
```
    $ python3 999.scanpy_run.py
```

# 5. Instructions - Identify canonical promoter boundaries with mapping of ChIP-Seq peaks
## 5.1 Map of ChIP-Seq data to promoter regions 
From 21,988 ChIP-Seq experiments. Determine occupancy and optimal boundaries for searching for TF binding events. See script description below regarding \~8GB of data that needs to be manually downloaded from [http://gtrd.biouml.org:8888/downloads/current/intervals/chip-seq](http://gtrd.biouml.org:8888/downloads/current/intervals/chip-seq/).

```
    $ python3 999.gtrd_peaking.3prime_5prime.py
```



---
---
# 6. Script file descriptions

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
* (/data/promoter_variants) Summary table of MH vs. Neanderthal variants (e.g., SNPs) inside defined promoter region (e.g., +/- 2,500 bp of TSSs) of target transcript type (e.g., protein-coding).





---
### __002.tfbs_footprinter_run.human.chunks.py__
**DESCRIPTION:** 
Run a transcription factor binding site (TFBS) prediction using tfbs_footprinter on the human version of the human vs. Neanderthal SNPs, which have been identified within 2,500 bp of a human protein-coding transcript transcription start site (TSS). The number of threads to run simultaneously can be set in the script file ('thread_c=', default is 2); each thread will use 1-5GB of RAM. With 5 instances running, the total analysis can take 20+ hours. The primary time constraint is downloading of sequence and transcript data for each transcript from Ensembl REST.

__INPUTS:__
* (/data/TFBS_footprinter_analysis/human) TFBS_footprinter run friendly table.

__OUTPUTS:__
* (/data/TFBS_footprinter_analysis/human/tfbs_results) TFBS_footprinter run human results.





---
### __003.human2neandertal.py__
**DESCRIPTION:** 
Downloaded files from the human analysis are copied to a neanderthal dir, and edited to match the variant cataloged in the CatalogOfChanges data. TFBS_footprinter will in the next step be run on these created Neanderthal sequence versions.

__INPUTS:__
* (/data/TFBS_footprinter_analysis/human/tfbs_results) TFBS_footprinter run results.

__OUTPUTS:__
* (/data/TFBS_footprinter_analysis/neanderthal/tfbs_results) TFBS_footprinter downloaded data copied to Neanderthal directory, and modified based on SNP data.





---
### __004.tfbs_footprinter_run.neanderthal.chunks.py__
**DESCRIPTION:** 
Run a TFBS prediction using tfbs_footprinter on the Neanderthal version of all human vs. Neanderthal SNPs, which have been identified within 2,500 bp of a human protein-coding transcript transcription start site. The number of threads to run simultaneously can be set in the script file ('thread_c=', default is 2); each thread will use 1-5GB of RAM. With 5 instances running, the total analysis can take 5+ hours. The primary time constraint observed in the human analysis is not applicable here as the downloaded human data is copied and modified directly, and thus the Neanderthal analysis only involves calculation of TFBS affinity.

__INPUTS:__
* (/data/TFBS_footprinter_analysis/neanderthal) TFBS_footprinter run friendly table.

__OUTPUTS:__
* (/data/TFBS_footprinter_analysis/neanderthal/tfbs_results) TFBS_footprinter run Neanderthal results.





---
### __005.tf_frame_score_changes.py__
**DESCRIPTION:** 
Match and analyze differences in predicted TFBSs at modern human vs. Neanderthal SNPs. Produce a list of those TFs which have a statistically significant difference in binding between species; differentially binding TFs (DB TFs).

__INPUTS:__
* (/data/TFBS_footprinter_analysis/human/tfbs_results) TFBS_footprinter run human results.
* (/data/TFBS_footprinter_analysis/neanderthal/tfbs_results) TFBS_footprinter run Neanderthal results.

__OUTPUTS:__
* (/data/DBTFs) Tables of statistical significant differentially binding TFs, for each species. Tables cataloging all significant TFBS_footprinter putative TFBSs overlapping each SNP, for each species.





---
### __006.clustermap.py__
**DESCRIPTION:** 
Using the list of DB TFs, extract expression data from the FANTOM dataset and perform cluster analysis and produce clustermap figure. Using data from The Protein Atlas, identify which DB TFs have tissue-specific expression and produce matrix figure.

__INPUTS:__
* (/data_raw/FANTOM) Bulk RNA-Seq gene expression data from the FANTOM 5 project. Custom parquet file has been created to limit the FANTOM 5 data to expression of only those TF genes described in the JASPAR database. Likewise, disease tissues have been removed and remaining healthy tissues grouped to top-level umbrella tissues (\~220MB \*.txt.gz to \~1.1MB \*.parquet)
* (/data_raw/ProteinAtlas) ProteinAtlas annotations of gene to tissue types.

__OUTPUTS:__
* (/output/figures/FANTOM_clustermap) Clustermap depicting correlations of expression for 97 DB TFs in 100 tissues. 





---
### __999.scanpy_run.py__
**DESCRIPTION (IMPORTANT: Run can require up to 70+GB of RAM to run to completion):** 
Analyze single cell RNA-Seq data from human cortex by [Hodge et al., Conserved cell types with divergent features in human versus mouse cortex](https://www.nature.com/articles/s41586-019-1506-7). In order to run this analysis the single cell data must be downloaded manually and placed in the '/data_raw/Allen_Brain_Atlas_cortical' directory. This data can be found at [https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq](https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq), and the relevant three files are:
 * [Gene expression matrix (~5 GB)](https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv),
 * [Table of cell metadata](https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv), and 
 * [2D coordinates](https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/tsne.csv).

__INPUTS:__
* (/data/TFBS_footprinter_analysis/human/tfbs_results) TFBS_footprinter run human results.
* (/data_raw/Allen_Brain_Atlas_cortical) Single-cell RNA-Seq data from human cortex -- data needs to be manually downloaded as described above and placed here.

__OUTPUTS:__
* (/output/figures/cortical_single_cell) Force directed graph plots of cell type clusters and matrix plot figures and tables.





---
### __999.gtrd_peaking.3prime_5prime.py__
**DESCRIPTION (IMPORTANT: Run can require up to 70+GB of RAM to run to completion):**
Map high-scoring peaks from 21,988 human ChIP-Seq experiments to all human promoter regions to identify occupancy regions worth targeting for TFBS_footprinter analyses, i.e. what distance upstream and downstreatm of transcription start sites. In order to run this analysis the meta-analysis of ChIP-Seq experiment data compiled by the GTRD project [https://gtrd.biouml.org](https://gtrd.biouml.org) needs to be manually downloaded and placed in the '/data_raw/GTRD' directory. Only the most expressed transcripts from each human gene are used to calculate the occupancy [https://gtexportal.org/home/datasets](https://gtexportal.org/home/datasets), and this data should be placed in the '/data_raw/GTEx' directory.
The relevant files are:
 * [GTRD ChIP-Seq MACS peaks human (\~8 GB)](http://gtrd.biouml.org:8888/downloads/current/intervals/chip-seq/Homo_sapiens_ChIP-seq_peaks_MACS2.zip)
 * [GTEx transcript counts (\~3.5 GB)](https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz)
 
__INPUTS:__
* (/data_raw/Ensembl) Ensembl GTF data.
* (/data_raw/GTRD/Homo_sapiens_ChIP-seq_peaks_MACS2) ChIP-Seq data from GTRD project -- data needs to be manually downloaded as described above and placed here.
* (/data_raw/GTEx) Trascript level expression counts data from GTEx project -- data needs to be manually downloaded as described above and placed here.

__OUTPUTS:__
* (/output/figures) Histogram of promoter occupancy counts.


## TFBS_footprinter manuscript analyses

## 1 Background
This repository provides Python scripts used to perform analyses and generate images for the pre-print article found at [biorxiv](https://www.biorxiv.org/content/10.1101/2020.09.04.282749v2)
Full documentation available at: [ReadTheDocs](https://tfbs-footprinting.readthedocs.io/en/latest/index.html)

## 2 Instructions
We will give instructions to set up a miniconda environment to contain the dependencies needed and to run the manuscript's analysis scripts.

### 2.1 Install miniconda
Miniconda is the barebones version of the larger Conda package. We will use this so that we choose only the dependencies that are needed and therefore reduce the installation size and time. The Miniconda installation instructions are here: [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

### 2.2 Create environment
```
    $ conda create -n tfbs_footprinter_analyses python=3.7
```

### 2.3 Install analyses dependencies
```
    $ conda activate tfbs_footprinter_analyses
    $ conda install -c anaconda seaborn
    $ conda install -c bioconda pyliftover
```

## 3 Script files
__000.get_data.py__ - Will download Ensembl GTF data and Neanderthal/Human/Chimp variant data files, each ~<50MB. Variant data will be unzipped, and its contents are ~250MB.


__001.neanderthal_variants.py__ - Identify modern human vs. Neanderthal single-nucleotide polymorphisms (SNPs) which occur within 2,500 bp (upstream/downstream) of a human protein-coding transcript transcription start site.


__002.tfbs_footprinter_run.human.chunks.py__ - Run a transcription factor binding site (TFBS) prediction using tfbs_footprinter on the human version of the human vs. Neanderthal SNPs, which have been identified within 2,500 bp of a human protein-coding transcript transcription start site (TSS).


__003.human2neandertal_chimp.py__ - Downloaded files from the human analysis are copied to a neanderthal dir, and edited to match the variant cataloged in the CatalogOfChanges data. TFBS_footprinter will in the next step be run on these created Neanderthal sequence versions.


__004.tfbs_footprinter_run.neanderthal.chunks.py__ - Run a TFBS prediction using tfbs_footprinter on the Neanderthal version of all human vs. Neanderthal SNPs, which have been identified within 2,500 bp of a human protein-coding transcript transcription start site.


__005.tf_frame_score_changes.py__ - Match and analyze differences in predicted TFBSs at modern human vs. Neanderthal SNPs. Produce a list of those TFs which have a statistically significant difference in binding between species; differentially binding TFs (DBTFs).


__006.clustermap.py__ - Using the list of DBTFs, extract expression data from the FANTOM dataset and perform cluster analysis and produce clustermap figure. Using data from The Protein Atlas, identify which DBTFs have tissue-specific expression and produce matrix figure.
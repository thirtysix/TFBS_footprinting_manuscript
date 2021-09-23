

## TFBS_footprinter manuscript analyses

## 1 Background
This repository provides Python scripts used to perform analyses and generate images for the pre-print article found at [biorxiv](https://www.biorxiv.org/content/10.1101/2020.09.04.282749v2)
Full documentation available at: [ReadTheDocs](https://tfbs-footprinting.readthedocs.io/en/latest/index.html)

## 2 Instructions
We will give instructions to set up a miniconda environment to contain the dependencies needed and to run the manuscript's analysis scripts.

### 2.1 Install miniconda
Miniconda is the barebones version of the larger Conda package.  We will use this so that we choose only the dependencies that are needed and therefore reduce the installation size and time.  The Miniconda installation instructions are here: [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

### 2.2 Create environment
    `$ conda create -n tfbs_footprinter_analyses python=3.7

### 2.3 Install analyses dependencies
    `$ conda activate tfbs_footprinter_analyses
    `$ conda install -c anaconda seaborn
    `$ conda install -c bioconda pyliftover




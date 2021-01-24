# HRDecipher

**Table of Contents**

1. [Introduction] (#introduction)
2. [Installation] (#installation)
3. [Usage] (#usage)



## Introduction

HRDecipher is a simple Python script that computes HRD genomic scars and plot each scar on chromosomes. 

## Installation 

HRDecipher can be installed via git

```git clone https://github.com/ZKai0801/HRDecipher.git```

Python package `pandas` is required:

```pip install pandas```

R package `argparse` and `karyoploteR` is also required to visualise genomic scars

```R
install.packages('argparse');

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("karyoploteR")
```



## Usage

```bash
[kai@admin HRDecipher]$ python3 HRDecipher.py -h
usage: HRDecipher.py [-h] [-o OUTPUT] input

positional arguments:
  input                 sampleID.pre_hrd.tsv, must contain following columns:
                        Chromosome, Start_position, End_position, total_cn,
                        A_cn, B_cn, ploidy

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output file prefix
```

`Sequenza` is by far the most popular CNV caller that used in HRD calculation. But of course, any CNV callers that incorporates celluarity and ploidy calculation will do. 

Here is an example pipeline:  [Sequenza_HRD.sh](https://github.com/ZKai0801/HRDecipher/blob/master/Sequenza_HRD.sh)

Beware that although NGS data from many targeted-sequencing panels can be used to calculate HRD scars, ideally only WES/WGS data or panels that specifically designed that uniformly enriched on heterozygous sites (similar with SNP-array technology) are suitable for this aim.
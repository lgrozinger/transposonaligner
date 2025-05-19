# Overview

TnAtlas is a Python package for identifying and annotating transposon integration events into genomes.

Given a set of sequencing reads, transposon sequences, and genomes, the TnAtlas package can:

* Looks for reads which contain genomic DNA preceded by transposon DNA.
* Annotate the reads with corresponding features from the genome.
* Produces a summary for a set of reads in excel format.

The package ships with 2 utilities, `tnfind` and `tnmeta`, which can be used to run analysis from the command line. 
'tnfind' uses Blastn to align sequencing reads to a given transposon plasmid sequence to identify the transposon end and subsequently aligns the reads to a given genome. The annotations from the genome file are included in the output file results.xlsx. 
'tnmeta' adds metadata to your results.xlsx file. A usual metadata added is the plate layout to identify each sequencing read. 

# Installing

## Dependencies

* Python >= 3.8 
* blastn >= 2.12

### Optionally
Some parts of the pipeline also require

* fastqc (for sequencing quality control reports)
* sickle (for trimming based on sequencing quality)

## From source code

1. Get the code:
   `git clone https://github.com/biocomputationlab/transposonaligner`
3. Install using pip:
   
   `python3 -m pip install ./transposonaligner`

## From PyPI (using pip)

`python3 -m pip install tnatlas`

## Using Docker (recommended for Windows users to use trimming function with sickle)

`python3 -m pip install tnatlas`

# Usage

`tnfind sequencing_data path_to_results_folder -transposon transposon_file.gb -genome genome_file.gb -trim -sam` 

i.e (in data folder) `tnfind . ./results/results.xlsx -transposon transposons.gb -genome pputidakt2240.gb -trim -sam` 

`tnmeta -o path_to_results_folder/results.xlsx 'PLATE-' '-WELL-premix' output_file_name.xlsx`

i.e (in data folder) `tnmeta -o ./results.xlsx 'PLATE' 'WELL-premix' ./results(results.meta`

# Contributing

Contributions of all kinds are welcomed, including:
* Code for new features and bug fixes
* Test code and examples
* Bug reports and suggestions for new features (e.g. opening a github issue)

If you plan to change the source code, please open an issue for discussing the proposed changes and fork the repository.

# Citing

If you use this package as part of a publication, please cite: 

# Acknowledgements
This work was funded by grants BioSinT-CM (Y2020/TCS-6555) and CONTEXT
(Atracción de Talento Program; 2019-T1/BIO-14053) projects of the Comunidad de
Madrid, MULTI-SYSBIO (PID2020-117205GA-I00) and Severo Ochoa Program for
Centres of Excellence in R&D (CEX2020-000999-S) funded by
MCIN/AEI/10.13039/501100011033, and the ECCO (ERC-2021-COG-101044360)
contract of the EU.

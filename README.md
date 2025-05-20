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

`docker pull biocomputationlab/tnatlas:latest`

# Usage
## tnfind

`tnfind sequencing_data path_to_results_folder -transposon transposon_file.gb -genome genome_file.gb -trim -sam` 

**usage**: 

tnfind [-h] [-v] [--sequencing-type {sanger,solexa,illumina}] [--input-type {fasta,fastq}] [--input-ext INPUT_EXT] [-o OUTPUT_FILE] [-qc [PATH_TO_FASTQC]] [-trim [PATH_TO_SICKLE]]
              [--trim-quality TRIM_QUALITY] [--trim-length TRIM_LENGTH] [--trim-save] [-blastn PATH_TO_BLASTN] [-sam] -transposon TRANSPOSON [TRANSPOSON ...] [--transposon-type TRANSPOSON_TYPE]
              [--transposon-save] [--transposon-word-size TRANSPOSON_WORD_SIZE] [--transposon-evalue TRANSPOSON_EVALUE] -genome GENOME [GENOME ...] [--genome-type GENOME_TYPE] [--genome-save]
              [--genome-word-size GENOME_WORD_SIZE] [--genome-evalue GENOME_EVALUE] [--genome-prefix GENOME_PREFIX] [-config CONFIG]
              input_dir output_dir

Use a set of sequencing reads to identify the positions of likely transposon integration events in genomes. The integration events are identified using target sequences from the files provided in the
`transposon` argument. Integration events are position in the genomes by alignment of the remainder of the sequencing read. Likely integration events annotated and saved in genbank format. A summary
table is also produced containing detail on each likely event.

**positional arguments**:

  input_dir             The directory containing the reads to be processed
  
  output_dir            The directory in which to save results

**options**:

  -h, --help            show this help message and exit
  
  -v
  --sequencing-type {sanger,solexa,illumina} The method from which the single end reads are obtained (default: sanger)
   
  --input-type {fasta,fastq} The file format in which the single end reads are saved (default: fastq)
                        
  --input-ext INPUT_EXT The file extension of the single end read files (default: ab1 if INPUT_TYPE is fastq, fasta otherwise)
                        
  -o OUTPUT_FILE        A filename for the output summary table (default: results.xlsx)
  
  -qc [PATH_TO_FASTQC]  If present, fastqc analysis of read quality is performed
  
  -trim [PATH_TO_SICKLE]  If present, sickle is used to trim the reads
                        
  --trim-quality TRIM_QUALITY The quality threshold passed to sickle for trimming (default: 20)
                        
  --trim-length TRIM_LENGTH The length threshold passed to sickle for trimming (default: 20)
                        
  --trim-save           If present, the trimmed reads will be saved in OUTPUT_DIR
  
  -blastn PATH_TO_BLASTN The blastn executable (default: blastn)
                        
  -sam                  If present, save the blastn output in SAM format in OUTPUT_DIR
  
  -transposon TRANSPOSON [TRANSPOSON ...] The file or files containing the sequences used to identify integrations
                        
  --transposon-type TRANSPOSON_TYPE The file format of the file(s) given in TRANSPOSON (default: genbank)
                        
  --transposon-save     If present, save reads annotated with transposon sequences in OUTPUT_DIR
  
  --transposon-word-size TRANSPOSON_WORD_SIZE The word size passed to blastn for the transposon alignment (default: 10)
                        
  --transposon-evalue TRANSPOSON_EVALUE The evalue passed to blastn for the transposon alignment (default: 0.01)
                        
  -genome GENOME [GENOME ...] The file or files containing genomes for alignment
                        
  --genome-type GENOME_TYPE The file format of the file(s) given in GENOME (default: genbank)
                        
  --genome-save         If present, save reads annotated with all genome alignments in OUTPUT_DIR
  
  --genome-word-size GENOME_WORD_SIZE The word size passed to blastn for the genome alignment (default: 10)
                        
  --genome-evalue GENOME_EVALUE The evalue passed to blastn for the genome alignment (default: 0.01)
                        
  --genome-prefix GENOME_PREFIX  The number of base pairs, after the integration site, of genomic DNA to include in the output summary table
                        
  -config CONFIG        A configuration file in TOML format. Command line arguments are overrided with value in the configuration file, if present
  

i.e (in data folder) `tnfind . ./results/results.xlsx -transposon transposons.gb -genome pputidakt2240.gb -trim -sam` 

## tnmeta


`tnmeta -o path_to_results_folder/results.xlsx 'PLATE-' '-WELL-premix' output_file_name.xlsx`

**usage**: 

tnmeta [-h] [-o OUTPUT_FILE] plate_regex well_regex result_file metadata [metadata ...]

Attach well metadata to a results table from a plate map

**positional arguments**:

  plate_regex     A regular expression which matches plates in record names
  
  well_regex      A regular expression which matches wells in record names
  
  result_file     The spreadsheet to annotate with metadata
  
  metadata        The spreadsheets containing the metadata map of a plate

**options**:

  -h, --help      show this help message and exit
  
  -o OUTPUT_FILE  If given, output will be saved to OUTPUT_FILE, if not, RESULT_FILE will be overwritten

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

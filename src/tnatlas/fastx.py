from Bio import SeqIO

import pathlib


def directory_to_fastq_records(directory, sequencing_type, extension):
    "Takes sequencing quality files in a directory and generates SeqRecords."
    directory = pathlib.Path(directory).resolve(strict=True)

    quality_format = extension
    if extension == "fastq" and sequencing_type == "solexa":
        quality_format = "fastq-solexa"
    if extension == "fastq" and sequencing_type == "illumina":
        quality_format = "fastq-illumina"
    if extension == "ab1":
        quality_format = "abi"

    for quality_file in directory.glob(f"*.{extension}"):
        quality_file = quality_file.resolve(strict=True)
        for record in SeqIO.parse(quality_file, quality_format):
            yield record

def directory_to_fasta_records(directory, extension):
    "Takes FASTA files in a directory and generates SeqRecords."
    directory = pathlib.Path(directory).resolve(strict=True)
    
    for fasta_file in directory.glob(f"*.{extension}"):
        fasta_file = fasta_file.resolve(strict=True)
        for record in SeqIO.parse(fasta_file, "fasta"):
            record.annotations["molecule_type"] = "DNA"
            yield record

def fastq_to_fasta(fastq_path, fasta_path):
    "Takes a FASTQ filepath, converts sequences to FASTA and returns a FASTA filepath."
    fastq_path = pathlib.Path(fastq_path).resolve(strict=True)
    fasta_path = pathlib.Path(fasta_path).resolve()
    SeqIO.write(SeqIO.parse(fastq_path, "fastq"), fasta_path, "fasta")

def fasta_to_fastq(fasta_path, fastq_path):
    "Takes a FASTA filepath, converts sequences to FASTQ and returns a FASTQ filepath."
    fastq_path = pathlib.Path(fastq_path).resolve(strict=True)
    fasta_path = pathlib.Path(fasta_path).resolve()
    SeqIO.write(SeqIO.parse(fasta_path, "fasta"), fastq_path, "fastq")

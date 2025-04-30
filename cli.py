#!/usr/local/bin/python3
from pathlib import Path
from Bio import SeqIO

import argparse
import pandas
import fastx
import qc
import blastn
import insertion

from config import CONFIG


DATADIR = Path(CONFIG["input_dir"]).resolve(strict=True)

records = []
if CONFIG["input_type"] == "fasta":
    records = list(fastx.directory_to_fasta_records(DATADIR, CONFIG["input_ext"]))
elif CONFIG["input_type"] == "fastq":
    records = list(fastx.directory_to_fastq_records(DATADIR, CONFIG["sequencing_type"], CONFIG["input_ext"]))

reads = {record.id: insertion.Read(record) for record in records}

print(f"{len(reads)} read{'s' if len(reads) > 1 else ''} from {DATADIR} loaded")

OUTDIR = Path(CONFIG["output_dir"]).resolve(strict=True)

if CONFIG["qc"] is not None:
    qc.quality_control(reads.values(), OUTDIR, CONFIG["qc"])

if CONFIG["trim"] is not None:
    trimmed_reads = qc.trimming(
        reads.values(),
        CONFIG["sequencing_type"],
        CONFIG["trim_quality"],
        CONFIG["trim_length"],
        CONFIG["trim"]
    )
    ntrimmed = len(reads) - len(trimmed_reads)
    if ntrimmed > 0:
        print(f"Trimmed {ntrimmed} read{'s' if ntrimmed > 1 else ''} for quality")
    else:
        print(f"Did not trim any reads for quality")
    reads = {read.id: insertion.Read(read) for read in trimmed_reads}

    if CONFIG["trim_save"]:
        SeqIO.write(reads.values(), OUTDIR / "records.trimmed.fastq", "fastq")


transposon_path = Path(CONFIG["transposon"])
if not transposon_path.is_absolute():
    transposon_path = DATADIR / transposon_path
transposon_path = transposon_path.resolve(strict=True)

transposons = SeqIO.parse(transposon_path, CONFIG["transposon_type"])
transposons = {transposon.id: transposon for transposon in transposons}

print(f"{len(transposons)} transposon{'s' if len(transposons) > 1 else ''} loaded")

genome_path = Path(CONFIG["genome"])
if not genome_path.is_absolute():
    genome_path = DATADIR / genome_path
genome_path = genome_path.resolve(strict=True)

genomes = SeqIO.parse(genome_path, CONFIG["genome_type"])
genomes = {genome.id: genome for genome in genomes}

print(f"{len(genomes)} genome{'s' if len(genomes) > 1 else ''} loaded")

insertion.insertion_search(
    reads.values(),
    transposons.values(),
    evalue=CONFIG["transposon_evalue"],
    word_size=CONFIG["transposon_word_size"],
)

print(f"{sum(read.has_insertion for read in reads.values())} reads have transposon present")

if CONFIG["transposon_save"]:
    for read in reads.values():
        SeqIO.write(read, OUTDIR / f"{read.id}.transposon.aligned.gb", "genbank")

for read in reads.values():
    read.choose_insertion()

insertion.genome_search(
    reads.values(),
    genomes.values(),
    evalue=CONFIG["genome_evalue"],
    word_size=CONFIG["genome_word_size"],
)

print(f"{sum(read.has_genome for read in reads.values())} reads have genomic DNA present")
print(f"{sum(read.has_genome and read.has_insertion for read in reads.values())} reads have both transposon and genomic DNA present")


if CONFIG["genome_save"]:
    for read in reads.values():
        SeqIO.write(read, OUTDIR / f"{read.id}.genome.aligned.gb", "genbank")

for read in reads.values():
    read.choose_genome()

for read in reads.values():
    SeqIO.write(read, OUTDIR / f"{read.id}.processed.gb", "genbank")

data = pandas.concat([read.dataframe for read in reads.values()])
data.to_csv(OUTDIR / CONFIG["o"])
print(f"Saving result table to {OUTDIR / CONFIG['o']}")



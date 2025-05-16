#!/usr/local/bin/python3
from pathlib import Path
from Bio import SeqIO


import argparse
import pandas
from tnatlas import utils
from tnatlas import blastn
from tnatlas import qc
from tnatlas import insertion
from tnatlas import config


def find_transposons(reads, config):
    transposons = utils.load_transposons(config)
    insertion.insertion_search(
        reads.values(),
        transposons.values(),
        command=config["blastn"],
        evalue=config["transposon_evalue"],
        word_size=config["transposon_word_size"],
    )

    if config["transposon_save"]:
        outdir = Path(config["output_dir"]).resolve(strict=True)
        for read in reads.values():
            SeqIO.write(read, outdir / f"{read.id}.transposon.aligned.gb", "genbank")

    for read in reads.values():
        read.choose_insertion()

def find_genomes(reads, config):
    genomes = utils.load_genomes(config)

    insertion.genome_search(
        reads.values(),
        genomes.values(),
        config["genome_prefix"],
        command=config["blastn"],
        evalue=config["genome_evalue"],
        word_size=config["genome_word_size"],
    )

    if config["genome_save"]:
        outdir = Path(config["output_dir"]).resolve(strict=True)
        for read in reads.values():
            SeqIO.write(read, config["output_dir"] / f"{read.id}.genome.aligned.gb", "genbank")

    if config["sam"]:
        outdir = Path(config["output_dir"]).resolve(strict=True)
        with open(outdir / "samgenomes.sam", "w") as f:
            print(blastn.sam(reads.values(), genomes.values()), file=f)
        
    for read in reads.values():
        read.choose_genome()


def main():
    try: 
        CONFIG = config.get_configuration(config.get_parser())
        DATADIR = Path(CONFIG["input_dir"]).resolve(strict=True)
        OUTDIR = Path(CONFIG["output_dir"]).resolve(strict=True)
    except FileNotFoundError as err:
        print(f"Could not find the file: {err.filename}")
        print("Please check the filenames given in arguments and whether they exist.")
        return None

    reads = utils.load_reads(CONFIG)

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
        reads = {read.id: insertion.Read(read) for read in trimmed_reads}

    if CONFIG["trim_save"]:
        SeqIO.write(reads.values(), OUTDIR / "records.trimmed.fastq", "fastq")

    find_transposons(reads, CONFIG)
    find_genomes(reads, CONFIG)


    for read in reads.values():
        SeqIO.write(read, OUTDIR / f"{read.id}.processed.gb", "genbank")

    data = pandas.concat([read.dataframe for read in reads.values()])
    data.to_excel(OUTDIR / CONFIG["o"])
    print(f"Saving summary table to {OUTDIR / CONFIG['o']}")

    return None

if __name__ == "__main__":
    main()

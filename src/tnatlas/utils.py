from pathlib import Path
from Bio import SeqIO

from tnatlas.insertion import Read
from tnatlas.fastx import directory_to_fasta_records
from tnatlas.fastx import directory_to_fastq_records
from tnatlas.qc import trimming


def load_reads(config):
    indir = config["input_dir"]
    ext = config["input_ext"]

    records = []
    if config["input_type"] == "fasta":
        records = directory_to_fasta_records(indir, ext)
    elif config["input_type"] == "fastq":
        records = directory_to_fastq_records(indir, config["sequencing_type"], ext)

    reads = {record.id: Read(record) for record in records} 
    print(f"{len(reads)} read{'' if len(reads) < 2 else 's'} from {indir} loaded")
    return reads

def load_transposons(config):
    ts = {}
    for path in config["transposon"]:
        path = Path(path).resolve(strict=True)
        print(f"Loading transposons from {path}")
        records = SeqIO.parse(path, config["transposon_type"])
        ts |= {record.id: record for record in records}

    print(f"{len(ts)} transposon{'' if len(ts) < 2 else 's'} loaded")
    return ts

def load_genomes(config):
    gs = {}
    for path in config["genome"]:
        path = Path(path).resolve(strict=True)
        print(f"Loading genomes from {path}")
        records = SeqIO.parse(path, config["genome_type"])
        gs |= {record.id: record for record in records}

    print(f"{len(gs)} genome{'' if len(gs) < 2 else 's'} loaded")
    return gs

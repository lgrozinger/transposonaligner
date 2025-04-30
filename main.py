from Bio import SeqIO

from pathlib import Path

import pandas
import blastn
import fastx
import insertion

DATADIR = Path("/home/lewis/transposonaligner/data")
DATADIR = DATADIR.resolve(strict=True)

seqrecords = fastx.directory_to_fasta_records(DATADIR, "fasta")
reads = {record.id: insertion.Read(record) for record in seqrecords}
vectors = {v.id: v for v in SeqIO.parse(DATADIR / "vectors.gb", "genbank")}
vectors = {
    "pBAMD1-2-MarC9": vectors["pBAMD1-2-MarC9"][1182:1210],
    "pBAMD1-2-pBAD-YFP_(4A)": vectors["pBAMD1-2-pBAD-YFP_(4A)"][2533:2552],
}

SeqIO.write(vectors.values(), DATADIR / "transposons.gb", "genbank")

insertion.insertion_search(reads.values(), vectors.values(), evalue=0.01)

genome = {g.id: g for g in SeqIO.parse(DATADIR / "pputidakt2440.gb", "genbank")}

insertion.genome_search(reads.values(), genome.values(), evalue=0.01, word_size=12)

for read in reads.values():
   read.choose_insertion()
   read.choose_genome()

for read in reads.values():
    SeqIO.write(read, DATADIR / f"results/{read.id}.processed.gb", "genbank")

data = pandas.concat([read.dataframe for read in reads.values()])
data.to_csv("/tmp/data.csv")

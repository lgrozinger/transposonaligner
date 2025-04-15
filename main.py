from Bio import SeqIO

from pathlib import Path

import sequences
import blastn
import fastx
import insertion

DATADIR = Path("/home/lewis/transposonaligner/data")
DATADIR = DATADIR.resolve(strict=True)

seqrecords = fastx.directory_to_fasta_records(DATADIR, "fasta")
insertions = {record.id: insertion.Insertion(record) for record in seqrecords}
vectors = {v.id: v for v in SeqIO.parse(DATADIR / "vectors.gb", "genbank")}
genome = {g.id: g for g in SeqIO.parse(DATADIR / "pputidakt2440.gb", "genbank")}

vector_sequences = {
    "pBAMD1-2-MarC9": sequences.region_between(vectors["pBAMD1-2-MarC9"], 1065, 1210),
    "pBAMD1-2-pBAD-YFP_(4A)": sequences.region_between(vectors["pBAMD1-2-pBAD-YFP_(4A)"], 2416, 2552),
}

insertion.transposon_search(insertions.values(), vector_sequences.values())

suffixes = [i.suffix for i in insertions.values()]
blastn.align(suffixes, genome.values(), word_size=8)

for i in suffixes:
    insertions[i.id].merge_features_with(i)
    
for insertion in insertions.values():
    SeqIO.write(insertion, DATADIR / f"results/{insertion.id}.processed.gb", "genbank")

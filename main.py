from Bio import SeqIO

from pathlib import Path

import sequences
import blastn
import fastx
import insertion

DATADIR = Path("/home/lewis/transposonaligner/data")
DATADIR = DATADIR.resolve(strict=True)

seqrecords = fastx.directory_to_fasta_records(DATADIR, "fasta")
reads = {record.id: insertion.Read(record) for record in seqrecords}
vectors = {v.id: v for v in SeqIO.parse(DATADIR / "vectors.gb", "genbank")}
vectors = {
    "pBAMD1-2-MarC9": vectors["pBAMD1-2-MarC9"][1065:1210],
    "pBAMD1-2-pBAD-YFP_(4A)": vectors["pBAMD1-2-pBAD-YFP_(4A)"][2416:2552],
}

insertion.transposon_search(reads.values(), vectors.values())

genome = {g.id: g for g in SeqIO.parse(DATADIR / "pputidakt2440.gb", "genbank")}

# suffixes = [i.suffix for i in reads.values()]
# blastn.align(suffixes, genome.values(), word_size=8)

# for i in suffixes:
#     reads[i.id].merge_features_with(i)
    
for read in reads.values():
    SeqIO.write(read, DATADIR / f"results/{read.id}.processed.gb", "genbank")

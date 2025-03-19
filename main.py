from Bio import SeqIO

from pathlib import Path

import itertools
import sequences
import blastn
import fastx
from insertion import Insertion
from insertion import transposon_search

DATADIR = Path("/home/lewis/transposonaligner/data")
DATADIR = DATADIR.resolve(strict=True)

seqrecords = fastx.directory_to_fasta_records(DATADIR, "fasta")
insertions = {record.id: Insertion(record) for record in seqrecords}
vectors = {v.id: v for v in SeqIO.parse(DATADIR / "vectors.gb", "genbank")}
genome = {g.id: g for g in SeqIO.parse(DATADIR / "pputidakt2440.gb", "genbank")}

vector_sequences = {
    "pBAMD1-2-MarC9": sequences.region_between(vectors["pBAMD1-2-MarC9"], 1065, 1210),
    "pBAMD1-2-pBAD-YFP_(4A)": sequences.region_between(vectors["pBAMD1-2-pBAD-YFP_(4A)"], 2416, 2552),
}

transposon_search(insertions.values(), vector_sequences.values())
number_transposon = len(list(filter(lambda x: x.has_transposon, insertions.values())))
blastn.align(insertions.values(), genome.values())

for insertion in insertions.values():
    SeqIO.write(insertion, DATADIR / f"results/{insertion.id}.processed.gb", "genbank")

# alignment_annotation(records.values(), vector_sequences.values())

# for record in records.values():
#     keep_only_best_hsp(record)

# alignment_annotation(map(sequences.suffix, records.values()), genome.values())

# for record in records.values():
#     SeqIO.write(record, DATADIR / f"results/{record.id}.processed.gb", "genbank")

# for record in vectors.values():
#     SeqIO.write(record, DATADIR / f"results/{record.id}.processed.gb", "genbank")



# def annotate_record_from_target(record, feature, target):
#     logger.info(f"Appending feature from {target.name} to {record.id}")
#     record.features.append(feature)

#     query_location = feature.location
#     target_location = hsp_to_location(feature.hsp, target=True)
#     target_features = truncated_features(target.features, target_location)
#     for target_feature in target_features:
#         target_feature = translate_feature(target_feature, target_location, query_location)
#         record.features.append(target_feature)
        
#     return record


# def keep_best_hsp_feature(record):
#     hsps = [feature.hsp for feature in record.features if feature.type == "hsp"]
#     if not hsps:
#         return record
    
#     best_score = hsp_score(best_hsp(hsps))
#     def f(feature):
#         return (not feature.type == "hsp") or (hsp_score(feature.hsp) >= best_score)
#     record.features = list(filter(f, record.features))
#     return record

# def align_and_annotate(records, queries, targets):
#     for alignment in blastn(queries, targets):
#         record = records[alignment.query.id]
#         for feature in alignment_to_features(alignment, target=False):
#             target = next(filter(lambda x: x.name == feature.hsp.target.name, targets))
#             annotate_record_from_target(record, feature, target)
#     return None


# for record in records.values():
#     keep_best_hsp_feature(record)

# align_and_annotate(records, map(suffix, records.values()), list(genome()))

# def find_unknown_sources(record):
#     nfeatures = len(record.features)
#     for location in unannotated_locations(record):
#         unknown_region = region_between(record, location.start, location.end)
#         align_and_annotate(records, [unknown_region], list(vectors()) + list(genome()))
#         if len(record.features) > nfeatures:
#             return True
#     return False

# identified_regions = []
# for record_id in records:
#     record = records[record_id]
#     logger.info(f"Searching for missing regions in {record_id}")
#     while find_unknown_sources(record):
#         identified_regions.append(record.id)
#         logger.info(f"Found missing region in {record.id}")
        
                                                       


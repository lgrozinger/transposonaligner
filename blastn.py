"Adds SeqFeatures to SeqRecords where there are BLASTN HSPs."

import subprocess
import tempfile
import pathlib
import psutil
import io

from Bio.SeqFeature import SimpleLocation
from Bio.SeqFeature import SeqFeature
from Bio import Blast
from Bio import SeqIO

def nthreads():
    "Returns half of the number of (virtual) threads available on system."
    n = psutil.cpu_count()
    if n is None:
        n = 1
    return max(1, n // 2)

def blastn(query_records, target_records, **kwargs):
    "Does a BLASTN query for query_records against target_records."
    with tempfile.TemporaryDirectory() as directory:
        query_fasta = pathlib.Path(directory).resolve() / "query.fasta"
        SeqIO.write(query_records, query_fasta, "fasta")

        target_fasta = pathlib.Path(directory).resolve() / "target.fasta"
        SeqIO.write(target_records, target_fasta, "fasta")

        cmd = ["blastn", "-query", str(query_fasta)]
        cmd += ["-subject", str(target_fasta)]
        cmd += ["-parse_deflines"]
        cmd += ["-num_threads", str(nthreads())]
        cmd += ["-mt_mode", str(1)]
        cmd += ["-outfmt", "5"]
        for argument in kwargs:
            cmd += [f"-{argument}", str(kwargs[argument])]

        print(" ".join(cmd))

        result = subprocess.run(cmd, capture_output=True, check=True)
        return Blast.parse(io.BytesIO(result.stdout))

def hsp_location(hsp, target=True):
    "Converts a HSP from a blast result to a Bio.SeqFeature.SimpleLocation."
    sequence = 0 if target else 1
    x = int(hsp.coordinates[sequence, 0])
    y = int(hsp.coordinates[sequence, -1])
    strand = 1 if x <= y else -1
    return SimpleLocation(min(x, y), max(x, y), strand=strand)

def hsp_feature(hsp, target=True):
    "Converts a HSP from a blast result to a Bio.SeqFeature.SeqFeature."
    location = hsp_location(hsp, target=target)
    qualifiers = {"label": hsp.query.id if target else hsp.target.name}
    feature = SeqFeature(location, type="misc_feature", qualifiers=qualifiers)
    feature.hsp = hsp
    return feature

def get_query(hit, records):
    return next(filter(lambda x: x.id == hit[0].query.id, records))

def get_target(hit, records):
    return next(filter(lambda x: x.name == hit.target.name, records))

def align(query_records, target_records, **kwargs):
    for alignment in blastn(query_records, target_records, **kwargs):
        for hit in alignment:
            query_record = get_query(hit, query_records)
            target_record = get_target(hit, target_records)
            for hsp in hit:
                query_feature = hsp_feature(hsp, target=False)
                query_record.features.append(query_feature)
                hsp.query = query_record
                
                target_feature = hsp_feature(hsp, target=True)
                target_record.features.append(target_feature)
                hsp.target = target_record

    return None

def hsp_score(hsp):
    "Returns the bit score of hsp, with a tie breaker of the hsp length."
    return (hsp.annotations["bit score"], hsp.length)

def best_hsp(hsps):
    "Returns the HSP with the best bit score."
    return sorted(hsps, key=hsp_score)[-1] if hsp else None    

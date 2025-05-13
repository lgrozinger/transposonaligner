"Does BLASTN between SeqRecords and works with BLASTN HSPs."

from Bio.SeqFeature import SimpleLocation
from Bio.SeqFeature import SeqFeature
from Bio import Blast
from Bio import SeqIO

from io import BytesIO
from subprocess import run
from tempfile import TemporaryDirectory
from pathlib import Path
from psutil import cpu_count

def nthreads():
    "Returns half of the number of (virtual) threads available on system."
    n = cpu_count() if cpu_count() is not None else 1
    return max(1, n // 2)

def sam(query_records, target_records, command="blastn", **kwargs):
    "Does a BLASTN query and return SAM output."
    with TemporaryDirectory() as directory:
        query_fasta = Path(directory).resolve() / "query.fasta"
        SeqIO.write(query_records, query_fasta, "fasta")

        target_fasta = Path(directory).resolve() / "target.fasta"
        SeqIO.write(target_records, target_fasta, "fasta")

        outputfile = Path(directory).resolve() / "output.sam"

        cmd = [command, "-query", str(query_fasta)]
        cmd += ["-subject", str(target_fasta)]
        cmd += ["-parse_deflines"]
        cmd += ["-out", str(outputfile)]
        cmd += ["-outfmt", str(17)]
        for argument in kwargs:
            cmd += [f"-{argument}", str(kwargs[argument])]

        print(" ".join(cmd))
        result = run(cmd, capture_output=True, check=True)

        with open(outputfile, "r") as f:
            return f.read()

def blastn(query_records, target_records, command="blastn", **kwargs):
    "Does a BLASTN query for query_records against target_records."
    with TemporaryDirectory() as directory:
        query_fasta = Path(directory).resolve() / "query.fasta"
        SeqIO.write(query_records, query_fasta, "fasta")

        target_fasta = Path(directory).resolve() / "target.fasta"
        SeqIO.write(target_records, target_fasta, "fasta")

        cmd = [command, "-query", str(query_fasta)]
        cmd += ["-subject", str(target_fasta)]
        cmd += ["-parse_deflines"]
        cmd += ["-outfmt", str(5)]
        for argument in kwargs:
            cmd += [f"-{argument}", str(kwargs[argument])]

        print(" ".join(cmd))

        result = run(cmd, capture_output=True, check=True)
        return Blast.parse(BytesIO(result.stdout))

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

def hsp_score(hsp):
    "Returns the bit score of hsp, with a tie breaker of the hsp length."
    return (hsp.annotations["bit score"], hsp.length)

def best_hsp(hsps):
    "Returns the HSP with the best bit score."
    return sorted(hsps, key=hsp_score)[-1] if hsp else None    

from Bio import SeqIO

from io import StringIO
from subprocess import run
from pathlib import Path
from tempfile import TemporaryDirectory

def quality_control(records, outdir, program="fastqc"):
    "Takes SeqRecords and runs them through fastqc."
    outdir = Path(outdir).resolve(strict=True)
    with TemporaryDirectory() as directory:
        directory = Path(directory).resolve()
        fastq_path = directory / "records.fastq"
        SeqIO.write(records, fastq_path, "fastq")

        cmd = [program]
        cmd += [str(fastq_path)]
        cmd += ["--extract", "--delete", "--quiet"]
        cmd += ["--outdir", str(outdir)]
        
        return run(cmd, capture_output=True, check=True)

def trimming(records, sequencing_type="sanger", quality=20, length=20, program="sickle"):
    "Trims FASTQ SeqRecords with sickle, and returns a list of the trimmed SeqRecords."
    with TemporaryDirectory() as directory:
        directory = Path(directory).resolve()
        fastq_path = directory / "records.fastq"
        output_path = directory / "sickleout.fastq"
        SeqIO.write(records, fastq_path, "fastq")
        
        cmd = [program, "se"]
        cmd += ["--fastq-file", str(fastq_path)]
        cmd += ["--qual-type", str(sequencing_type)]
        cmd += ["--output-file", str(output_path)]
        cmd += ["-q", str(quality)]
        cmd += ["-l", str(length)]
        
        output = run(cmd, capture_output=True, check=True)
        print(output.stdout.decode("utf-8"))
        print(output.stderr.decode("utf-8"))

        with open(output_path, "r") as output_file:
            records = []
            for record in SeqIO.parse(output_file, "fastq"):
                record.annotations["molecule_type"] = "DNA"
                records.append(record)
                
        return records

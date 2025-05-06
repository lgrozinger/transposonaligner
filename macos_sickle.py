from transposonaligner import fastx

from Bio import SeqIO

import subprocess

import pathlib


def main():
    directory = pathlib.Path(__file__).parent / "data"
    records = fastx.directory_to_fastq_records(directory, "sanger", "ab1")
    fastq_path = pathlib.Path(directory) / "sicklein.fastq"
    SeqIO.write(records, fastq_path , "fastq")

    cmd = ["sickle", "se"]
    cmd += ["--fastq-file", str(fastq_path)]
    cmd += ["--output-file", pathlib.Path(directory) / "sickeout.fastq"]
    cmd += ["--qual-threshold", "20"]
    cmd += ["--length-threshold", "20"]

    output = run(cmd, capture_output=True)

    print(output.stdout.decode("utf-8"))
    print(output.stderr.decode("utf-8"))

if __name__ == "__main__":
    main()

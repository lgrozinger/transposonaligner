import os
import platform
import argparse
import tomllib
from pathlib import Path

ROOT = Path(__file__).parent.resolve(strict=True)
DEFAULTCONFIG = Path(ROOT / "default.toml").resolve(strict=True)

def load_config(path_or_string=DEFAULTCONFIG):
    path = Path(path_or_string).resolve(strict=True)
    with open(path, "rb") as f:
        return tomllib.load(f)

def get_config_property(config, *args):
    result = config
    for key in args:
        result = result.get(key, None)
        if result is None:
            return None
    return result

def add_config_property(to, key, config, *args):
    if get_config_property(config, *args) is not None:
        to[key] = get_config_property(config, *args)
    return None

def config_to_dict(config):
    d = {}
    add_config_property(d, "v", config, "output", "verbosity")
    add_config_property(d, "input_type", config, "input", "type")
    add_config_property(d, "input_ext", config, "input", "extension")
    add_config_property(d, "sequencing_type", config, "input", "sequencing_type")
    add_config_property(d, "o", config, "output", "outfile")
    add_config_property(d, "qc", config, "fastqc", "path")
    add_config_property(d, "trim", config, "sickle", "path")
    add_config_property(d, "trim_quality", config, "sickle", "quality")
    add_config_property(d, "trim_length", config, "sickle", "length")
    add_config_property(d, "trim_save", config, "sickle", "output")
    add_config_property(d, "blastn", config, "blastn", "path")
    add_config_property(d, "sam", config, "blastn", "sam")
    add_config_property(d, "transposon_type", config, "transposon", "type")
    add_config_property(d, "transposon_save", config, "transposon", "output")
    add_config_property(d, "transposon_word_size", config, "transposon", "word_size")
    add_config_property(d, "transposon_evalue", config, "transposon", "evalue")
    add_config_property(d, "genome_type", config, "genome", "type")
    add_config_property(d, "genome_save", config, "genome", "output")
    add_config_property(d, "genome_word_size", config, "genome", "word_size")
    add_config_property(d, "genome_evalue", config, "genome", "evalue")
    add_config_property(d, "genome_prefix", config, "genome", "prefix")
    add_config_property(d, "genome_window", config, "genome", "window")
    return d

def get_parser():
    parser = argparse.ArgumentParser(
        prog="tnfind",
        description = """Use a set of sequencing reads to identify the positions of likely transposon integration events in genomes.

        
The integration events are identified using target sequences from the files provided in the `transposon` argument.

        
Integration events are position in the genomes by alignment of the remainder of the sequencing read.

        
Likely integration events annotated and saved in genbank format. A summary table is also produced containing detail on each likely event.""",
    )

    parser.add_argument("-v", action="count", default=0)

    parser.add_argument(
        "input_dir",
        help="The directory containing the reads to be processed",
    )

    parser.add_argument(
        "--sequencing-type",
        default="sanger",
        choices=["sanger", "solexa", "illumina"],
        help="The method from which the single end reads are obtained (default: sanger)",
    )

    parser.add_argument(
        "--input-type",
        default="fastq",
        choices=["fasta", "fastq"],
        help="The file format in which the single end reads are saved (default: fastq)",
    )

    parser.add_argument(
        "--input-ext",
        default="ab1",
        help="The file extension of the single end read files (default: ab1 if INPUT_TYPE is fastq, fasta otherwise)",
    )

    parser.add_argument(
        "output_dir",
        help="The directory in which to save results",
    )
    
    parser.add_argument(
        "-o",
        default="results.xlsx",
        metavar="OUTPUT_FILE",
        help="A filename for the output summary table (default: results.xlsx)",
    )

    parser.add_argument(
        "-qc",
        nargs="?",
        const="fastqc.exe" if platform.system() == "Windows" else "fastqc",
        default=None,
        metavar="PATH_TO_FASTQC",
        help="If present, fastqc analysis of read quality is performed",
    )

    parser.add_argument(
        "-trim",
        nargs="?",
        const="sickle.exe" if platform.system() == "Windows" else "sickle",
        default=None,
        metavar="PATH_TO_SICKLE",
        help="If present, sickle is used to trim the reads",
        
    )

    parser.add_argument(
        "--trim-quality",
        default=20,
        type=int,
        help="The quality threshold passed to sickle for trimming (default: 20)",
    )
    
    parser.add_argument(
        "--trim-length",
        default=20,
        type=int,
        help="The length threshold passed to sickle for trimming (default: 20)"
    )

    parser.add_argument(
        "--trim-save",
        action="store_true",
        help="If present, the trimmed reads will be saved in OUTPUT_DIR", 
    )
    
    parser.add_argument(
        "-blastn",
        default="blastn.exe" if platform.system() == "Windows" else "blastn",
        metavar="PATH_TO_BLASTN",
        help="The blastn executable (default: blastn)",
    )

    parser.add_argument(
        "-sam",
        action="store_true",
        help="If present, save the blastn output in SAM format in OUTPUT_DIR",
    )


    parser.add_argument(
        "-transposon",
        required=True,
        nargs="+",
        help="The file or files containing the sequences used to identify integrations", 
    )
    
    parser.add_argument(
        "--transposon-type",
        default="genbank",
        help="The file format of the file(s) given in TRANSPOSON (default: genbank)",
    )
    
    parser.add_argument(
        "--transposon-save",
        action="store_true",
        help="If present, save reads annotated with transposon sequences in OUTPUT_DIR"
    )
    
    parser.add_argument(
        "--transposon-word-size",
        default=10,
        type=int,
        help="The word size passed to blastn for the transposon alignment (default: 10)",
    )
    
    parser.add_argument(
        "--transposon-evalue",
        default=0.01,
        type=float,
        help="The evalue passed to blastn for the transposon alignment (default: 0.01)",
    )
    
    parser.add_argument(
        "-genome",
        required=True,
        nargs="+",
        help="The file or files containing genomes for alignment",
    )
    
    parser.add_argument(
        "--genome-type",
        default="genbank",
        help="The file format of the file(s) given in GENOME (default: genbank)",
    )
    
    parser.add_argument(
        "--genome-save",
        action="store_true",
        help="If present, save reads annotated with all genome alignments in OUTPUT_DIR",
    )
    
    parser.add_argument(
        "--genome-word-size",
        default=10,
        type=int,
        help="The word size passed to blastn for the genome alignment (default: 10)",
    )
    
    parser.add_argument(
        "--genome-evalue",
        default=0.01,
        type=float,
        help="The evalue passed to blastn for the genome alignment (default: 0.01)",
    )
    
    parser.add_argument(
        "--genome-prefix",
        default=9,
        type=int,
        help="The number of base pairs, after the integration site, of genomic DNA to include in the output summary table", 
    )

    parser.add_argument(
        "-config",
        help="A configuration file in TOML format. Command line arguments are overrided with value in the configuration file, if present"
    )

    return parser

def absolute_path(path, directory):
    path = Path(path)
    if not path.is_absolute():
        directory = Path(directory).resolve(strict=True)
        path = directory / path
    return path
    
def make_paths_absolute(config):
    indir = absolute_path(config["input_dir"], os.getcwd())
    config["input_dir"] = indir

    outdir = absolute_path(config["output_dir"], os.getcwd())
    config["output_dir"] = outdir

    config["transposon"] = [absolute_path(t, os.getcwd()) for t in config["transposon"]]
    config["genome"] = [absolute_path(t, os.getcwd()) for t in config["genome"]]
    return config
    
def get_configuration(parser):
    cliargs = vars(parser.parse_args())
    if cliargs["input_type"] == "fasta":
        cliargs["input_ext"] = "fasta"
        
    configuration = {}
    
    if cliargs["config"] is not None:
        configuration = load_config(cliargs["config"])
        
    return make_paths_absolute(cliargs | config_to_dict(configuration))



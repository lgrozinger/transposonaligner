import os
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
        prog="XX",
        description = "",
        epilog = "",
    )

    parser.add_argument("-v", action="count", default=0)

    parser.add_argument("input_dir")

    parser.add_argument(
        "--sequencing-type",
        default="sanger",
        choices=["sanger", "solexa", "illumina"]
    )

    parser.add_argument(
        "--input-type",
        default="fastq",
        choices=["fasta", "fastq"]
    )

    parser.add_argument(
        "--input-ext",
        default="ab1",
    )

    parser.add_argument("output_dir")
    
    parser.add_argument("-o", default="results.csv", metavar="OUTPUT_FILE")

    parser.add_argument(
        "-qc",
        nargs="?",
        const="fastqc",
        default=None,
        metavar="PATH_TO_FASTQC"
    )

    parser.add_argument(
        "-trim",
        nargs="?",
        const="sickle",
        default=None,
        metavar="PATH_TO_SICKLE"
    )

    parser.add_argument("--trim-quality", default=20, type=int)
    parser.add_argument("--trim-length", default=20, type=int)
    parser.add_argument("--trim-save", action="store_true")
    parser.add_argument("-blastn", default="blastn", metavar="PATH_TO_BLASTN")
    parser.add_argument("-transposon", required=True, nargs="+")
    parser.add_argument("--transposon-type", default="genbank")
    parser.add_argument("--transposon-save", action="store_true")
    parser.add_argument("--transposon-word-size", default=10, type=int)
    parser.add_argument("--transposon-evalue", default=0.01, type=float)
    parser.add_argument("-genome", required=True, nargs="+")
    parser.add_argument("--genome-type", default="genbank")
    parser.add_argument("--genome-save", action="store_true")
    parser.add_argument("--genome-word-size", default=10, type=int)
    parser.add_argument("--genome-evalue", default=0.01, type=float)
    parser.add_argument("--genome-prefix", default=9, type=int)
    parser.add_argument("--genome-window", default=18, type=int)
    parser.add_argument("-config")

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

    config["transposon"] = [absolute_path(t, indir) for t in config["transposon"]]

    return config
    
def get_configuration(parser):
    cliargs = vars(parser.parse_args())
    configuration = {}
    
    if cliargs["config"] is not None:
        configuration = load_config(cliargs["config"])
        
    return make_paths_absolute(cliargs | config_to_dict(configuration))



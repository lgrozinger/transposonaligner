#!/usr/local/bin/python3
import argparse
import pandas
import re
import pathlib


def get_args():
    parser = argparse.ArgumentParser(
        prog="tnmeta",
        description = "Attach well metadata to a results table from a plate map"
    )

    parser.add_argument(
        "well_plate_regex",
        help="A regular expression which matches the plates and wells in record names",
    )

    parser.add_argument(
        "result_file",
        help="The spreadsheet to annotate with metadata"
    )

    parser.add_argument(
        "metadata",
        nargs="+",
        help="The spreadsheets containing the metadata map of a plate",
    )

    parser.add_argument(
        "-o",
        default=None,
        metavar="OUTPUT_FILE",
        help="If given, output will be saved to OUTPUT_FILE, if not, RESULT_FILE will be overwritten",
    )

    return parser.parse_args()

class Meta:
    def __init__(self, path_or_filename, regex):
        path = pathlib.Path(path_or_filename).resolve(strict=True)
        path_regex = re.compile(r"^([A-Za-z\d]+)_([A-Za-z\d]+).*")
        self.plate = re.match(path_regex, path.stem).group(1)
        self.column_name = re.match(path_regex, path.stem).group(2)
        self.regex = re.compile(regex.replace("PLATE", self.plate).replace("WELL", r"([A-Za-z0-9]+)"))

        with open(path, "rb") as metadata:
            self.data = pandas.read_excel(metadata, index_col=0)

    @property
    def size(self):
        return (len(self.data), len(self.data.columns))

    def to_alphanumeric(self, i):
        number = (i - 1) // len(self.data) + 1
        letter = self.data.index[(i - 1) % len(self.data)]
        return f"{letter}{number}"

    def get_well_data(self, well):
        alphanum = re.compile(r"([A-Za-z]+)(\d+)")
        numeric = re.compile(r"(\d+)")

        if re.match(alphanum, well) is not None:
            m = re.match(alphanum, well)
            wellletter = m.group(1)
            wellnumber = m.group(2)

            try:
                return self.data.loc[wellletter, wellnumber]
            except KeyError:
                return self.data.loc[wellletter, str(int(wellnumber))]
        elif re.match(numeric, well) is not None:
            wellnumber = int(re.match(numeric, well).group(1))
            return self.get_well_data(self.to_alphanumeric(wellnumber))
        else:
            print(f"The well '{well}' could not be found in the metadata map.")
            return ""

    def add_to_row(self, series):
        print(f"Looking at read {series.name[1]}")
        if re.search(self.regex, series.name[1]) is not None:
            well = re.search(self.regex, series.name[1]).group(1)
            data = self.get_well_data(well)
            print(f"On plate {self.plate} in well {well} adding metadata: {data}")
            series[self.column_name] = data
        else:
            print(f"No metadata for read {series.name[1]}")
        return series
        
    def add_as_column(self, result_table):
        if self.column_name not in result_table:
            result_table.insert(
                len(result_table.columns),
                self.column_name,
                None,
                False
            )
        return result_table.apply(lambda s: self.add_to_row(s), axis=1)
    
def main():
    ARGS = get_args()
    result_path = pathlib.Path(ARGS.result_file).resolve(strict=True)
    with open(result_path, "rb") as results_in:
        result_table = pandas.read_excel(
            results_in,
            index_col=[0, 1],
            dtype={"multiple candidates": bool},
        )

    metadata_files = [pathlib.Path(m).resolve(strict=True) for m in ARGS.metadata]
    metadata = [Meta(metadata_file, ARGS.well_plate_regex) for metadata_file in metadata_files]

    for data in metadata:
        result_table = data.add_as_column(result_table)


    outpath = result_path
    if ARGS.o is not None:
        outpath = pathlib.Path(ARGS.o).resolve()

    with open(outpath, "wb") as results_out:
        result_table.to_excel(results_out)

if __name__ == "__main__":
    main()



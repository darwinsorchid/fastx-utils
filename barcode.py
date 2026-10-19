"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

[file_name]

Description:


User-defined functions:
Non-standard modules:


Procedure:
    1.

Input: FASTQ file
Output:

Usage: python barcode.py input_file.fastq [-o sample]

Version: 1.00
Date:
Name: Alexandra Bekou
--------------------------------

Limitations:
    •

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""

"""This script takes as input a FASTQ file, which contains 4-line records:

@HEADER
SEQUENCE (1 line) 
+ 
QUALITY CHARACTERS (1 line)

Sequence length == Quality line length 
"""

import argparse

from file_io import read_fastq, write_file


def get_args(argv=None):

    parser = argparse.ArgumentParser(
        prog="barcode.py", description="Filter and remove barcodes from FASTQ file."
    )
    parser.add_argument("input_file", help="FASTQ file containing barcoded reads.")
    parser.add_argument(
        "-o",
        "--output_title",
        default="sample",
        help="[optional] Sample name for output files.",
    )

    return parser.parse_args(argv)


def barcode_filter(
    data: dict[str, tuple[str, str]],
    barcode_list=None,
) -> dict[str, list[str]]:

    # Handle barcodes list as optional argument
    if barcode_list is None:
        barcode_list = ["TATCCTCT", "GTAAGGAG", "TCTCTCCG"]

    # Check if provided barcodes list is empty
    elif not barcode_list:
        raise ValueError("Empty barcode list.")

    # Check that all barcodes are of equal length
    len_0 = len(barcode_list[0])
    valid_barcodes = [len(item) == len_0 for item in barcode_list]
    if not all(valid_barcodes):  # Returns True only if every item is True
        raise ValueError(f"Invalid barcode list: {barcode_list}")

    # Initialize sample dictionary with barcodes as keys
    samples_dict = {code: [] for code in barcode_list}
    # Add key for reads without barcode
    samples_dict["undetermined"] = []

    for head, info in data.items():
        # Get DNA chunk from sequence (check both ends)
        seq = info[0]
        start = len(seq) - len(barcode_list[0])  # check 3' end
        end = len(barcode_list[0])  # check 5' end
        chunk1 = seq[start:]
        chunk2 = seq[:end]

        # Check if chunk is exact match with any of the barcodes
        # Store samples with right format to write to output FASTQ file
        if chunk1 in barcode_list:
            samples_dict[chunk1].append(
                f"{head}\n{seq[:start]}\n+\n{info[1][:start]}\n"
            )
        elif chunk2 in barcode_list:
            samples_dict[chunk2].append(f"{head}\n{seq[end:]}\n+\n{info[1][end:]}\n")
        else:
            samples_dict["undetermined"].append(f"{head}\n{seq}\n+\n{info[1]}\n")

    return samples_dict


def write_files(data: dict[str, list[str]], out_name: str = "sample") -> None:

    for i, (code, samples_list) in enumerate(data.items()):
        if code == "undetermined":
            write_file(file_name="undetermined.fastq", lines=samples_list)
        else:
            write_file(file_name=f"{out_name}{i + 1}.fastq", lines=samples_list)


def main(argv=None):

    # Get arguments
    args = get_args(argv)

    # Read FASTQ file
    data_dict = read_fastq(args.input_file)

    # Filter samples
    samples_dict = barcode_filter(data_dict)

    # Write to files
    write_files(data=samples_dict, out_name=args.output_title)


if __name__ == "__main__":
    main()

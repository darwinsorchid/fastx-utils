"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

barcode.py

Description: This script filters barcoded reads in a FASTQ file and writes unique samples into different
            output files, including sequences with no matches for the given barcodes, named as "undetermined"
            reads. The barcode search looks for exact matches and covers both 5' and 3' ends of the FASTQ reads.
            Both the sequences containing barcode matches and their corresponding quality lines are trimmed to
            remove the oligonucleotide before writing to the output files. Undetermined sequences are written
            as is.

User-defined functions: None
Non-standard modules:   None
Internal modules consumed: file_io.py


Procedure:
    1. Parse command-line arguments.
    2. Read FASTQ file and store records as a dictionary of {header: (sequence, quality)}.
    3. Search sequences (5' and 3' ends) for exact matches to given barcodes.
    4. Write matches to each barcode in a different output FASTQ file after trimming sequence
    and quality lines, including non-matches which are written to "undetermined.fastq" untrimmed.


Input: FASTQ file containing DNA sequences
Output: FASTQ files [maximum (barcodes + 1) per input file]

Usage: python barcode.py input_file.fastq [-o sample]


Version: 1.00
Date: 2026-10-20
Name: Alexandra Bekou

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""

# Import modules
import argparse
import os

from file_io import read_fastq, write_file


def get_args(argv=None) -> argparse.Namespace:
    """Parse command-line arguments:
        • input_file: Path to input FASTQ file
        • -o: (optional) Title to be used for output files naming
    Args:
        - argv: (optional) If not None, the passed-in list of arguments is used
                instead of command-line args - used for testing.

    Returns:
        - List of arguments passed in as command-line args or directly passed in
          the function during testing.
    """

    # Create parser instance for script
    parser = argparse.ArgumentParser(
        prog="barcode.py", description="Filter and remove barcodes from FASTQ file."
    )
    # Required input FASTQ file
    parser.add_argument("input_file", help="FASTQ file containing barcoded reads.")

    # Optional title for output files
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
    """Filter FASTQ records for barcode matches.

    Search 5' and 3' end of FASTQ reads for exact barcode matches and trim barcodes from DNA sequence
    and the corresponding quality line before storing. FASTQ records with no matches are marked as
    "undetermined" and are stored untrimmed.

    Args:
        - data: Dictionary containing FASTQ records as {header: (sequence, quality)}.
        - barcode_list: (optional) List of barcodes to be used as query for matches in FASTQ reads.

    Returns:
        - Dictionary containing FASTQ records categorized as different samples (one for every barcode)
        stored as : {barcode: [header, sequence, +, quality]}. Dictionary values are already formatted
        to be written in output files in FASTQ format.

    Raises:
        - ValueError: when an empty barcode list is given.
        - ValueError: when an invalid barcode list is given (unequal lengths).
    """

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
        start = len(seq) - len(barcode_list[0])
        end = len(barcode_list[0])
        chunk1 = seq[start:]
        chunk2 = seq[:end]

        # Check if chunk is exact match with any of the barcodes
        # Store samples with right format to write to output FASTQ file
        # Check 3' end
        if chunk1 in barcode_list:
            samples_dict[chunk1].append(
                f"{head}\n{seq[:start]}\n+\n{info[1][:start]}\n"
            )
        # Check 5' end
        elif chunk2 in barcode_list:
            samples_dict[chunk2].append(f"{head}\n{seq[end:]}\n+\n{info[1][end:]}\n")
        # Handle non-matches
        else:
            samples_dict["undetermined"].append(f"{head}\n{seq}\n+\n{info[1]}\n")

    return samples_dict


def write_files(data: dict[str, list[str]], out_name: str = "sample") -> None:
    """Write filtered FASTQ records to output files.

    Trimmed FASTQ records are split into different files, each holding a set of different
    sample/barcoded records. Undetermined FASTQ records are written to a separate file.
    Maximum number of output files per input file: number of barcodes + 1.

    Args:
        - data: Dictionary containing trimmed and undetermined FASTQ records as:
                {barcode: (header, sequence, +, quality)}.

    Returns:
        - None.

    """

    for i, (code, samples_list) in enumerate(data.items()):
        # Write undetermined FASTQ records to a separate file
        if code == "undetermined":
            write_file(file_name="undetermined.fastq", lines=samples_list)

        # Write each sample to a different file
        else:
            write_file(file_name=f"{out_name}{i + 1}.fastq", lines=samples_list)


def main(argv=None) -> None:
    """Execute main logic.

    Parse arguments, read and filter FASTQ records from input file and create FASTQ output files
    (max barcodes + 1 per input file).

    Args:
        - argv: (optional) List of arguments for input file and title of output files.

    Returns:
        - None.

    Raises:
        - FileNonFoundError: if FASTQ input file path does not exist.
    """
    # Get arguments
    args = get_args(argv)

    # Check that input file exists
    if not os.path.isfile(args.input_file):
        raise FileNotFoundError(f"File {args.input_file} does not exist.")

    # Read FASTQ file
    data_dict = read_fastq(args.input_file)

    # Filter samples
    samples_dict = barcode_filter(data_dict)

    # Write to files
    write_files(data=samples_dict, out_name=args.output_title)


if __name__ == "__main__":
    main()

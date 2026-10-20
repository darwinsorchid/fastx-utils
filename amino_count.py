"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

amino_count.py

Description: This script calculates absolute amino acid abundance from a set of amino acid sequences in
            a FASTA file and creates a tab-separated file with amino acids and their absolute abundance.
            Non-valid amino acids are counted as a separate category and denoted by "X" in the output.

User-defined functions: None
Non-standard modules: None
Internal modules consumed: file_io.py


Procedure:
    1. Parse command-line arguments
    2. Read FASTA file and store records as dictionary: {header: sequence}
    3. Iterate over amino acid sequences and store absolute amino acid abundance as dictionary: {AA: abundance}
    4. Write tab-separated output file with AA and absolute abundance columns


Input: FASTA file containing amino acid sequences
Output: Tab-separated file with amino acid in single-letter format and their absolute abundance including
        undetermined non-valid amino acids denoted by "X"

Usage: python amino_count.py amino.faa [-o output_file.txt]


Version: 1.00
Date: 2026-10-20
Name: Alexandra Bekou

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""

# Import modules
import argparse
import os

from file_io import read_file, write_file


def get_args(argv=None) -> argparse.Namespace:
    """Parse command-line arguments:
        • input_file: Path to input FASTA file
        • -o: (optional) Output file name
    Args:
        - argv: (optional) If not None, the passed-in list of arguments is used
                instead of command-line args - used for testing.

    Returns:
        - List of arguments passed in as command-line args or directly passed in
          the function during testing.
    """

    # Create parser instance for script
    parser = argparse.ArgumentParser(
        prog="amino_count.py",
        description="Count absolute amino acid abundances in a FASTA file with translated sequences.",
    )
    # Required input FASTA file
    parser.add_argument(
        "input_file", help="FASTA file containing translated polypeptides."
    )
    # Optional output file name
    parser.add_argument(
        "output_file",
        nargs="?",
        default="output_file.txt",
        help="[optional] Name for output file.",
    )

    return parser.parse_args(argv)


def count_aa(sequences: list[str]) -> dict[str, int]:
    """Count abundance of amino acids in multiple protein sequences.

    Args:
        - sequences: list of amino acid sequences.

    Returns:
        - Dictionary of amino acids and their abundances in all the given sequences.
    """

    # Define valid amino acids
    valid_aa = "ACDEFGHIKLMNPQRSTVWY"

    # Initialize dictionary with valid amino acids as keys and 0 initial values
    aa_count = dict.fromkeys(valid_aa, 0)

    # Add key for invalid amino acid count
    aa_count["X"] = 0

    for seq in sequences:
        # Handle lowercase letters
        seq = seq.upper()

        for a in seq:
            # Count occurences for each amino acid in all sequences
            if a in aa_count and a != "X":
                aa_count[a] += 1
            else:
                # Count non-valid amino acids
                aa_count["X"] += 1

    return aa_count


def format_counts(counts: dict[str, int]) -> list[str]:
    """Turn {"A" : 2, ...} into a list of tab-separated lines.

    Args:
        - counts: Dictionary containing amino acids and their absolute abundance in a file's sequences.

    Returns:
        - List of strings with amino acid and abundance tab-separated.
    """

    return [f"{a}\t{count}\n" for a, count in counts.items()]


def main(argv=None) -> None:
    """Execute main logic.

    Parse arguments, read file with FASTA records of amino acid sequences, count absolute abundance
    for each amino acid in all given sequences and write to tab-separated file. Non-valid amino acids
    are counted as a separate undetermined amino acid denoted by "X".

    Args:
        - argv: (optional) List of arguments for input and output file names.

    Returns:
        - None.

    Raises:
        - ValueError: if empty sequence list created from FASTA records dictionary.
        - FileNonFoundError: if FASTA input file path does not exist.
    """

    # Get arguments
    args = get_args(argv)

    # Check that input file exists
    if not os.path.isfile(args.input_file):
        raise FileNotFoundError(f"File {args.input_file} does not exist.")

    # Read FASTA file with amino acid sequences
    fasta_dict = read_file(args.input_file)

    # Store sequences to a list
    seq_list = [seq for seq in fasta_dict.values()]

    # Check that sequences exist in list
    if not seq_list:
        raise ValueError("Empty sequence list.")

    # Count absolute abundances for all sequences in the file
    aa_dict = count_aa(seq_list)

    # Format amino acid count dictionary for writing to file
    aa_list = format_counts(aa_dict)

    # Write to output file
    write_file(args.output_file, aa_list)


if __name__ == "__main__":
    main()

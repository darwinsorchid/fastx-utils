"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

[file_name]

Description:


User-defined functions:
Non-standard modules:


Procedure:
    1.

Input 1:
Input 2:
Output:

Usage: python amino_count.py amino.faa output_file.txt [optional]

Version: 1.00
Date:
Name: Alexandra Bekou
--------------------------------

Limitations:
    •

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""

import argparse

from file_io import read_file, write_file


def get_args(argv=None):

    parser = argparse.ArgumentParser(
        prog="amino_count.py",
        description="Count absolute amino acid abundances in a FASTA file with translated sequences.",
    )
    parser.add_argument(
        "input_file", help="FASTA file containing translated polypeptides."
    )
    parser.add_argument(
        "-o",
        "--output_file",
        default="output_file.txt",
        help="[optional] Name for output file.",
    )

    return parser.parse_args(argv)


def count_aa(sequences: list[str]) -> dict[str, int]:
    """Count amino acids in multiple protein sequences.

    Args:
        - sequences: amino acid sequences

    Returns:
        - Dictionary of amino acids and their abundances in the sequences.
    """

    valid_aa = "ACDEFGHIKLMNPQRSTVWY"

    # Initialize dictionary with valid amino acids as keys and 0 initial values
    aa_count = dict.fromkeys(valid_aa, 0)

    # Add key for invalid amino acid count
    aa_count["X"] = 0

    for seq in sequences:
        # Handle lowercase letters
        seq = seq.upper()

        for a in seq:
            if a in aa_count and a != "X":
                aa_count[a] += 1
            else:
                aa_count["X"] += 1

    return aa_count


def format_counts(counts: dict[str, int]) -> list[str]:
    """Turn {"A" : 2, ...} into a list of tab-separated lines.

    Args:
        - counts: Dictionary containing amino acids and their absolute abundance in a file's sequences

    Returns:
        - list of strings with amino acid and abundance tab-separated
    """

    return [f"{a}\t{count}\n" for a, count in counts.items()]


def main(argv=None) -> None:

    # Get arguments
    args = get_args(argv)

    fasta_dict = read_file(args.input_file)

    seq_list = []

    for seq in fasta_dict.values():
        seq_list.append(seq)

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

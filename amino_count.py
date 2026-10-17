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

Usage:

Version: 1.00
Date:
Name: Alexandra Bekou
--------------------------------

Limitations:
    •

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""

from file_io import read_file, write_file


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


def main(in_file: str, out_file: str = "output.txt") -> None:

    fasta_dict = read_file(in_file)

    seq_list = []

    for seq in fasta_dict.values():
        seq_list.append(seq)

    # Count absolute abundances for all sequences in the file
    aa_dict = count_aa(seq_list)

    # Format amino acid count dictionary for writing to file
    aa_list = format_counts(aa_dict)

    # Write to output file
    write_file(out_file, aa_list)


if __name__ == "__main__":
    seqs = ["EHCFYPGGVWLLLLVVVPHENQVSPT", "ARVPKATRWPFRCSPPICIT*FDPVL"]
    in_file = "out_test.fna"
    main(in_file)

    aa_count = count_aa(seqs)

    print(aa_count)

"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

file_reader.py

Description: This script is designed to work as a file-reading module consumed by all other scripts that
            perform analyses on DNA sequences. It parses through a FASTA file and outputs a dictionary of
            header: sequence key-value pairs. It accounts for sequences wrapping over one or multiple lines.
            Since this is developed to be used as a module, all error handling and input QC are to be handled
            in the main script.

User-defined functions: None
Non-standard modules:   None


Procedure:
    1. Parse FASTA file (manages interleaved sequences)
    2. Store header: sequence pairs in dictionary


Input:  File name (str)
Output: Dictionary of header: sequence pairs

Usage: (internal) from file_reader import *

Version: 1.00
Date: 2026-10-17
Name: Alexandra Bekou

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""


def read_file(file_name: str) -> dict:
    """Read a FASTA file and return dictionary with headers:sequences.
    Accounts for sequences wrapping over more than one lines.

    Args:
        - file_name: String containing name of the file to be read

    Returns:
        - Dictionary populated by header: sequence from parsed FASTA file
    """

    with open(file_name, "r") as f:
        # Initialize objects for storing values when parsing FASTA lines
        fasta_dict = {}  # for: Whole FASTA content (header:sequence)
        closed_headers = []  # for: Headers with their whole sequence added to the dictionary
        open_header = (
            ""  # for: Header with still incomplete sequence due to multi-line wrapping
        )
        seqs = []  # for: Subsequences for each multi-line wrapping sequence

        for line in f:
            # Manage header lines
            if line.startswith(">"):
                # I. Already existing open_header (not empty string)
                if line != open_header and open_header != "":
                    # Add header with complete sequence to dictionary as key:value pair
                    fasta_dict[open_header.rstrip("\n")] = "".join(seqs)
                    # Consider header with complete sequence as closed
                    closed_headers.append(open_header)
                    # Initialize empty sequence list for new header
                    seqs = []
                    # Update open header to new header
                    open_header = line

                # II. First header
                else:
                    open_header = line

            # Manage sequence lines; add to subsequence list without newline characters
            else:
                if open_header:
                    seqs.append(line.rstrip("\n"))

    return fasta_dict


if __name__ == "__main__":
    file_name = "DNA.fna"

    fasta_dict = read_file(file_name)
    print(fasta_dict)

"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

file_io.py

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

import os


def read_file(file_name: str) -> dict[str, str]:
    """Read a FASTA file and return dictionary with headers:sequences.
    Accounts for sequences wrapping over more than one lines.

    Args:
        - file_name: String containing name of the file to be read

    Returns:
        - Dictionary populated by header: sequence from parsed FASTA file
    """

    if not os.path.isfile(file_name):
        raise FileNotFoundError(f"File {file_name} does not exist.")

    with open(file_name, "r") as f:
        # Initialize objects for storing values when parsing FASTA lines
        fasta_dict = {}  # for: Whole FASTA content (header:sequence)
        closed_headers = []  # for: Headers with their whole sequence added to the dictionary
        open_header = (
            ""  # for: Header with still incomplete sequence due to multi-line wrapping
        )
        seqs = []  # for: Subsequences for each multi-line wrapping sequence

        for line in f:
            # Extract clean line
            line = line.strip()
            # Skip blank lines
            if not line:
                continue
            # Manage header lines
            if line.startswith(">"):
                # I. Already existing open_header (not empty string)
                if line != open_header and open_header:
                    # Check is sequence for header exists
                    if not seqs:
                        raise ValueError(f"Header {open_header} has no sequence.")
                    else:
                        # Add header with complete sequence to dictionary as key:value pair
                        fasta_dict[open_header] = "".join(seqs)
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
                    seqs.append(line)

        if not open_header:
            raise ValueError(f"{file_name}: No FASTA header found.")

        # Add last header record
        if open_header in fasta_dict:
            raise ValueError("Duplicate header found.")
        else:
            fasta_dict[open_header] = "".join(seqs) if open_header != "" else None

    return fasta_dict


def write_file(file_name: str, lines) -> None:
    """Iterate over provided lines and write them to specified output file.

    Args:
        - file_name: Name of file to be used for writing and saving output
        - lines: Iterable string object containing generated data

    Returns:
        - None
    """

    if type(lines) not in [list, str, dict]:
        raise TypeError("Unsupported data type.")

    with open(file_name, "w") as f:
        if type(lines) in [list, str]:
            f.writelines(lines)
        elif type(lines) == dict:
            f.writelines(f"{key}{value}" for key, value in lines.items())


def read_fastq(file_name: str) -> dict[str, tuple[str, str]]:

    if not os.path.isfile(file_name):
        raise FileNotFoundError(f"File {file_name} does not exist.")

    fastq_dict = {}

    with open(file_name, "r") as f:
        for n, line in enumerate(f):
            line = line.rstrip("\n")
            # Check iteration position based on the FASTQ 4-line record format
            pos = n % 4

            # Check that line is not empty
            if not line:
                raise ValueError("Empty FASTQ line.")

            else:
                # Store record's header line
                if pos == 0:
                    if not line.startswith("@"):
                        raise ValueError(f"Line {n + 1}: header must start with @")
                    else:
                        header = line

                # Store record's sequence line
                elif pos == 1:
                    if not line:
                        raise ValueError("Sequence line is empty.")
                    else:
                        seq = line

                elif pos == 2:
                    if not line.startswith("+"):
                        raise ValueError(
                            f"Line {n + 1}: record's '+' line missing or wrong file format."
                        )

                # Store record's quality line
                elif pos == 3:
                    quality = line

                    # Check sequence and quality lines lengths match
                    if len(seq) != len(quality):
                        raise ValueError(
                            f"Sequence length does not match quality characters: \n{seq}\n{quality}"
                        )

                    # Save full 4-line record
                    if header in fastq_dict:
                        raise ValueError("Duplicate header found.")
                    else:
                        fastq_dict[header] = (seq, quality)

    return fastq_dict


if __name__ == "__main__":
    file_name = "DNA.fna"

    fasta_dict = read_file(file_name)
    print(fasta_dict)

"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

file_io.py

Description: This script is designed to work as a file-reading module consumed by all other scripts that
            perform analyses on DNA sequences. It contains the following functions:

                • read_file: Parse through a FASTA file and output a dictionary of
                            {header: sequence} key-value pairs. It supports sequences
                            wrapping over one or multiple lines.

                • write_file: Iterate over str, list or dict object values and write them
                            to an output file. Content is expected to be formatted before-hand.
                            Supports FASTA and FASTQ file formats.

                • read_fastq: Parse through a FASTQ file and output a dictionary of
                            {header: (sequence, quality)}.

            Since this is developed to be used as a module, all error handling regarding input QC are to
            be handled in the main script and inputs / outputs are separately defined for each function.


User-defined functions: None
Non-standard modules:   None


Usage: (internal) from file_io import *


Version: 1.00
Date: 2026-10-17
Name: Alexandra Bekou

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""


def read_file(file_name: str) -> dict[str, str]:
    """Iterate over FASTA file and return dictionary with headers:sequences.
    Accounts for sequences wrapping over more than one lines.

    Args:
        - file_name: String containing name of the file to be read

    Returns:
        - Dictionary populated by header: sequence from parsed FASTA file

    Raises:
        - ValueError: if header corresponds to empty sequence
        - ValueError: if file does not contain FASTA headers (content compatibility check)
        - ValueError: if duplicate header is found in file
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

        # Check format - raise if no compatible FASTA header found
        if not open_header:
            raise ValueError(f"{file_name}: No FASTA header found.")

        # Check for duplicate header entries
        if open_header in fasta_dict:
            raise ValueError("Duplicate header found.")
        else:
            # Add last header record
            fasta_dict[open_header] = "".join(seqs) if open_header != "" else None

    return fasta_dict


def write_file(file_name: str, lines) -> None:
    """Iterate over provided lines and write them to specified output file.

    Args:
        - file_name: String containing name of file to be used for writing and saving output
        - lines: Iterable string object containing already formatted content

    Returns:
        - None

    Raises:
        - TypeError: if lines does not correspond to a supported type (list, str, dict)
    """

    # Check that content is stored in a supported data type
    if type(lines) not in [list, str, dict]:
        raise TypeError("Unsupported data type.")

    with open(file_name, "w") as f:
        # Write string / list content
        if type(lines) in [list, str]:
            f.writelines(lines)
        # Write dict content
        elif type(lines) == dict:
            f.writelines(f"{key}{value}" for key, value in lines.items())


def read_fastq(file_name: str) -> dict[str, tuple[str, str]]:
    """Iterate over FASTQ file and store 4-line records as a dictionary.

    Args:
        - file_name: String containing name of file to be read.

    Returns:
        - Dictionary of FASTQ records as: {header: (sequence, quality)}.

    Raises:
        - FileNotFoundError: if input file does not exist
    """

    fastq_dict = {}

    with open(file_name, "r") as f:
        for n, line in enumerate(f):
            # Keep clean info (no newline characters)
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
                        # FASTQ format check
                        raise ValueError(f"Line {n + 1}: header must start with @")
                    else:
                        header = line

                # Store record's sequence line
                elif pos == 1:
                    if not line:
                        raise ValueError("Sequence line is empty.")
                    else:
                        seq = line

                # FASTQ format check
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

                    # Check for duplicate headers
                    if header in fastq_dict:
                        raise ValueError("Duplicate header found.")
                    else:
                        # Save full 4-line record
                        fastq_dict[header] = (seq, quality)

    return fastq_dict

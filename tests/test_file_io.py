"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

test_file_io.py

Description: Driver script for testing file_io.py. Tests the accuracy and robustness of functions:
            - read_file
            - read_fastq
            - write_file


Usage: python -m tests.test_file_io

Version: 1.00
Date: 2026-10-20
Name: Alexandra Bekou
--------------------------------

Limitations:
    • Hard-coded paths for test data
    • Simple print statements instead of pytest

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""

# Import modules
import os

from file_io import read_fastq, read_file, write_file


def test_read_file(in_file):
    """Read a FASTA file and print the {header: sequence} dictionary."""

    try:
        fasta_dict = read_file(in_file)
        print(f"✅ File {in_file} read: \nFASTA content as dict:\n {fasta_dict}.")

    except ValueError as e:
        print(f"✅{in_file}: Error handled -> {e}")


def test_write_file(out_file, lines):
    """Write content to a file and check the file exists."""

    write_file(f"tests/{out_file}", lines)

    if os.path.exists(f"tests/{out_file}"):
        print(f"✅ File saved at: tests/{out_file}.")


def test_read_fastq(in_file):
    """Read a FASTQ file and print the number of records and the first record."""

    try:
        fastq_dict = read_fastq(in_file)
    except ValueError as e:
        print(f"✅{in_file}: Error handled -> {e}")
        # End loop
        return

    # Get first dictionary record to print
    first_header = ""

    for header, (seq, quality) in fastq_dict:
        # Store first header
        first_header = header
        break

    seq, quality = fastq_dict[first_header]

    print(
        f"✅ FASTQ records in {in_file} read. First record: \n{first_header}\n{seq}\n{quality}"
    )


if __name__ == "__main__":
    """Run main test logic."""

    # Test read_file (Handles FASTA files)
    fasta_files = {
        "empty_file": "tests/test_data/empty_file.txt",
        "multi_line_file": "tests/test_data/fasta_multiline_wrapped.txt",
        "single_record_file": "tests/test_data/fasta_single_record.txt",
        "fasta_trailing_blank": "tests/test_data/fasta_trailing_blank_line.txt",
    }

    for path in fasta_files.values():
        test_read_file(path)

    # Test read_fastq (Handles FASTQ files)
    fastq_files = {
        "empty_file": "tests/test_data/empty_file.txt",
        "fastq_length_mismatch": "tests/test_data/fastq_seq_qual_length_mismatch.txt",
        "fastq_trailing_blank": "tests/test_data/fastq_trailing_blank_line.txt",
        "valid_2_reads": "tests/test_data/fastq_valid_two_reads.txt",
    }

    for path in fastq_files.values():
        test_read_fastq(path)

    # Test write_file
    lines_str = ">seq1\nAAGCTNATNAT\n>seq2\nGCGATAGCGGTATTATAT\n"
    lines_dict = {">seq1\n": "AAGCTNATNAT\n", ">seq2\n": "GCGATAGCGGTATTATAT\n"}
    lines_list = [">seq1\nAAGCTNATNAT\n", ">seq2\nGCGATAGCGGTATTATAT\n"]
    out_file_str = "test_output_str.txt"
    out_file_list = "test_output_list.txt"
    out_file_dict = "test_output_dict.txt"

    test_write_file(out_file_str, lines_str)
    test_write_file(out_file_list, lines_list)
    test_write_file(out_file_dict, lines_dict)

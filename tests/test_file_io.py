"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

[file_name]

Description:


User-defined functions:
Non-standard modules:


Procedure:
    1.

Input:
Output:

Usage: python -m tests.test_file_io

Version: 1.00
Date: 2026-10-17
Name: Alexandra Bekou
--------------------------------

Limitations:
    •

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""

import os

from file_io import read_file, write_file


def test_read_file(in_file):

    fasta_dict = read_file(in_file)

    print(f"✅ File read: \nFASTA content as dict:\n {fasta_dict}")


def test_write_file(out_file, lines):

    write_file(f"tests/{out_file}", lines)

    if os.path.exists(f"tests/{out_file}"):
        print(f"✅ File saved at: tests/{out_file}")


if __name__ == "__main__":
    in_file = "DNA.fna"
    lines_str = ">seq1\nAAGCTNATNAT\n>seq2\nGCGATAGCGGTATTATAT\n"
    lines_dict = {">seq1\n": "AAGCTNATNAT\n", ">seq2\n": "GCGATAGCGGTATTATAT\n"}
    lines_list = [">seq1\nAAGCTNATNAT\n", ">seq2\nGCGATAGCGGTATTATAT\n"]
    out_file_str = "test_output_str.txt"
    out_file_list = "test_output_list.txt"
    out_file_dict = "test_output_dict.txt"

    test_read_file(in_file)
    test_write_file(out_file_str, lines_str)
    test_write_file(out_file_list, lines_list)
    test_write_file(out_file_dict, lines_dict)

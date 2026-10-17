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

import dna2aa as transl


def test_translation(seq):

    polypeptide = transl.main(seq)

    print(f"Original sequence: {seq}\nTranslated sequence: {polypeptide}")


if __name__ == "__main__":
    seq = "ATAGTCGTGaTGTA"
    false_seq = "ATAGTCGTAaNnAT"
    rna_seq = "AUAGUCGUGAUGUA"

    true_protein = test_translation(
        seq
    )  # Expected output: printed original & translated sequence
    rna_protein = test_translation(
        rna_seq
    )  # Expected output: ValueError (Expected DNA got RNA)
    false_protein = test_translation(
        false_seq
    )  # Expected output: ValueError (invalid codon found)

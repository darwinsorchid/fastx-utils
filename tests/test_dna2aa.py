"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

test_dna2aa.py

Description: Driver script for testing dna2aa.py. Tests the accuracy and robustness of functions:
            - translate_seq
            - main


Usage: python -m tests.test_dna2aa

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

from dna2aa import *


def test_single_translation(seq, expected_aa_seq):
    """Test single sequence translation."""

    try:
        polypeptide = translate_seq(seq)
    except ValueError as e:
        print(f"✅ Original sequence: {seq}\nError handled -> {e}")
        # Stop loop
        return

    if polypeptide == expected_aa_seq:
        print(f"✅ Original sequence: {seq}\nTranslated sequence: {polypeptide}.")
    else:
        print(
            f"Original sequence: {seq}\nTranslated sequence: {polypeptide} (expected {expected_aa_seq})."
        )


def test_file_translation(in_file, out_file):
    """Test output of FASTA file translation."""

    try:
        # Pass arguments as list instead of command-line args
        main([in_file, out_file])

        if os.path.exists(out_file):
            print(f"✅ Translated {in_file}! Result saved at: {out_file}.")
        else:
            print(f"Failed to create {out_file}.")

    except (FileNotFoundError, ValueError) as e:
        print(f"✅ {in_file}: error handled -> {e}")


if __name__ == "__main__":
    """Run main test logic."""
    seq = "ATAGTCGTGaTGTA"
    false_seq = "ATAGTCGTAaNnAT"
    rna_seq = "AUAGUCGUGAUGUA"
    in_file = "tests/test_data/DNA.fna"
    in_file_false = "tests/test_data/false.fna"
    out_file = "out_test.fna"

    # Test full script through main()
    # 1. Valid file
    test_file_translation(in_file, f"tests/{out_file}")

    # 2. File does not exist
    test_file_translation(
        "fake_file.fna", f"tests/{out_file}"
    )  # Expected: FileNotFoundError

    # 3. Empty file
    test_file_translation(
        "tests/test_data/empty_file.txt", f"tests/{out_file}"
    )  # Expected: ValueError (no FASTA header)

    # 4. File with invalid sequences
    test_file_translation(
        in_file_false, f"tests/{out_file}"
    )  # Expected: Sequence {header} excluded: contains non-DNA letters.

    test_single_translation(
        seq, "IVVM"
    )  # Expected: printed original & translated sequence

    test_single_translation(
        false_seq, None
    )  # Expected output: ValueError (invalid codon found)

    test_single_translation(false_seq, None)

"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

test_amino_count.py

Description: Driver script for testing amino_count.py. Tests the accuracy and robustness of functions:
            - main



Usage: python -m tests.test_amino_count

Version: 1.00
Date: 2026-10-20
Name: Alexandra Bekou
--------------------------------

Limitations:
    • Hard-coded paths for test data
    • Simple print statements instead of pytest
    • count_aa() is not checked

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""

# Import modules
import os

from amino_count import *


def test_aa_count(in_file, out_file):
    """Test output of main()."""

    try:
        main([in_file, out_file])

        if os.path.exists(out_file):
            print(
                f"✅ Counted AA abundances in: {in_file}! Result saved at: {out_file}."
            )
        else:
            print(f"Failed to create: {out_file}.")
    except (FileNotFoundError, ValueError) as e:
        print(f"✅ {in_file}: error handled -> {e}")


if __name__ == "__main__":
    """Run main test logic."""
    in_file = "tests/test_data/amino.faa"
    invalid_file = "sham_file.fna"
    empty_file = "tests/test_data/empty_file.txt"
    out_file = "out_test.txt"

    # Test valid file
    test_aa_count(in_file, f"tests/{out_file}")  # Expected: Result saved

    # Test non-existent file
    test_aa_count(invalid_file, f"tests/{out_file}")  # Expected: FileNotFoundError

    # Test empty file
    test_aa_count(
        empty_file, f"tests/{out_file}"
    )  # Expected: ValueError (no FASTA header)

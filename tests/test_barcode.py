"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

test_barcode.py

Description: Driver script for testing barcode.py. Tests the accuracy and robustness of functions:
            - barcode_filter
            - main


Usage: python -m tests.test_barcode.py

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
from barcode import *


def test_barcode(in_file):
    """Test filtering and trimming of a FASTQ file into sample files."""

    try:
        main([in_file])

        if os.path.exists("sample1.fastq"):
            print(
                f"✅ Filtered {in_file}! Samples saved as: sample1.fastq, sample2.fastq, ..."
            )
        else:
            print("Failed to create sample1.fastq.")

    except (FileNotFoundError, ValueError) as e:
        print(f"✅ {in_file}: error handled -> {e}")


if __name__ == "__main__":
    """Run main test logic."""
    in_file = "tests/test_data/barcode.fastq"
    invalid_file = "sham_file.fastq"
    empty_file = "tests/test_data/empty_file.txt"
    fasta_file = "tests/test_data/fasta_single_record.txt"

    # Test valid file
    test_barcode(in_file)  # Expected: sample files created

    # Test non-existent file
    test_barcode(invalid_file)  # Expected: FileNotFoundError

    # Test empty file
    test_barcode(empty_file)  # Expected: ValueError (no records)

    # Test FASTA file
    test_barcode(fasta_file)  # Expected: ValueError (FASTA, not FASTQ)

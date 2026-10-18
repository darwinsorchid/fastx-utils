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

from amino_count import *


def test_aa_count(in_file, out_file="output.txt"):

    main(in_file)

    if os.path.exists(out_file):
        print(f"✅ Counted AA abundances in: {in_file}! Result saved at: {out_file}.")
    else:
        print(f"Failed to create: {out_file}.")


if __name__ == "__main__":
    in_file = "data/amino.faa"

    test_aa_count(in_file)

"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

dna2aa.py

Description: This script translates multiple DNA sequences extracted from a FASTA file to amino acid
            sequences using the standard genetic code. Translation is performed starting from the
            first base of the sequence (without searching for a start codon). Since sequences are
            split into codons (3-base strings), entries of length that is not a multiple of three
            leave trailing bases at the end. This program ignores 1-2 trailing bases in such cases
            during translation.Stop codons are denoted by "*" and translation continues after a
            stop codon is found. Sequences that contain any non-DNA letters are reported as invalid
            and excluded from the translation workflow.


User-defined functions: None
Non-standard modules: None
Internal modules consumed: file_io.py


Procedure:
    1. Parse command-line arguments
    2. Read FASTA file and store records as dictionary: {header: sequence}
    3. Iterate over sequences
    4. Check if sequence contains only valid DNA nucleotides
    5. Convert T -> U (DNA -> RNA)
    6. Split sequence into codons
    7. Translate each codon to the corresponding amino acid using the standard
      genetic code and store as list
    8. Join amino acids list to polypeptide string
    9. Write {header: AA sequence} to output file in FASTA format

Input: FASTA file with DNA sequences
Output: FASTA file containing

Usage: python dna2aa.py DNA.fna [-o output_file.txt]

Version: 1.00
Date: 2026-10-19
Name: Alexandra Bekou
------------------------------

Limitations:
    - Only standard genetic code implemented (other tables can be supported if added to helpers module)
    - Translation does not search for a start codon
    - Translation is performed using only 1 reading frame (instead of 6)

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""

# Imports
import argparse
import os

from file_io import read_file, write_file
from helpers import *


def get_args(argv=None) -> argparse.Namespace:
    """Parse command-line arguments:
        • input_file: Path to input FASTA file
        • -o: (optional) Output file name
    Args:
        - argv: (optional) If not None, the passed-in list of arguments is used
                instead of command-line args - used for testing.

    Returns:
        - List of arguments passed in as command-line args or directly passed in
          the function during testing.
    """

    # Create parser instance for script
    parser = argparse.ArgumentParser(
        prog="dna2aa.py",
        description="Translates DNA sequences from FASTA file to amino acid sequences.",
    )
    # Required input FASTA file
    parser.add_argument("input_file", help="FASTA file containing DNA sequences.")

    # Optional output file name
    parser.add_argument(
        "-o",
        "--output_file",
        default="output_file.txt",
        help="[optional] Name for output file.",
    )

    return parser.parse_args(argv)


def dna_to_rna(seq: str) -> str:
    """Take in a DNA sequence and replace T with U.

    Args:
        - seq: DNA sequence.

    Returns:
        - DNA sequence with U instead of T.
    """
    return seq.replace("T", "U")


def check_dna(seq: str) -> bool:
    """Check if a given sequence contains any non-DNA letters.

    A boolean variable is used to control if a sequence will be included for
    translation and written to the output file or if it will be reported as
    invalid (containing non-DNA letters) and excluded from the remaining steps
    of the translation workflow.

    Args:
        -seq: DNA sequence in string format.

    Returns:
        - Boolean variable for reporting sequence as invalid or
        including it for translation.
    """

    # Valid nucleotide letters
    nucleotides = "ATGC"

    for n in seq:
        if n not in nucleotides:
            # Break loop if a non-valid letter is found
            include_seq = False
            break
        else:
            include_seq = True

    return include_seq


def seq_to_codons(seq: str) -> list[str]:
    """Split sequence in codons (3-base strings) and return a list of the codons.
    Disregard 1-2 leftover bases at the end if sequence length is not multiple of 3.

    Args:
        -seq: Nucleotide sequence in string format.
    Returns:
        - List of 3-base codons.
    """

    # Initialize empty list to store codons
    codons = []

    # Split sequence in 3-base codon strings
    for i in range(0, len(seq), 3):
        # Skip last 1-2 nucleotides if sequence is not a multiple of 3
        if len(seq[i : i + 3]) < 3:
            continue
        # Store codons
        else:
            codons.append(seq[i : i + 3])

    return codons


def translate_codon(codon: str, translation_table: str = "standard") -> str:
    """Single codon translation.

    Translate a single codon to the corresponding amino acid using the
    specified genetic code, using the standard genetic code as default.

    Args:
        - codon: 3-base string.
        - translation_table: the genetic code selected for performing codon
          to amino acid translation.

    Returns:
        - Translated amino acid in single-letter format.

    Raises:
        - ValueError: if a 3-base string is not a valid codon.
    """

    # Retrieve appropriate genetic code table using helper function
    if translation_table == "standard":
        code = get_code()

    # Check if 3-base string is a valid codon
    if codon not in all_possible_codons():
        raise ValueError("Invalid codon found.")

    # Get amino acid from table for the given codon
    else:
        amino_acid = code["codons"][codon]

    return amino_acid


def codons_to_polypeptide(codons_list: list[str]) -> str:
    """Translate a list of codons to a list of amino acids and convert into string.

    Args:
        - codons_list: list of codons.

    Returns:
        - String of amino acids (polypeptide).
    """

    # Initialize empty list for storing amino acids
    aa_list = []

    # Translate codons one by one and add to list
    for codon in codons_list:
        aa_list.append(translate_codon(codon))

    return "".join(aa_list)


def translate_seq(seq: str) -> str:
    """Translate a single DNA sequence.

    Translation begins from the first codon (no looking for start codons)
    and continues after a stop codon.

    Args:
        -seq: DNA sequence in string format.

    Returns:
        - polypeptide sequence in string format.
    """

    # Account for lowercase nucleotides if any
    seq = seq.upper()

    # DNA -> RNA
    rna_seq = dna_to_rna(seq)

    # Split sequence in codons
    codons = seq_to_codons(rna_seq)

    # Translate sequence to polypeptide
    polypeptide = codons_to_polypeptide(codons)

    return polypeptide


def translate_fasta(fasta_dict: dict[str, str]) -> dict[str, str]:
    """Translate a full FASTA file of DNA sequences and
    save output to file with FASTA headers + polypeptides.

    Args:
        - fasta_dict: FASTA records dictionary with header: DNA sequence.

    Returns:
        - FASTA records dictionary as {header: amino acid sequence}.
    """

    protein_dict = {}

    for head, seq in fasta_dict.items():
        # Check validity of DNA sequence
        report_seq = check_dna(seq)

        # Report any sequences that have ANY non-DNA letter
        if not report_seq:
            print(f"Sequence {head} excluded: contains non-DNA letters.")
            continue

        else:
            # Convert a full sequence into a polypeptide and store as {FASTA header: AA sequence}
            polypeptide = translate_seq(seq)
            protein_dict[f"{head}\n"] = f"{polypeptide}\n"

    return protein_dict


def main(argv=None) -> None:
    """Execute main logic.

    Parse arguments, read file with FASTA records of DNA sequences and translate sequences
    containing only DNA letters. Sequences with any non-DNA letters are reported as invalid
    and are not translated. Translated sequences are written to an output file in FASTA format
    (no multi-line wrapping).

    Args:
        - argv: (optional) List of arguments for input and output file names.

    Returns:
        - None

    Raises:
        - FileNonFoundError: if FASTA input file path does not exist.
    """
    # Get arguments
    args = get_args(argv)

    # Check that input file exists
    if not os.path.isfile(args.input_file):
        raise FileNotFoundError(f"File {args.input_file} does not exist.")

    # Read and store FASTA records
    fasta_dict = read_file(args.input_file)

    # Translate DNA sequences and store as FASTA records
    protein_dict = translate_fasta(fasta_dict)

    # Write to output file with FASTA format (header - translated sequence)
    write_file(args.output_file, protein_dict)


if __name__ == "__main__":
    main()

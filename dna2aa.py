"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

dna2aa.py

Description:


User-defined functions: None
Non-standard functions: None


Procedure:


Input 1: FASTA file with DNA sequences
Output: FASTA file format headers + translated sequences

Usage: python dna2aa.py DNA.fna output_file.txt [optional]

Version: 1.00
Date: 2026-10-19
Name: Alexandra Bekou
------------------------------


Limitations:


Convert T→U and translate from the first
base (no start codon search).

Stop codons appear as *, and translation
continues after them.
•Ignore 1-2 leftover bases at the end.
•Report any sequence with non-DNA letters.
•Create your own DNA.fna to test with.
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""

# Imports
import argparse

from file_io import read_file, write_file
from helpers import *


def get_args(argv=None):

    parser = argparse.ArgumentParser(
        prog="dna2aa.py",
        description="Translates DNA sequences from FASTA file to amino acid sequences.",
    )
    parser.add_argument("input_file", help="FASTA file containing DNA sequences.")
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
        - seq: DNA sequence

    Returns:
        - DNA sequence with U instead of T
    """
    return seq.replace("T", "U")


def check_dna(seq: str) -> bool:
    """Check if a given sequence contains non-DNA letters.

    Args:
        -seq: DNA sequence in string format

    Returns:
        - boolean variable for reporting sequence as invalid or including it for translation.
    """

    # Valid nucleotide letters
    nucleotides = "ATGC"

    for n in seq:
        # Break loop if a non-valid letter is found
        if n not in nucleotides:
            include_seq = False
            break
        else:
            include_seq = True

    return include_seq


def seq_to_codons(seq: str) -> list[str]:
    """Split sequence to codons (3-base strings) and return a list of the codons.
    Disregard 1-2 leftover bases at the end if number of bases is not divisible by 3.

    Args:
        -seq: Nucleotide sequence in string format
    Returns:
        - List of 3-base codons
    """

    # Initialize empty list to store codons
    codons = []

    # Split sequence in 3-base strings
    for i in range(0, len(seq), 3):
        if len(seq[i : i + 3]) < 3:
            continue
        else:
            codons.append(seq[i : i + 3])

    return codons


def translate_codon(codon: str, translation_table: str = "standard") -> str:
    """
    Translate a single codon to the corresponding amino acid using the specified genetic code.

    Args:
        - codon: 3-base string
        - translation_table: the genetic code selected for performing codon to amino acid translation

    Returns:
        - translated amino acid in single-letter format
    """

    if translation_table == "standard":
        code = get_code()

    if codon not in all_possible_codons():
        raise ValueError("Invalid codon found.")

    else:
        amino_acid = code["codons"][codon]

    return amino_acid


def codons_to_polypeptide(codons_list: list[str]) -> str:
    """Translate a list of codons to a list of amino acids and convert that to a string.

    Args:
        - codons_list: List of codons

    Returns:
        - string of amino acids (polypeptide)
    """

    # Initialize empty list for storing amino acids
    aa_list = []

    for codon in codons_list:
        aa_list.append(translate_codon(codon))

    return "".join(aa_list)


def translate_seq(seq: str) -> str:
    """Translate one DNA sequence.

    Args:
        -seq: DNA sequence in string format

    Returns:
        - polypeptide sequence in string format

    Raises:
        - ValueError: when sequence is RNA instead of DNA
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
        - fasta_dict: FASTA records dictionary with header: DNA sequence

    Returns:
        - FASTA records dictionary with header: amino acid sequence
    """

    protein_dict = {}

    for head, seq in fasta_dict.items():
        # Check validity of DNA sequence
        report_seq = check_dna(seq)

        # Check whether sequence is actually DNA
        if "U" in seq:
            raise ValueError("Expected DNA sequence, got RNA.")

        elif not report_seq:
            print(f"Sequence {head} excluded: contains non-DNA letters.")
            continue

        else:
            polypeptide = translate_seq(seq)
            protein_dict[f"{head}\n"] = f"{polypeptide}\n"

    return protein_dict


def main(argv=None):

    # Get arguments
    args = get_args(argv)

    # Read and store FASTA records
    fasta_dict = read_file(args.input_file)

    # Translate DNA sequences and store as FASTA records
    protein_dict = translate_fasta(fasta_dict)

    # Write to output file with FASTA format (header - translated sequence)
    write_file(args.output_file, protein_dict)


if __name__ == "__main__":
    main()

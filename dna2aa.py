"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

dna2aa.py

Description:


User-defined functions: None
Non-standard functions: None


Procedure:


Input 1: FASTA file with DNA sequences
Output: FASTA file format headers + translated sequences

Usage: python dna2aa.py DNA.fna [output_file.txt]

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
from helpers import _all_possible_codons, _get_code


def dna_to_rna(seq: str) -> str:
    """Takes in a DNA sequence and replaces T with U.

    Args:
        - seq: DNA sequence

    Returns:
        - DNA sequence with U instead of T
    """
    return seq.replace("T", "U")


def check_dna(seq: str) -> bool:
    """Checks if a given sequence contains non-DNA letters.

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
            print("SequenceExcluded: contains non-DNA letters.")
            break
        else:
            include_seq = True

    return include_seq


def seq_to_codons(seq: str) -> list[str]:
    """Splits sequence to codons (3-base strings) and returns a list of the codons.
    Disregards 1-2 leftover bases at the end if number of bases is not divisible by 3.

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
    Translates a single codon to the corresponding amino acid using the specified genetic code.

    Args:
        - codon: 3-base string
        - translation_table: the genetic code selected for performing codon to amino acid translation

    Returns:
        - translated amino acid in single-letter format
    """

    if codon not in _all_possible_codons():
        print("Invalid codon found.")

    if translation_table == "standard":
        code = _get_code()

    amino_acid = code["codons"][codon]

    return amino_acid


def codons_to_polypeptide(codons_list: list[str]) -> str:
    """Translates a list of codons to a list of amino acids and converts that to a string.

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


def main(seq):
    """Executes the basic logic of the program.

    Args:
        -seq: DNA sequence in string format

    Returns:
        - polypeptide sequence in string format

    Raises:
        - ValueError: when sequence is RNA instead of DNA
    """

    # Account for lowercase nucleotides if any
    seq = seq.upper()

    # Check whether sequence is actually DNA
    if "U" in seq:
        raise ValueError("Expected DNA sequence, got RNA.")

    # DNA -> RNA
    rna_seq = dna_to_rna(seq)

    # Split sequence in codons
    codons = seq_to_codons(rna_seq)

    # Translate sequence to polypeptide
    polypeptide = codons_to_polypeptide(codons)

    return polypeptide


if __name__ == "__main__":
    seq = "ATAGTCGTGaTGTA"
    false_seq = "ATAGTCGTAaNnAT"
    rna_seq = "AUAGUCGUGAUGUA"

    # rna_seq = dna_to_rna(seq)
    # print(rna_seq)

    # codon = "UGA"
    # amino_acid = translate_codon(codon)
    # print(amino_acid)

    # codons = seq_to_codons(seq)
    # print(codons)

    # polypeptide = codons_to_polypeptide(codons)
    # print(polypeptide)

    polypeptide = main(seq)
    print(polypeptide)

    rna_polypeptide = main(rna_seq)

    report_seq_true = check_dna(seq)
    print(report_seq_true)
    report_seq_false = check_dna(false_seq)
    print(report_seq_false)


# Check that it's DNA --> Function #1
# Check if it's RNA or non-DNA letters --> Report sequence

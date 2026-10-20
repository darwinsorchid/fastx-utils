"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

helpers.py

Description: This script is made for defining helper functions and keeping code modular. It contains the
            following functions:

                • get_code: Returns the standard genetic code (NCBI Table 1) as a dictionary with start
                codons, stop codons and a codon → amino acid look-up table.


                • all_possible_codons: Returns a set of all 64 possible RNA codons, used to validate
                codons before translation.


User-defined functions: None
Non-standard modules: None


Usage: (internal) from helpers import *


Version: 1.00
Date:
Name: Alexandra Bekou
--------------------------------

Limitations:
    • Only standard genetic code is supported
    • Codons must be uppercase RNA

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""


def get_code():
    """Helper function to get the genetic code as a dictionary"""

    # NCBI Table 1
    standard_gen_code = {
        "start": ["AUG", "UUG", "CUG"],
        "stop": ["UAA", "UAG", "UGA"],
        "codons": {
            "AAA": "K",
            "AAU": "N",
            "AAG": "K",
            "AAC": "N",
            "AUA": "I",
            "AUU": "I",
            "AUG": "M",
            "AUC": "I",
            "AGA": "R",
            "AGU": "S",
            "AGG": "R",
            "AGC": "S",
            "ACA": "T",
            "ACU": "T",
            "ACG": "T",
            "ACC": "T",
            "UAA": "*",
            "UAU": "Y",
            "UAG": "*",
            "UAC": "Y",
            "UUA": "L",
            "UUU": "F",
            "UUG": "L",
            "UUC": "F",
            "UGA": "*",
            "UGU": "C",
            "UGG": "W",
            "UGC": "C",
            "UCA": "S",
            "UCU": "S",
            "UCG": "S",
            "UCC": "S",
            "GAA": "E",
            "GAU": "D",
            "GAG": "E",
            "GAC": "D",
            "GUA": "V",
            "GUU": "V",
            "GUG": "V",
            "GUC": "V",
            "GGA": "G",
            "GGU": "G",
            "GGG": "G",
            "GGC": "G",
            "GCA": "A",
            "GCU": "A",
            "GCG": "A",
            "GCC": "A",
            "CAA": "Q",
            "CAU": "H",
            "CAG": "Q",
            "CAC": "H",
            "CUA": "L",
            "CUU": "L",
            "CUG": "L",
            "CUC": "L",
            "CGA": "R",
            "CGU": "R",
            "CGG": "R",
            "CGC": "R",
            "CCA": "P",
            "CCU": "P",
            "CCG": "P",
            "CCC": "P",
        },
    }

    return standard_gen_code


def all_possible_codons():
    """Helper function that returns a set of all possible codons
    Args:
        - None
    Returns:
        - Set of all possible codons
    """
    return {
        "AAA",
        "AAC",
        "AAG",
        "AAU",
        "ACA",
        "ACC",
        "ACG",
        "ACU",
        "AGA",
        "AGC",
        "AGG",
        "AGU",
        "AUA",
        "AUC",
        "AUG",
        "AUU",
        "CAA",
        "CAC",
        "CAG",
        "CAU",
        "CCA",
        "CCC",
        "CCG",
        "CCU",
        "CGA",
        "CGC",
        "CGG",
        "CGU",
        "CUA",
        "CUC",
        "CUG",
        "CUU",
        "GAA",
        "GAC",
        "GAG",
        "GAU",
        "GCA",
        "GCC",
        "GCG",
        "GCU",
        "GGA",
        "GGC",
        "GGG",
        "GGU",
        "GUA",
        "GUC",
        "GUG",
        "GUU",
        "UAA",
        "UAC",
        "UAG",
        "UAU",
        "UCA",
        "UCC",
        "UCG",
        "UCU",
        "UGA",
        "UGC",
        "UGG",
        "UGU",
        "UUA",
        "UUC",
        "UUG",
        "UUU",
    }

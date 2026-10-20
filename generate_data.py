"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

generate_data.py

Description: This script generates a single file containing DNA sequences of various lengths in
             FASTA format. It also gives the option of creating sequences with false DNA letters
             or RNA nucleotides to some of the sequences to test the robustness of a program's QC.

User-defined functions: None
Non-standard modules: None

Procedure:
    1. Parses command-line arguments for:
        - output file name
        - number of sequences to be generated
        - valid DNA sequences or sequences with non-DNA letters/RNA nucleotides
    2. Uses header & sequence generators to generate headers +  valid/invalid sequences
        of user-specified number and random length
    3. Writes to FASTA file (.fna) with user-defined file name

Input 1: Output file name (e.g. DNA.fna)
Input 2: (optional) Number of sequences to be generated
Input 3: (optional) Perturb / Non-perturb mode

Output: File in FASTA format with generated headers and corresponding sequences

Usage:
    - MODE I. Generate valid DNA sequences:
        python generate_data.py out_file.fna 10

    - MODE II. Generate invalid DNA sequences:
        python generate_data.py out_file.fna 10 -p


Version: 1.00
Date: 2026-10-16
Name: Alexandra Bekou
------------------------

Limitations:
    • Currently perturbation works on all generated sequences with a uniform probability; Perturbing
    a portion of the file's sequences would be more realistic.

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""

# Import modules
import argparse
import random
from collections.abc import Generator


def parse_args(argv=None) -> argparse.Namespace:
    """Parse command-line arguments:
        • file_name: Output file name that will be used to write content.
        • num_seq: Number of sequences to be generated.
        • -p: (optional) If used it stores True as a boolean value and the script produces
            invalid sequences, otherwise it stores False.
    Args:
        - argv: (optional) If not None, the passed-in list of arguments is used
                instead of command-line args - used for testing.

    Returns:
        - List of arguments passed in as command-line args or directly passed in
          the function during testing.
    """
    parser = argparse.ArgumentParser(
        description="Generate FASTA file containing (valid/invalid) DNA sequences for testing.\n"
        "Defaults to generating valid DNA sequences unless otherwise specified by the user."
    )
    parser.add_argument("file_name", type=str, help="Output file name.")
    parser.add_argument(
        "num_seq",
        type=int,
        nargs="?",
        default=10,
        help="Number of sequences generated.",
    )
    parser.add_argument(
        "-p",
        "--perturb",
        action="store_true",
        help="Choose whether to produce valid or invalid DNA sequences for testing.",
    )

    return parser.parse_args(argv)


def header_generator(num_seq: int) -> Generator[str]:
    """Simple generator function for producing FASTA headers.

    Args:
        - num_seq: Number of headers to be generated. Equals to the number of sequences to be generated.

    Yields:
        - The generated header in a single line
    """
    for n in range(num_seq):
        yield f">sequenceID-{(n + 1):03d}"


def seq_generator(num_seq: int, perturb: bool) -> Generator[str]:
    """Simple generator function for producing FASTA sequences.

    Args:
        - num_seq: Number of sequences to be generated. Equals to the number of headers to be generated
        - perturb: The generator's mode; if False it produces valid DNA sequences, otherwise it produces
                   sequences with RNA nucleotides or non-DNA letters.

    Yields:
        - The generated sequence in a single line
    """

    # Valid DNA nucleotides
    nucleotides = "ATCG"
    # Valid DNA nucleotides + RNA + non-valid letters
    pert_nucleotides = "ATCGUNBDEFGHIJKLMOPQRSVWXYZ"

    # Select between modes I (Valid DNA sequences) and II (Invalid DNA sequences)
    vocab = nucleotides if not perturb else pert_nucleotides

    for n in range(num_seq):
        # Generate pseudo-random number for each sequence's length (between 10 and 200)
        seq_length = random.randint(10, 200)

        yield "".join(random.choices(vocab, k=seq_length))


def make_content(num_seq: int, perturb: bool) -> str:
    """Use header and sequence generators to create header + sequence pairs and
    and convert the generated strings to FASTA format taking into account sequences
    wrapping over more than one lines.

    Args:
        - num_seq: Number of sequences to be generated
        - perturb: Mode selection for generating valid (perturb=False) or invalid (perturb=True) DNA sequences

    Returns:
        - Formatted FASTA content of generated headers + sequences as a string
    """

    # Initialize empty lists for storing generated and formatted headers + sequences
    fasta_records = []
    headers = []
    sequences = []

    # Initialize generators
    header = header_generator(num_seq)
    sequence = seq_generator(num_seq, perturb)

    # Generate content for a specified number of headers and sequences
    for n in range(num_seq):
        headers.append(next(header))
        sequences.append(next(sequence))

    for head, seq in zip(headers, sequences):
        # Initialize empty subsequences list for every new header + sequence record
        subseqs = []
        # Append header followed by newline
        fasta_records.append(head + "\n")

        # Append sequence followed by newline if not over 60 nucleotides long
        if len(seq) <= 60:
            fasta_records.append(seq + "\n")

        # Split sequence into subsequences of 60 letters max length if original length over 60 nucleotides
        else:
            for ind in range(0, len(seq), 60):
                subseqs.append(seq[ind : ind + 60])
            # Join them into string while adding newline characters to follow FASTA format
            fasta_records.append("\n".join(subseqs) + "\n")

    return "".join(fasta_records)


def make_file(file_name: str, fasta_content: str) -> None:
    """Use formatted FASTA content (headers+ sequences) dictionary
    and write it to user-specified output file name.

    Args:
        - file_name: Output file name (e.g. DNA.fna)
        - formatted_cont: FASTA-formatted headers & sequences in string type

    Returns:
        - None
    """

    # Open output file
    with open(file_name, "w") as f:
        f.write(fasta_content)

    # Print message for user after closing and saving file
    print(f"FASTA file generated! \nSaved as: {file_name}.")


if __name__ == "__main__":
    args = parse_args()
    fasta_content = make_content(args.num_seq, args.perturb)
    make_file(args.file_name, fasta_content)

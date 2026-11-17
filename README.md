# fastx-utils

A small, dependency-free Python toolkit for everyday DNA sequence tasks: parsing FASTA/FASTQ files, translating DNA to protein, counting amino acids, demultiplexing barcoded reads, and generating synthetic test data.

Every script is a standalone command-line tool, and they all share one I/O module (`file_io.py`), so parsing and validation logic lives in one place.

## Contents

| File               | Type   | What it does                                                                                |
| ------------------ | ------ | ------------------------------------------------------------------------------------------- |
| `file_io.py`       | module | Reads FASTA (multi-line aware) and FASTQ files into dictionaries; writes formatted output   |
| `helpers.py`       | module | Standard genetic code (NCBI Table 1) and the set of all 64 codons                           |
| `dna2aa.py`        | CLI    | Translates DNA sequences in a FASTA file into protein sequences                             |
| `amino_count.py`   | CLI    | Counts absolute amino acid abundance across a protein FASTA file                            |
| `barcode.py`       | CLI    | Demultiplexes a FASTQ file by exact barcode match at the 5' or 3' end and trims the barcode |
| `generate_data.py` | CLI    | Generates random FASTA files of valid (or deliberately invalid) DNA for testing             |
| `tests/`           | tests  | Driver scripts that exercise each tool                                                      |

## Requirements

- Python **3.13+** (`generate_data.py` uses the single-argument `Generator[str]` type hint introduced in 3.13; the other scripts run on 3.9+)
- Standard library only, nothing to install

## Quick start

```bash
git clone git@github.com:darwinsorchid/DNASeq_funcs.git
cd DNASeq_funcs

# 1. Make some test data
python generate_data.py DNA.fna 20

# 2. Translate it
python dna2aa.py DNA.fna proteins.faa

# 3. Count amino acids in the translation
python amino_count.py proteins.faa aa_counts.tsv
```

## Usage

### `generate_data.py`: synthetic FASTA

```bash
python generate_data.py <output.fna> [num_seq] [-p]
```

- `num_seq`: number of sequences (default `10`); each one is 10–200 nt long, with random length
- `-p, --perturb`: mix RNA bases and non-DNA letters into the sequences, to test a program's input QC

Headers are written as `>sequenceID-001`, `>sequenceID-002`, …

### `dna2aa.py`: DNA → protein

```bash
python dna2aa.py <input.fna> [output_file]   # default output: output_file.txt
```

- Translates with the standard genetic code, starting at the **first base** (frame +1 only, no start-codon search)
- Stop codons are written as `*`, and translation **continues past them**
- Leftover 1–2 bases at the 3' end are ignored
- Input is upper-cased first; sequences containing any letter other than `ACGT` are reported on stdout and skipped
- Output is FASTA, one line per sequence

### `amino_count.py`: amino acid abundance

```bash
python amino_count.py <input.faa> [output_file]   # default output: output_file.txt
```

Writes a tab-separated table of the 20 standard amino acids and their absolute counts across all sequences. Any other character, including the stop symbol `*`, is counted under `X`.

```
A	12
C	3
...
Y	7
X	4
```

### `barcode.py`: FASTQ demultiplexing

```bash
python barcode.py <input.fastq> [-o TITLE]   # default TITLE: sample
```

- Looks for an **exact** match at either the 5' or 3' end of each read
- Default barcodes: `TATCCTCT`, `GTAAGGAG`, `TCTCTCCG` (all barcodes must be the same length)
- On a match, the barcode is trimmed from both the sequence and the quality line
- Output files are `TITLE1.fastq`, `TITLE2.fastq`, … (one per barcode, in order), plus `undetermined.fastq`, which holds unmatched reads untrimmed

To use your own barcodes from Python:

```python
from file_io import read_fastq
from barcode import barcode_filter, write_files

samples = barcode_filter(read_fastq("reads.fastq"), barcode_list=["ACGTACGT", "TTGGCCAA"])
write_files(samples, out_name="run42_")
```

## Using the modules in your own code

```python
from file_io import read_file, read_fastq
from dna2aa import translate_seq

records = read_file("DNA.fna")              # {">header": "ACGT..."}
proteins = {h: translate_seq(s) for h, s in records.items()}

reads = read_fastq("reads.fastq")           # {"@header": ("ACGT...", "IIII...")}
```

`read_file` raises `ValueError` for duplicate headers, headers with no sequence, files with no FASTA header, and empty files. `read_fastq` checks that records have 4 lines, that headers start with `@`, that there is a `+` separator, that sequence and quality lengths match, and that there are no duplicate headers.

## Running the tests

Run the tests as modules from the repository root so the imports resolve:

```bash
python -m tests.test_file_io
python -m tests.test_dna2aa
python -m tests.test_amino_count
python -m tests.test_barcode
```

The tests expect input files in `tests/test_data/` (for example `amino.faa`, `barcode.fastq`, `empty_file.txt`). That folder is not tracked in git, so create it locally first. You can generate FASTA inputs with `generate_data.py`. The tests print ✅/❌ messages rather than using a test framework.

## Known limitations

- Only the standard genetic code is supported (other tables can be added in `helpers.py`)
- Translation uses one reading frame: no start-codon search and no 6-frame translation
- Barcode matching is exact only, with no mismatch tolerance, and barcodes are only searched at the read ends
- `generate_data.py -p` perturbs every sequence rather than a subset
- Output overwrites existing files with the same name without warning

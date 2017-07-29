#!/usr/bin/env python3

"""
usage: random_lines.py [-h] file n

Retrieves a random sample of lines from a file.

positional arguments:
  file        Input file. If FASTA, will retrieve random header/sequence
              pairs.
  n           The size of the sample to retrieve.

optional arguments:
  -h, --help  show this help message and exit

"""
import sys
import argparse
import random
from collections import deque
import mmap

def fasta_parse(fasta, delimiter=">", separator="", trim_header=True):
    """
    Iterator which takes FASTA as input. Yields
    header/value pairs. Separator will be
    used to join the return value; use separator=
    None to return a list.

    If trim_header, parser will return the
    FASTA header up to the first space character.
    Otherwise, it will return the full, unaltered
    header string.

    """
    header, seq = None, []
    with open(fasta) as f:
        for line in f:
            if line.startswith(delimiter):
                if header is not None:
                    if separator is not None:
                        seq = separator.join(str(e) for e in seq)
                    yield header, seq
                header = line.strip().lstrip(delimiter)
                if trim_header:
                    try:
                        header = header.split()[0]
                    except IndexError:  # blank header
                        header = ""
                seq = []
            elif line.startswith('#'):
                continue
            else:
                if line.strip():  # don't collect blank lines
                    seq.append(line.rstrip('\n'))
        if separator is not None:
            seq = separator.join(str(e) for e in seq)
        yield header, seq


def check_fasta(f):
    """
    Checks for hallmark features of a FASTA file and returns a boolean.
    
    """
    with open(f) as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            elif line.startswith(">"):
                return True
            else:
                return False


# this version is slightly faster (memory-mapped file object optimized
# for line reading)
def count_lines(filename):
    fa = False
    if check_fasta(filename):
        fa = True
    with open(filename) as f:
        buf = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
        line_count = 0
        while True:
            line = buf.readline().decode()
            if not line:
                break
            if fa:
                if line.startswith('>'):
                    line_count += 1
            else:
                line_count += 1
    return line_count


def line_iterator(f):
    """
    Yield lines of a file, and if FASTA yields newline-joined
    header/sequence pairs
    
    """
    # Check if file is FASTA to determine correct parsing behavior
    is_fa = check_fasta(f)
    if is_fa:
        for h, s in fasta_parse(f):
            yield ">{}\n{}".format(h, s)
    else:
        with open(f) as infile:
            for l in infile:
                yield l.strip()


parser = argparse.ArgumentParser(description=(
    "Retrieves a random sample of lines from a file."))
parser.add_argument(
    'input_file',
    metavar='file',
    help='Input file. If FASTA, will retrieve random header/sequence pairs.')
parser.add_argument(
    'sample_size',
    metavar='n',
    type=int,
    help='The size of the sample to retrieve.')

args = parser.parse_args()

INFILE = args.input_file
SAMPLE_N = args.sample_size

# Get a count of all the lines (or header) in the file
TOTAL_LINES = count_lines(INFILE)

# deque object allows fast trimming of list-like data
RDM_SET = deque(sorted(random.sample(range(TOTAL_LINES), SAMPLE_N)))

# Iterate over lines or header/sequence pairs and print those in the random
# sample set
record_count = 0
for line in line_iterator(INFILE):
    try:
        if record_count == RDM_SET[0]:
            RDM_SET.popleft()
            print(line)
        record_count += 1
    except IndexError:
        sys.exit(0)

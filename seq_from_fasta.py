#!/usr/bin/env python3

"""
usage: seq_from_fasta.py [-h] [-f INFO_FILE] [--flank FLANK]
                         FASTA_file
                         [sequence_information [sequence_information ...]]

Retrieves sequences from a FASTA file using the following format: header
strand start stop [label]. If >label< is provided, will output in FASTA format
using >label< as header. Otherwise, outputs one sequence per line. If only
>header< is provided, will return the entire header/sequence pair.

positional arguments:
  FASTA_file            FASTA-formatted text file
  sequence_information  header strand start stop [label]

optional arguments:
  -h, --help            show this help message and exit
  -f INFO_FILE, --info_file INFO_FILE
                        File with information for sequences on separate lines.
  --flank FLANK         Size of any desired flanking region around specified
                        sequence(s)

"""
import sys
import argparse

def fasta_parse(fasta, delimiter=">", separator="", trim_header=False):
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

def demarcate(string, left_margin, right_margin, separator='|'):
    """
    Adds "|" characters to both ends of an input
    string, defining margins of length <margin>.

    """
    stringList = list(string)
    margins = [left_margin, -right_margin]
    for i in margins:
        if i != 0:
            stringList.insert(i, separator)
    return "".join(stringList)

def reverse_complement(seq):
    rev_list = []
    base_map = {
        'A':'T',
        'T':'A',
        'C':'G',
        'G':'C',
        'N':'N',
        'U':'T'
    }
    for c in seq:
        rev_list.append(base_map[c.upper()])
    rev_seq = ''.join(reversed(rev_list))

    return rev_seq

def get_subseq(parent_seq, start, stop, flank=0):
    """
    Returns a subsequence from <<parent_seq>>, adjusting
    indexing to 1-based (e.g. start=1, stop=5 --> parent_seq[0:5])

    """
    start = start - 1
    max_left_flank = len(parent_seq[:start])
    max_right_flank = len(parent_seq[stop:])
    l_flank = r_flank = flank
    if max_left_flank <= flank:
        l_flank = max_left_flank
    if max_right_flank <= flank:
        r_flank = max_right_flank
    out_string = parent_seq[start-l_flank:stop+r_flank]
    return out_string, l_flank, r_flank

def line_from_args(args):
    """
    Formats cmd-line arguments to match file line
    formatting

    """
    return "\t".join(args)

def line_to_info(line):
    bits = line.strip().split()
    loc, strand, start, stop = bits[:4]
    start, stop = min(int(start), int(stop)), max(int(start), int(stop))
    if len(bits) > 4:
        label = bits[-1]
    else:
        label = None
    line_info = {"loc": loc,
                 "strand": strand,
                 "start": int(start),
                 "stop": int(stop),
                 "label": label}
    return line_info

def seq_from_line(parent_seq, line):
    info = line_to_info(line)
    seq, l_flank, r_flank = get_subseq(
        parent_seq, info["start"], info["stop"], flank=FLANK)
    if info["strand"] == "-":
        seq = reverse_complement(seq)
    if FLANK:
        seq = demarcate(seq, l_flank, r_flank)
    return seq

def info_from_file(list_file):
    records = {}
    with open(list_file) as f:
        for l in f:
            l = l.strip()
            loc = l.split()[0]
            if loc not in records:
                records[loc] = []
            records[loc].append(l)
    return records

def info_from_args(args):
    loc = args[0]
    return {loc: [line_from_args(args)]}

def seqs_from_fasta(fasta, line_dict):
    num_keys = len(line_dict.keys())
    for h, s in fasta_parse(fasta, separator=SEP_CHAR):
        try:
            seq_lines = line_dict[h]
        except KeyError:
            continue
        num_keys -= 1
        if seq_lines[0] == h:  # they want whole seq
            yield s, {'label': h}
            if num_keys == 0:
                break
        else:
            for line in seq_lines:
                seq = seq_from_line(s, line)
                yield seq, line_to_info(line)
            if num_keys == 0:
                break

parser = argparse.ArgumentParser(
    description=(
        'Retrieves sequences from a FASTA file using the following '
        'format: header strand start stop [label]. If >label< is provided, '
        'will output in FASTA format using >label< as header. Otherwise, '
        'outputs one sequence per line. If only >header< is provided, will '
        'return the entire header/sequence pair.')
)
parser.add_argument(
    'FASTA_file',
    help='FASTA-formatted text file'
)
parser.add_argument(
    'sequence_information',
    nargs='*',
    help='header strand start stop [label]'
)
parser.add_argument(
    '-f',
    '--info_file',
    help='file with information for sequences on separate lines.'
)
parser.add_argument(
    '--flank',
    type=int,
    help='size of any desired flanking region around specified sequence(s)',
    default=0
)
parser.add_argument(
    '--unformatted',
    '-u',
    action='store_true',
    help=("leave original file formatting intact (e.g. don't join "
          "multi-line entries into a single line)")
)

if len(sys.argv) == 1:
    sys.exit(parser.print_help())

args = parser.parse_args()

FASTA = args.FASTA_file
FLANK = args.flank

if args.unformatted:
    SEP_CHAR = '\n'
else:
    SEP_CHAR = ''

# Parse information differently depending on source
if args.info_file:
    list_file = args.info_file
    info_dict = info_from_file(list_file)
else:
    direct_args = args.sequence_information
    info_dict = info_from_args(direct_args)

for seq, seq_info in seqs_from_fasta(FASTA, info_dict):
    if seq_info["label"]:
        print(">{}".format(seq_info["label"]))
    print(seq)

sys.exit(0)

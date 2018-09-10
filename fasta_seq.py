#!/usr/bin/python3

"""
usage: fasta_seq.py [-h] [-f INFO_FILE] [--flank FLANK] [-s SEPARATOR] [-u]
                    [-e] [--full_header]
                    FASTA_file [header strand start stop [label]
                    [header strand start stop [label] ...]]

Retrieves sequences from a FASTA file using the following format: header
strand start stop [label]. If >label< is provided, will output in FASTA format
using >label< as header. Otherwise, outputs one sequence per line. If only
>header< is provided, will return the entire header/sequence pair.

positional arguments:
  FASTA_file            FASTA-formatted text file
  header strand start stop [label]
                        header strand start stop [label]

optional arguments:
  -h, --help            show this help message and exit
  -f INFO_FILE, --info_file INFO_FILE
                        file with information for sequences on separate lines.
  --flank FLANK         size of any desired flanking region around specified
                        sequence(s)
  -s SEPARATOR, --separator SEPARATOR
                        Character to use to separate flanking sequence (if
                        any) from main sequence, which defaults to {tab}.
  -u, --unformatted     leave original file formatting intact (e.g. don't join
                        multi-line entries into a single line)
  -e, --exclude_header  exclude the header line from the output (does not
                        apply to explicitly labeled queries)
  --full_header         match on full header string (including whitespace)


"""
import sys
import argparse
from codecs import decode
from biogl import fasta_parse, flex_open, rev_comp


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
        seq = rev_comp(seq)
    if FLANK:
        seq = demarcate(seq, l_flank, r_flank, separator=FLANK_SEPARATOR)
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


def seqs_from_fasta(fasta, line_dict, exclude_header=False):
    num_keys = len(line_dict.keys())
    for h, s in fasta_parse(fasta, separator=SEP_CHAR, trim_header=False):
        complete_header = h
        if TRIM:
            h = h.split()[0]
        try:
            seq_lines = line_dict[h]
        except KeyError:
            continue
        num_keys -= 1
        if seq_lines[0] == h:  # they want whole seq
            if exclude_header:
                label = None
            else:
                label = complete_header
            yield s, {'label': label}
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
    metavar='header strand start stop [label]',
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
    '-s',
    '--separator',
    type=str,
    help='Character to use to separate flanking sequence (if any) from main '
    'sequence, which defaults to {tab}.',
    default='\t'
)
parser.add_argument(
    '-u',
    '--unformatted',
    action='store_true',
    help=('leave original file formatting intact (e.g. don\'t join '
          'multi-line entries into a single line)')
)
parser.add_argument(
    '-e',
    '--exclude_header',
    action='store_true',
    help='exclude the header line from the output (does not apply to '
    'explicitly labeled queries)'
)
parser.add_argument(
    '--full_header',
    action='store_true',
    help='match on full header string (including whitespace)'
)

args = parser.parse_args()

FASTA = args.FASTA_file
FLANK = args.flank
FLANK_SEPARATOR = decode(args.separator, 'unicode escape')
EXCLUDE_HEADER = args.exclude_header

if args.full_header:
    TRIM = False
else:
    TRIM = True

if args.unformatted:
    SEP_CHAR = '\n'
else:
    SEP_CHAR = ''

# Parse information differently depending on source
if args.info_file:
    list_file = args.info_file
    info_dict = info_from_file(list_file)
else:
    info_length = len(args.sequence_information)
    if 2 <= info_length < 4:
        sys.exit('Insufficient sequence information provided. Exiting.')
    direct_args = args.sequence_information
    info_dict = info_from_args(direct_args)

FASTA = flex_open(FASTA)

for seq, seq_info in seqs_from_fasta(FASTA, info_dict, EXCLUDE_HEADER):
    if seq_info["label"]:
        print(">{}".format(seq_info["label"]))
    print(seq)

sys.exit(0)

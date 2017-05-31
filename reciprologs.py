#!/usr/bin/env python3

"""
usage: reciprologs.py [-h] [-t THREADS] [-p [PARALLEL_PROCESSES]]
                      query subject {blastn,blastp,blastx,tblastn,tblastx}

Uses BLAST to find reciprocal best hits between two files.

positional arguments:
  query                 query file to be BLASTed
  subject               subject file to be BLASTed
  {blastn,blastp,blastx,tblastn,tblastx}
                        type of BLAST program to run

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        number of CPU threads to use (overridden by -p)
                        (default: 4)
  -p [PARALLEL_PROCESSES], --parallel_processes [PARALLEL_PROCESSES]
                        run the BLAST step using multiple parallel processes;
                        without specific input will use half available system
                        CPUs (default: None)

Depends on blast.py

"""
import sys
import subprocess
import os
import time
from operator import itemgetter
import argparse
from multiprocessing import cpu_count


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


def get_runtime(start_time, p=3):
    """
    Takes a start time and optional decimal precision p,
    returns a string of the total run-time until current
    time with appropriate units.

    """
    total = time.time() - start_time  # start with seconds
    divided = total/60.0
    if divided < 2:
        run_time = total
        units = "seconds"
    elif divided < 60:
        run_time = divided
        units = "minutes"
    else:
        run_time = divided/60.0
        units = "hours"
    rounded = round(run_time, p)
    return "{} {}".format(rounded, units)


def parse_blast_line(bl, *args):
    """
    Returns info from certain columns in a tab-separated BLAST
    output file. $args may be: query, subject, length, e, bitscore

    """
    columns = bl.strip().split("\t")
    query, subject, length, e_value, bitscore = itemgetter(0, 1, 4, 10, 11)(columns)
    arg_map = {
        "query": query,
        "subject": subject,
        "length": length,
        "e": float(e_value),
        "bitscore": float(bitscore)
    }
    results = []
    for a in args:
        results.append(arg_map[a])
    if len(results) == 1:
        results = results[0]
    return results


def get_top_hits(blast, paralogs=False):
    results = {}
    with open(blast) as blst:
        for l in blst:
            if l.startswith("#"):
                continue
            q, s, score = parse_blast_line(l, "query", "subject", "bitscore")
            # do not consider hits to self if BLASTing against self
            if paralogs and q == s:
                continue
            if q in results:
                # Check if this hit's score is better
                defender = results[q]["score"]
                if score > defender:
                    results[q] = {"name": s, "score": score}
            else:
                results[q] = {"name": s, "score": score}

    return results


def get_reciprocals(d1, d2):
    """
    Takes two dictionaries of top BLAST hits,
    returns a list of tuples of all pairs that were
    reciprocal best hits.

    """
    reciprologs = set()
    blast_directions = [d1, d2]
    for first, second in [(d1, d2), (d2, d1)]:
        for query, hit_info in first.items():
            best_hit = hit_info["name"]
            if best_hit in second:
                reciprocal_hit = second[best_hit]["name"]
                if query == reciprocal_hit:  # best hit refers back to query
                    hit_pair = sorted([query, best_hit])
                    reciprologs.add(tuple(hit_pair))
    return sorted(reciprologs)


def abbreviate(name, delimiter="."):
    name = name.split("/")[-1]  # in case of non-local file path
    abbreviation = name.split(delimiter)[0]
    return abbreviation


def concatenate(outname, file_list, clean=True):
    with open(outname, 'w') as outfile:
        for fn in file_list:
            with open(fn) as f:
                for l in f:
                    outfile.write(l)
    if clean:
        [os.remove(fn) for fn in file_list]
    

parser = argparse.ArgumentParser(
    description='Uses BLAST to find reciprocal best hits between two files.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    'query',
    help='query file to be BLASTed'
)
parser.add_argument(
    'subject',
    help='subject file to be BLASTed'
)
parser.add_argument(
    'blast_type',
    choices=[
        'blastn',
        'blastp',
        'blastx',
        'tblastn',
        'tblastx'],
    help='type of BLAST program to run'
)
parser.add_argument(
    '-t',
    '--threads',
    type=int,
    help='number of CPU threads to use (overridden by -p)',
    default=4
)
parser.add_argument(
    '-p',
    '--parallel_processes',
    help=('run the BLAST step using multiple parallel processes; '
          'without specific input will use half available system CPUs'),
    type=int,
    const=round(cpu_count() / 2),
    nargs='?'
)

if len(sys.argv) == 1:
    sys.exit(parser.print_help())

t_start = time.time()

args = parser.parse_args()

BLAST_TYPE = args.blast_type
THREADS = args.threads
PARALLEL = args.parallel_processes
SUBJECT = args.subject
QUERY = args.query

BLAST = 'blast.py'

# special case of BLASTing against self
if QUERY == SUBJECT:
    PARALOGS = True
else:
    PARALOGS = False

# BLAST in both directions (unless PARALOGS)

fw_names = [abbreviate(QUERY), abbreviate(SUBJECT)]
rv_names = fw_names[::-1]

fw_fn = '{}-{}.{}.blast_results'.format(*fw_names, BLAST_TYPE)
rv_fn = '{}-{}.{}.blast_results'.format(*rv_names, BLAST_TYPE)

# create a list with the flags/options to pass to the subsequence
# subprocess calls
call_options = {
    '-p': PARALLEL,
    '-t': THREADS,
}
optional = []
for k, v in call_options.items():
    if v:
        optional.extend([str(k), str(v)])

forward_args = [BLAST, SUBJECT, QUERY, BLAST_TYPE, '-n', fw_fn] + optional
subprocess.run(forward_args)
top_forward = get_top_hits(fw_fn, PARALOGS)
if PARALOGS:
    # filter for reciprocal best hits within the same file
    top_reverse = top_forward
    # reciprologs = get_reciprocals(top_forward, top_forward)
else:
    reverse_args = [BLAST, QUERY, SUBJECT, BLAST_TYPE, '-n', rv_fn] + optional
    subprocess.run(reverse_args)
    top_reverse = get_top_hits(rv_fn)

# Filter for reciprocal best hits
reciprologs = get_reciprocals(top_forward, top_reverse)

out_file = "{}-{}.{}.reciprologs".format(
    abbreviate(QUERY), abbreviate(SUBJECT), BLAST_TYPE)
with open(out_file, 'w') as out:
    for pair in reciprologs:
        out.write("\t".join(pair) + "\n")

runtime = get_runtime(t_start)

print('[#] Job finished in {}; {} pairs found. '
      'Results are in \'{}\''.format(runtime, len(reciprologs), out_file),
      file=sys.stderr)

sys.exit(0)

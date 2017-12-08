#!/usr/bin/env python3

"""
usage: reciprologs.py [-h] [-t THREADS] [-p [PARALLEL_PROCESSES]]
                      file_1 file_2 {blastn,blastp,blastx,tblastn,tblastx}

Uses BLAST to find reciprocal best hits between two files.

positional arguments:
  file_1                file to BLAST
  file_2                file to BLAST
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

NOTE: Depends on pblast.py

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
    (query, subject, length, 
    e_value, bitscore, q_start,
    q_stop) = itemgetter(0, 1, 3, 10, 11, 6, 7)(columns)
    arg_map = {
        "query": query,
        "subject": subject,
        "length": int(length) - 1,  # seem to be off by 1 in BLAST output
        "e": float(e_value),
        "bitscore": float(bitscore)
    }
    results = []
    for a in args:
        results.append(arg_map[a])
    if len(results) == 1:
        results = results[0]
    return results


def get_top_hits(blast, paralogs=False, query_match=None):
    results = {}
    with open(blast) as blst:
        for l in blst:
            if l.startswith("#"):
                continue
            (q, s, score, length) = parse_blast_line(
                l, "query", "subject", "bitscore", "length")
            # do not consider hits to self if BLASTing against self
            if paralogs and q == s:
                continue
            if query_match:
                # use query_match dictionary to compare query lengths to
                # match lengths to exclude matches where query percentage 
                # is below query_match_threshold key
                fraction = (length / query_match[q]) * 100
                if fraction < query_match['query_match_threshold']:
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
    'file_1',
    help='file to BLAST',
)
parser.add_argument(
    'file_2',
    help='file to BLAST'
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
    '-p',
    '--parallel_processes',
    help=(
        'run the BLAST step using multiple parallel processes; '
        'without specific input will use half available system CPUs'),
    type=int,
    const=round(cpu_count() / 2),
    nargs='?'
)
parser.add_argument(
    '-t',
    '--query_percentage_threshold',
    help=(
        'require a specified fraction of the query length to match in '
        'order for a hit to qualify (lowest allowable percentage'),
    type=float,
    default=None
)
parser.add_argument(
    '--overwrite',
    help='overwrite existing BLAST files (instead of using to bypass BLAST)',
    action='store_true'
)

if len(sys.argv) == 1:
    sys.exit(parser.print_help())

t_start = time.time()

args, EXTRA_ARGS = parser.parse_known_args()

BLAST_TYPE = args.blast_type
PARALLEL = args.parallel_processes
QUERY = args.file_1
SUBJECT = args.file_2
QUERY_PERCENTAGE = args.query_percentage_threshold
OVERWRITE = args.overwrite

BLAST = 'pblast.py'

# special case of BLASTing against self
if QUERY == SUBJECT:
    PARALOGS = True
else:
    PARALOGS = False

q_lengths, s_lengths = {}, {}
if QUERY_PERCENTAGE:
    # we need to get sequence lengths for each file
    for fa, target in zip([QUERY, SUBJECT], [q_lengths, s_lengths]):
        target['query_match_threshold'] = QUERY_PERCENTAGE
        for h, s in fasta_parse(fa):
            target[h] = len(s)


# BLAST in both directions (unless PARALOGS)

fw_names = [abbreviate(QUERY), abbreviate(SUBJECT)]
rv_names = fw_names[::-1]

fw_fn = '{}-{}.{}.blast_results'.format(*fw_names, BLAST_TYPE)
rv_fn = '{}-{}.{}.blast_results'.format(*rv_names, BLAST_TYPE)

# create a list with the flags/options to pass to the subsequence
# subprocess calls
call_options = {
    '-p': PARALLEL
}
optional = EXTRA_ARGS
for k, v in call_options.items():
    if v:
        optional.extend([str(k), str(v)])

# use existing BLAST output unless --overwrite is specified
if not os.path.isfile(fw_fn) or OVERWRITE:
    forward_args = [BLAST, QUERY, SUBJECT, BLAST_TYPE, '-n', fw_fn] + optional
    subprocess.run(forward_args)
else:
    print('[#] Using existing BLAST output \'{}\''.format(fw_fn))
top_forward = get_top_hits(
    fw_fn, PARALOGS, query_match=q_lengths)
if PARALOGS:
    # don't need to run BLAST again; process the same file
    top_reverse = get_top_hits(
        fw_fn, PARALOGS, query_match=s_lengths)
else:
    if not os.path.isfile(rv_fn) or OVERWRITE:
        reverse_args = [BLAST, SUBJECT, QUERY, BLAST_TYPE, '-n', rv_fn] + optional
        subprocess.run(reverse_args)    
    else:
        print('[#] Using existing BLAST output \'{}\''.format(rv_fn))
    top_reverse = get_top_hits(rv_fn, query_match=s_lengths)

# Filter for reciprocal best hits
reciprologs = get_reciprocals(top_forward, top_reverse)

out_file = "{}-{}.{}.reciprologs".format(
    abbreviate(QUERY), abbreviate(SUBJECT), BLAST_TYPE)
with open(out_file, 'w') as out:
    for pair in reciprologs:
        out.write("\t".join(pair) + "\n")

runtime = get_runtime(t_start)

print(
    '[#] Job finished in {}; {} pairs found: '
    '\'{}\''.format(runtime, len(reciprologs), out_file),
      file=sys.stderr)

sys.exit(0)

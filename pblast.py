#!/usr/bin/env python3

"""
usage: pblast.py [-h] [-t THREADS] [-p [PARALLEL_PROCESSES]] [-n NAME]
                query subject {blastn,blastp,blastx,tblastn,tblastx}

BLAST one file against another

positional arguments:
  query                 query file to be BLASTed
  subject               subject file to be BLASTed
  {blastn,blastp,blastx,tblastn,tblastx}
                        type of BLAST program to run

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        number of CPU threads to use (overridden if -p)
                        (default: 4)
  -p [PARALLEL_PROCESSES], --parallel_processes [PARALLEL_PROCESSES]
                        run the BLAST step using multiple parallel processes;
                        without specific input will use half available system
                        CPUs (default: None)
  -n NAME, --name NAME  filename for results (otherwise, automatic) (default:
                        None)

"""
import sys
import subprocess
import os
import time
import argparse
from multiprocessing import Pool, cpu_count
from functools import partial
from math import ceil

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


def abbreviate(name, delimiter="."):
    name = name.split("/")[-1]  # in case of non-local file path
    abbreviation = name.split(delimiter)[0]
    return abbreviation


def make_blast_db(fasta, db_type="nucleotide"):
    """
    Creates a local BLAST+ database using a
    FASTA file. $db_type may be either "protein"
    or "nucleotide".

    """
    print("[#] Creating local BLAST database", file=sys.stderr)
    type_map = {
        'protein': 'prot',
        'nucleotide': 'nucl'
    }
    cmd_args = [
        "makeblastdb",
        "-in",
        fasta,
        "-parse_seqids",
        "-dbtype",
        type_map[db_type]
    ]
    subprocess.call(cmd_args)
    print("\n[#] BLAST database for '{}' created"
          .format(fasta), file=sys.stderr)


def local_blast(
        db_file, 
        blast_type, 
        query_file, 
        filename,
        out_fmt=6, 
        threads=1, 
        e_value=1e-10):
    """
    Runs BLAST+ sub-program (blastp, blastn, blastx, tblastn, tblastx)
    with the given query on the given database.

    Optional: out_fmt = type as per BLAST+ documentation; threads = number
    of threads to use for job; e = e-value cutoff for hits

    """
    start_time = time.time()
    db_abbrev = abbreviate(db_file)
    query_abbrev = abbreviate(query_file)
    if filename is None:
        filename = "{}-{}.blast_results".format(query_abbrev, db_abbrev)
    cmd_args = [
        blast_type,
        "-db",
        db_file,
        "-query",
        query_file,
        "-evalue",
        e_value,
        "-outfmt",
        out_fmt,
        "-out",
        filename,
        "-num_threads",
        threads,
    ]
    cmd_args = [str(c) for c in cmd_args]
    cmd_string = ' '.join(cmd_args)
    print("[#] Starting BLAST search on '{}' vs. '{}':\n{}"
          .format(query_file, db_file, cmd_string), file=sys.stderr)
    subprocess.call(cmd_args)
    run_time = get_runtime(start_time)
    print("[#] BLAST finished in {}".format(run_time), file=sys.stderr)
    return filename


def is_fasta(some_file):
    with open(some_file) as f:
        for l in f:
            if l.startswith("#"):
                continue
            if l.startswith(">"):
                return True
            else:
                return False


def seq_type(fasta):
    # just get first sequence
    _, s = next(fasta_parse(fasta))
    if any(c not in 'ATCGN' for c in s):
        return 'protein'
    else:
        return 'nucleotide'


def make_fasta(some_file):
    fasta_name = os.path.basename(some_file) + ".fasta"
    if not os.path.isfile(fasta_name):
        with open(some_file) as infile, open(fasta_name, 'w') as outfile:
            for i, l in enumerate(infile, start=1):
                outfile.write('>{}\n{}\n'.format(i, l.strip()))
    return fasta_name


def make_pro(fasta, ftype=''):
    fasta_name = os.path.basename(fasta)
    out_name = "{}.pro".format(fasta_name)
    if ftype:
        ftype = ' for {}'.format(ftype)
    if os.path.isfile(out_name):
        print('[#] Using exising file \'{}\'{}'
              .format(out_name, ftype), file=sys.stderr)
    else:
        print("[#] Starting translation of '{}'{}"
              .format(fasta, ftype), file=sys.stderr)
        subprocess.call(["translator.py", fasta],
                        stdout=open(out_name, 'w'))
        print("[#] Translation finished", file=sys.stderr)
    return out_name


def prep_blast(subject, query, blast_type):
    db_type_map = {
        'blastn': {
            'query': 'nucleotide',
            'subject': 'nucleotide'},
        'blastp': {
            'query': 'protein',
            'subject': 'protein'},
        'blastx': {
            'query': 'nucleotide',
            'subject': 'protein'},
        'tblastn': {
            'query': 'protein',
            'subject': 'nucleotide'},
        'tblastx': {
            'query': 'nucleotide',
            'subject': 'nucleotide'}
    }

    blast_files = {'query': query, 'subject': subject}

    # determine what format the files should be (pro or nucl)
    # based upon the BLAST type
    target_format = {}
    for ftype, fn in blast_files.items():
        target_format[ftype] = db_type_map[blast_type][ftype]

    # convert the appropriate files to protein if needed
    for ftype, fn in blast_files.items():
        fmt = target_format[ftype]
        if fmt == 'protein' and seq_type(fn) != 'protein':
            fn = make_pro(fn, ftype)
            # add new fn to dict
            target_format[ftype] = fmt
            blast_files[ftype] = fn
    
    subject = blast_files['subject']
    query = blast_files['query']

    # make BLAST database
    make_blast_db(subject, db_type=target_format['subject'])

    return subject, query


def parallel_blast(
        subject, 
        query, 
        blast_type, 
        e_value=1e-10, 
        out_name=None):
    pool = Pool(PARALLEL)
    blast = partial(local_blast, subject, blast_type, e_value=e_value)
    count = 0
    for h, s in fasta_parse(query):
        count += 1
    block_size = ceil(count / PARALLEL)
    chunk_name = '{}.{}.chunk'.format(query, 1)
    chunk = open(chunk_name, 'w')
    chunked = [chunk_name]
    tally = 1
    for index, (h, s) in enumerate(fasta_parse(query), start=1):   
        if index % block_size == 0:
            tally += 1
            # order matters here
            chunk.close()
            chunk_name = '{}.{}.chunk'.format(query, tally)
            chunked.append(chunk_name)
            chunk = open(chunk_name, 'w')
            chunk.write('>{}\n{}\n'.format(h, s))
        else:
            chunk.write('>{}\n{}\n'.format(h, s))
    chunk.close()

    chunk_list = '\n'.join(chunked)    
    print('[#] Query split into {} files:\n{}\n'
          .format(len(chunked), chunk_list))
    if out_name is None:
        out_name = '{}-{}.{}.blast_results'.format(
            abbreviate(query), abbreviate(subject), blast_type)
    filenames = [
        '{}.{}.tmp'.format(out_name, i) 
        for i in range(1, len(chunked) + 1)
        ]
    results = pool.starmap(blast, zip(chunked, filenames))

    # clean up chunk files
    [os.remove(f) for f in chunked]

    return results


def concatenate(outname, file_list, clean=True):
    with open(outname, 'w') as outfile:
        for fn in file_list:
            with open(fn) as f:
                for l in f:
                    outfile.write(l)
    if clean:
        [os.remove(fn) for fn in file_list]


parser = argparse.ArgumentParser(
    description='BLAST one file against another',
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
    help='number of CPU threads to use (overridden if -p)',
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
parser.add_argument(
    '-n',
    '--name',
    type=str,
    help='filename for results (otherwise, automatic)'
)
parser.add_argument(
    '-e',
    '--e_value',
    help='e-value threshold to use for search',
    type=float,
    default=1e-10
)

if len(sys.argv) == 1:
    sys.exit(parser.print_help())

args = parser.parse_args()

BLAST_TYPE = args.blast_type
THREADS = args.threads
PARALLEL = args.parallel_processes
OUT_NAME = args.name
E_VALUE = args.e_value

run_files = {'query': args.query, 'subject': args.subject}

for ftype, fn in run_files.items():
    if not is_fasta(fn):
        run_files[ftype] = make_fasta(fn)

SUBJECT = run_files['subject']
QUERY = run_files['query']

if not OUT_NAME:
    OUT_NAME = '{}-{}.{}.blast_results'.format(
        abbreviate(QUERY), 
        abbreviate(SUBJECT),
        BLAST_TYPE)

SUBJECT, QUERY = prep_blast(SUBJECT, QUERY, BLAST_TYPE)

if PARALLEL:
    pblast_out = parallel_blast(
        SUBJECT, 
        QUERY, 
        BLAST_TYPE, 
        e_value=E_VALUE, 
        out_name=OUT_NAME)
    concatenate(OUT_NAME, pblast_out)
else:
    blast_out = local_blast(
        SUBJECT, 
        BLAST_TYPE, 
        QUERY, 
        filename=OUT_NAME, 
        threads=THREADS,
        e_value=E_VALUE)

sys.exit(0)

#!/usr/bin/env python3

"""
usage: pblast.py [-h] [-p [PARALLEL_PROCESSES]] [-s] [-o OUTPUT_FORMAT]
                 [-n NAME] [-e E_VALUE] [--clobber_db]
                 query subject {blastn,blastp,blastx,tblastn,tblastx}

BLAST one file against another. Any arguments not listed here will be passed
to BLAST unmodified.

positional arguments:
  query                 query file to be BLASTed
  subject               subject file to be BLASTed against
  {blastn,blastp,blastx,tblastn,tblastx}
                        type of BLAST program to run

optional arguments:
  -h, --help            show this help message and exit
  -p [PARALLEL_PROCESSES], --parallel_processes [PARALLEL_PROCESSES]
                        run the BLAST step using multiple parallel processes;
                        without specific input will use half of available
                        system CPUs
  -s, --single          disable parallel processing (default: False)
  -o OUTPUT_FORMAT, --output_format OUTPUT_FORMAT
                        integer output format for BLAST results (default: 6)
  -n NAME, --name NAME  filename for results (otherwise, automatic based on
                        input) (default: None)
  -e E_VALUE, --e_value E_VALUE
                        e-value threshold to use for search (default: 1e-10)
  --clobber_db          create new database even if one already exists
                        (default: False)
"""
import sys
import subprocess
import os
import time
import argparse
from multiprocessing import Pool, cpu_count
from functools import partial
from math import ceil, log

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
    name = os.path.basename(name)  # in case of non-local file path
    abbreviation = name.split(delimiter)[0]

    return abbreviation


def unique_filenames(*file_list):
    abbreviated = [abbreviate(f) for f in file_list]
    if len(set(abbreviated)) < len(abbreviated):  # not all are unique
        return [os.path.basename(f) for f in file_list]
    else:
        return abbreviated


def db_check(db_filename):
    db_directory = os.path.dirname(os.path.abspath(db_filename))
    db_dir_files = os.listdir(db_directory)
    db_endings = ('sq', 'si', 'sd', 'og', 'in', 'hr')
    db_files = [f.rsplit('.', 1) for f in db_dir_files]
    # get just the file endings
    db_name = os.path.basename(db_filename)    
    db_files = [f[1] for f in db_files if f[0] == db_name]
    present_endings = [f[-2:] for f in db_files]
    if all(dbe in present_endings for dbe in db_endings):
        previous_db = True
    else:
        previous_db = False
    
    return previous_db


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
        e_value=1e-10,
        extra_blast_args=None):
    """
    Runs BLAST+ sub-program (blastp, blastn, blastx, tblastn, tblastx)
    with the given query on the given database.

    Optional: out_fmt = type as per BLAST+ documentation; threads = number
    of threads to use for job; e = e-value cutoff for hits

    extra_blast_args is a list of additional arguments to pass to BLAST

    """
    if extra_blast_args is None:
        extra_blast_args = []  # make unpacking syntax work
    start_time = time.time()
    db_abbrev = abbreviate(db_file)
    query_abbrev = abbreviate(query_file)
    if filename is None:
        filename = "{}-vs-{}.{}.blast".format(
            query_abbrev, db_abbrev, blast_type)
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
        *extra_blast_args
    ]
    cmd_args = [str(c) for c in cmd_args]
    cmd_string = ' '.join(cmd_args)
    vs = "'{}' vs. '{}'".format(query_file, db_file)
    print("[#] Starting BLAST on {}:\n{}"
          .format(vs, cmd_string), file=sys.stderr)
    subprocess.call(cmd_args)
    run_time = get_runtime(start_time)
    print("[#] BLAST on {} finished in {}".format(vs, run_time), file=sys.stderr)
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


def prep_blast(subject, query, blast_type, overwrite=True):
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

    # check for already-created database
    if db_check(subject):
        if overwrite:
            make_blast_db(subject, db_type=target_format['subject'])
        else:
            print('[#] Using existing BLAST database for \'{}\''
            .format(subject), file=sys.stderr)
    else:
        make_blast_db(subject, db_type=target_format['subject'])

    return subject, query


def parallel_blast(
        subject, 
        query, 
        blast_type,
        out_fmt=6, 
        e_value=1e-10, 
        out_name=None,
        extra_blast_args=None):
    pool = Pool(PARALLEL)
    blast = partial(
        local_blast, 
        subject, 
        blast_type,
        out_fmt=out_fmt,
        e_value=e_value, 
        extra_blast_args=extra_blast_args)
    count = sum([1 for p in fasta_parse(query)])
    block_size = ceil(count / PARALLEL)
    # change padding depth according to total number of chunks
    zero_pad = ceil(log(PARALLEL + 1, 10))
    chunk_name_scheme = '{0}.{1:0{2}}.chunk'
    chunk_name = chunk_name_scheme.format(query, 1, zero_pad)
    chunk = open(chunk_name, 'w')
    chunked = [chunk_name]
    tally = 1
    for index, (h, s) in enumerate(fasta_parse(query), start=1):   
        if index > block_size and index % block_size == 1:
            tally += 1
            # order matters here
            chunk.close()
            chunk_name = chunk_name_scheme.format(query, tally, zero_pad)
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
        out_name = '{}-vs-{}.{}'.format(
            abbreviate(query), abbreviate(subject), blast_type)
    filenames = [
        '{0}.{1:0{2}}.tmp'.format(out_name, i, zero_pad) 
        for i in range(1, len(chunked) + 1)
        ]
    # use apply_async instead of starmap to allow messages
    # to print to screen without clobbering each other
    results = []
    for pair in zip(chunked, filenames):
        results.append(pool.apply_async(blast, args=pair))
        time.sleep(.01)
    pool.close()
    pool.join()
    results = [r.get() for r in results]

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
    description=(
        'BLAST one file against another. '
    'Any arguments not listed here will be passed to BLAST unmodified.'),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
    allow_abbrev=False)

parser.add_argument(
    'query',
    help='query file to be BLASTed'
)
parser.add_argument(
    'subject',
    help='subject file to be BLASTed against'
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
    help=('run the BLAST step using multiple parallel processes; '
          'without specific input will use half of available system CPUs'),
    type=int,
    const=round(cpu_count() / 2),
    nargs='?',
    default=round(cpu_count() / 2)
)
parser.add_argument(
    '-s',
    '--single',
    help='disable parallel processing',
    action='store_true'
)
parser.add_argument(
    '-o',
    '--output_format',
    type=int,
    help='integer output format for BLAST results',
    default=6
)
parser.add_argument(
    '-n',
    '--name',
    type=str,
    help='filename for results (otherwise, automatic based on input)'
)
parser.add_argument(
    '-e',
    '--e_value',
    help='e-value threshold to use for search',
    type=float,
    default=1e-10
)
parser.add_argument(
    '--clobber_db',
    help='create new database even if one already exists',
    action='store_true'
)

if len(sys.argv) == 1:
    sys.exit(parser.print_help())

args, EXTRA_ARGS = parser.parse_known_args()

BLAST_TYPE = args.blast_type
# THREADS = args.threads
PARALLEL = args.parallel_processes
SINGLE = args.single
OUT_NAME = args.name
E_VALUE = args.e_value
OUT_FORMAT = args.output_format
OVERWRITE = args.clobber_db

run_files = {'query': args.query, 'subject': args.subject}

for ftype, fn in run_files.items():
    if not is_fasta(fn):
        run_files[ftype] = make_fasta(fn)

SUBJECT = run_files['subject']
QUERY = run_files['query']

if not OUT_NAME:
    subj, quer = unique_filenames(SUBJECT, QUERY)
    OUT_NAME = '{}-vs-{}.{}'.format(
        quer, 
        subj,
        BLAST_TYPE)

SUBJECT, QUERY = prep_blast(SUBJECT, QUERY, BLAST_TYPE, overwrite=OVERWRITE)

if not SINGLE and PARALLEL > 1: #PARALLEL:
    # run multiple processes, and then join the output afterward
    pblast_out = parallel_blast(
        SUBJECT, 
        QUERY, 
        BLAST_TYPE,
        out_fmt=OUT_FORMAT,
        e_value=E_VALUE, 
        out_name=OUT_NAME,
        extra_blast_args=EXTRA_ARGS)
    concatenate(OUT_NAME, pblast_out)
else:
    blast_out = local_blast(
        SUBJECT, 
        BLAST_TYPE, 
        QUERY,
        out_fmt=OUT_FORMAT,
        filename=OUT_NAME, 
        e_value=E_VALUE,
        extra_blast_args=EXTRA_ARGS)

sys.exit(0)

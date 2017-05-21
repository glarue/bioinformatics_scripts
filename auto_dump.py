#!/usr/bin/env python3

"""
usage: auto_dump.py [-h] [-f FILE] [-k] [accessions [accessions ...]]

A program to run fastq-dump sequentially on a list of accession numbers. Will
automatically detect if reads are single- or paired-end and will run fastq-
dump accordingly, adding proper ID suffixes as needed. Accession numbers may
be provided directly, or in a file using the -f option. If provided directly,
accessions may be comma or space separated, or hyphen-separated to specify a
range, e.g. SRR123455, SRR123456, SRR123457 || SRR123455-SRR123457

positional arguments:
  accessions            List of SRA accession numbers (e.g. SRR123456
                        SRR123457)

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  File with SRA accession numbers on separate lines
  -k, --keep_sra_files  Keep SRA files after dumping to FASTQ format

"""
import sys
import os
import subprocess
from collections import namedtuple
from itertools import zip_longest
import shlex
import time
import argparse

# taken directly from https://docs.python.org/3/library/itertools.html#recipes
def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def fastq_parse(fastq):
    FqRec = namedtuple("FqRec", ["id", "seq", "meta", "quality", "full_record", "cleaned"])
    for group in grouper(fastq, 4):  # blocks of 4 lines
        if None in group:
            continue
        # assign fastq lines to a namedtuple output
        record_block = [g.strip() for g in group]
        full_record = "\n".join(record_block)
        fq_id, fq_seq, fq_meta, fq_quality = record_block
        fq_id = fq_id.split()[0]
        fq_meta = fq_meta[0]
        cleaned = "\n".join([fq_id, fq_seq, fq_meta, fq_quality])
        record_info = FqRec(fq_id, fq_seq, fq_meta, fq_quality, full_record, cleaned)
        yield record_info

def paired_check(acc):
    """
    Runs fastq-dump on the first 2 reads of an accession number
    <acc> to determine whether the data is single- or paired-end.

    Returns "paired", "single" or None.

    """
    test_dump = "fastq-dump {} -X 6 --split-files -I --stdout".format(acc)
    test_args = shlex.split(test_dump)
    # stderr of subprocess sent to /dev/null to silence fastq-dump progress output
    test_process = subprocess.Popen(test_args, stdout=subprocess.PIPE,
                                    stderr=subprocess.DEVNULL, universal_newlines=True)
    suffixes = set()

    for rec in fastq_parse(test_process.stdout):
        # it's safe to just take the last character because of the -I flag
        # in the fastq-dump call
        suffixes.add(rec.id[-1])
    if len(suffixes) == 2:
        read_type = "paired"
    elif len(suffixes) == 1:
        read_type = "single"
    else:
        read_type = None
    return read_type


def fetch_sra(acc):
    """
    Grabs an SRA file from NCBI.
    
    """
    alpha = acc[:3]
    first_six = acc[:6]
    sra_filename = '{}.sra'.format(acc)

    # build file path to SRA data
    base_path = 'ftp://ftp-trace.ncbi.nlm.nih.gov'

    # pattern of path:
    #  /sra/sra-instant/reads/ByRun/sra/{SRR|ERR|DRR}/<first 6 characters of accession>/<accession>/<accession>.sra
    path_bits = [alpha, first_six, acc, sra_filename]
    relative_path = '/'.join(path_bits)
    sra_path = '/sra/sra-instant/reads/ByRun/sra/{}'.format(relative_path)
    full_path = base_path + sra_path

    subprocess.run(['wget', full_path, '-O', sra_filename])

    return sra_filename

def acc_list_from_file(acc_file):
    accs = []
    with open(acc_file) as f:
        for l in f:
            accs.append(l.strip())
    return accs

def clean_args(arguments):
    stage_1 = []
    for arg in arguments:
        new_arg = arg.replace(',', ' ')
        cleaned = new_arg.split()
        stage_1 += cleaned
    final_args = []
    for c in stage_1:
        if '-' not in c:
            final_args.append(c)
            continue
        begin, end = c.split('-')
        prefix = begin[:3]
        start, stop = list(map(int, [e[3:] for e in [begin, end]]))
        numeric_range = list(range(start, stop+1))
        acc_range = [prefix+str(r) for r in numeric_range]
        final_args += acc_range

    return final_args

parser = argparse.ArgumentParser(
    description='A program to run fastq-dump sequentially on a list of '
    'accession numbers. Will automatically detect if reads are single- or '
    'paired-end and will run fastq-dump accordingly, adding proper ID '
    'suffixes as needed. Accession numbers may be provided directly, or '
    'in a file using the -f option. If provided directly, accessions may '
    'be comma or space separated, or hyphen-separated to specify a range, '
    'e.g. SRR123455, SRR123456, SRR123457 || SRR123455-SRR123457'
)
parser.add_argument(
    'accessions',
    nargs='*',
    help='List of SRA accession numbers (e.g. SRR123456 SRR123457)'
)
parser.add_argument(
    '-f',
    '--file',
    help='File with SRA accession numbers on separate lines'
)
parser.add_argument(
    '-k',
    '--keep_sra_files',
    help='Keep SRA files after dumping to FASTQ format',
    action='store_true'
)

if len(sys.argv) == 1:
    sys.exit(parser.print_help())

args = parser.parse_args()

if args.file:
    accs = acc_list_from_file(args.file)
else:
    accs = clean_args(args.accessions)


if len(accs) > 25:
    choice = input('Are you sure you want to run this on {} files? (y/n): '.format(len(accs)))
    if choice.lower() != 'y':
        sys.exit('Aborting run.')


cmds = {
    'paired': 'fastq-dump --defline-seq \'@$ac.$si:$sn[_$rn]/$ri\' --split-files ',
    'single': 'fastq-dump '
}

sra_files = []

for acc in accs:
    acc_filename = fetch_sra(acc)
    read_type = paired_check(acc_filename)    
    print('[#] Detected read type for {}: {}'.format(acc, read_type), 
          file=sys.stderr)
    sra_files.append(acc_filename)
    cmd = cmds[read_type] + acc_filename
    dump_args = shlex.split(cmd)
    print('[#] Running command \'{}\''.format(cmd), file=sys.stderr)
    subprocess.run(dump_args)
    # allow time for fastq-dump to print to stdout
    time.sleep(2)
    if not args.keep_sra_files:
        os.remove(acc_filename)

# if not args.keep_sra_files:
#     for sra in sra_files:
#         os.remove(sra)

print('[#] All commands finished. Exiting now.', file=sys.stderr)

sys.exit(0)

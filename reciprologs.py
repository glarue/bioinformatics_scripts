#!/usr/bin/env python3

"""
usage: reciprologs.py [-h] [-p [PARALLEL_PROCESSES]] [-t PERCENTAGE]
                      [--overwrite] [--one_to_one] [--logging]
                      file_1 file_2 ... [file_1 file_2 ... ...]
                      {blastn,blastp,blastx,tblastn,tblastx}

Uses BLAST to find reciprocal best hits between two or more files.

positional arguments:
  file_1 file_2 ...     files to use to build reciprolog sets (space
                        separated)
  {blastn,blastp,blastx,tblastn,tblastx}
                        type of BLAST program to run

optional arguments:
  -h, --help            show this help message and exit
  -p [PARALLEL_PROCESSES], --parallel_processes [PARALLEL_PROCESSES]
                        run the BLAST step using multiple parallel processes;
                        without specific input will use half available system
                        CPUs (default: None)
  -t PERCENTAGE, --query_percentage_threshold PERCENTAGE
                        require a specified fraction of the query length to
                        match in order for a hit to qualify (lowest allowable
                        percentage (default: None)
  --overwrite           overwrite existing BLAST files (instead of using to
                        bypass BLAST) (default: False)
  --one_to_one          remove any many-to-one reciprolog relationships in
                        each pairwise set, such that each member of each
                        pairwise comparison is only present exactly one time
                        in output (default: False)
  --logging             output a log of BLAST best-hit choice criteria
                        (default: False)

NOTE: Depends on pblast.py

"""
import sys
import subprocess
import os
import time
import argparse
import re
from operator import itemgetter
from multiprocessing import cpu_count
from collections import defaultdict
from itertools import combinations, permutations


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


def is_better(challenger, defender, seq_lengths=None):
    """
    Compares attributes of two dictionaries of BLAST
    hits for a given query to determine which is better.

    Returns the winning dictionary and reason if it's 
    better, otherwise False.
    
    """
    cbs = challenger['score']
    dbs = defender['score']
    # criteria: bitscore
    if cbs < dbs:
        return False
    elif cbs > dbs:
        return challenger, 'bitscore'
    elif cbs == dbs:
        # criteria --> e-value
        cev = challenger['evalue']
        dev = defender['evalue']
        if cev < dev: # lower is better
            return challenger, 'e-value'
        elif seq_lengths is not None:
            # criteria --> length
            # if scores are equal, check if sequence lenths
            # have been provided as an additional tiebreaking
            # criteria and look up the subject length to
            # see if there's a difference
            dn = defender['name']
            cn = challenger['name']
            if seq_lengths[cn] > seq_lengths[dn]:
                return challenger, 'length'
        else:
            return False
    else:
        return False


def get_top_hits(blast, paralogs=False, query_match=None, seq_lengths=None):
    results = {}
    # dictionary to store tie-broken matches
    win_ledger = defaultdict(lambda: defaultdict(set))
    with open(blast) as blst:
        for l in blst:
            new_best_hit = False
            if l.startswith("#"):
                continue
            (q, s, score, length, evalue) = parse_blast_line(
                l, "query", "subject", "bitscore", "length", "e")
            challenger = {
                'name': s,
                'score': score,
                'evalue': evalue
            }
            # do not consider hits to self if BLASTing against self,
            # but allow query/subject names to be the same
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
                defender = results[q]
                challenger_wins = is_better(
                    challenger, defender, seq_lengths)
                if challenger_wins:  # new hit is better
                    new_best_hit = True
                    defender_name = results[q]['name']
                    reason = challenger_wins[1]
                    loser_info = (defender_name, reason)
                    win_ledger[q]['losers'].add(loser_info)
                    win_ledger[q]['best'] = s
                    
            else:
                new_best_hit = True

            if new_best_hit is True:
                results[q] = {"name": s, "score": score, "evalue": evalue}

    return results, win_ledger


def get_reciprocals(d1, d2):
    """
    Takes two dictionaries of top BLAST hits,
    returns a list of tuples of all pairs that were
    reciprocal best hits, along with their bitscore
    values.

    """
    reciprologs = set()
    blast_combos = [(d1, d2), (d2, d1)]
    for first, second in blast_combos:
        for query, hit_info in first.items():
            best_hit = hit_info["name"]
            score = hit_info["score"]
            if best_hit in second:
                reciprocal_hit = second[best_hit]["name"]
                if query == reciprocal_hit:  # best hit refers back to query
                    r_score = second[best_hit]["score"]
                    hit_pair = sorted([query, best_hit])
                    score_tuple = tuple(sorted([score, r_score]))
                    hit_pair.append(score_tuple)
                    reciprologs.add(tuple(hit_pair))
    return sorted(reciprologs)


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


def concatenate(outname, file_list, clean=True):
    with open(outname, 'w') as outfile:
        for fn in file_list:
            with open(fn) as f:
                for l in f:
                    outfile.write(l)
    if clean:
        [os.remove(fn) for fn in file_list]


def make_ortho_dict(*orthos):
    """
    IN:
    [
        [('a', 'b'), ('a', 'c'), ('a', 'd')],
        [('b', 'c'), ('b', 'e'), ('b', 'f')],
        [('c', 'e'), ('c', 'f'), ('c', 'g')],
        [('z', 'x'), ('z', 'y'), ('z', 'w')]
    ]
    OUT:
    defaultdict(set,
            {'a': {'b', 'c', 'd'},
             'b': {'a', 'c', 'e', 'f'},
             'c': {'a', 'b', 'e', 'f', 'g'},
             'd': {'a'},
             'e': {'b', 'c'},
             'f': {'b', 'c'},
             'g': {'c'},
             'w': {'z'},
             'x': {'z'},
             'y': {'z'},
             'z': {'w', 'x', 'y'}})
    
    """
    collector = defaultdict(set)
    for o_list in orthos:
        for pair in o_list:
            for a, b in permutations(pair, 2):
                collector[a].add(b)
    return collector


def aggregate_dict_chained(ortho_dict):
    """
    IN:
    defaultdict(set,
            {'a': {'b', 'c', 'd'},
             'b': {'a', 'c', 'e', 'f'},
             'c': {'a', 'b', 'e', 'f', 'g'},
             'd': {'a'},
             'e': {'b', 'c'},
             'f': {'b', 'c'},
             'g': {'c'},
             'w': {'z'},
             'x': {'z'},
             'y': {'z'},
             'z': {'w', 'x', 'y'}})
    OUT:
    defaultdict(set, {'a': {'b', 'c', 'd', 'e', 'f', 'g'}, 'z': {'w', 'x', 'y'}})
    
    """
    changed = False
    processed = []
    master = defaultdict(set)
    for k, v in ortho_dict.items():
        if k in processed:
            continue
        processed.append(k)
        for v2 in v:
            if v2 == k:
                continue
            master[k].add(v2)
            processed.append(v2)
            if v2 not in ortho_dict:
                continue
            changed = True
            master[k].update(ortho_dict[v2])
    if changed is True:
        master = aggregate_dict_chained(master)
    return master


def aggregate_orthos_chained(orthos):
    """
    IN:
    [
        [('a', 'b'), ('a', 'c'), ('a', 'd')],
        [('b', 'c'), ('b', 'e'), ('b', 'f')],
        [('c', 'e'), ('c', 'f'), ('c', 'g')],
        [('z', 'x'), ('z', 'y'), ('z', 'w')]
    ]
    OUT:
    [['a', 'b', 'c', 'd', 'e', 'f', 'g'], ['w', 'x', 'y', 'z']]
    
    """
    o_dict = make_ortho_dict(*orthos)
    aggregated = aggregate_dict_chained(o_dict)
    ortho_groups = []
    for k, v in aggregated.items():
        combined = tuple(v) + (k,)
        ortho_groups.append(sorted(combined))
    
    return sorted(ortho_groups)


def every_member_match(members, m_dict):
    all_match = True
    for m in members:
        others = [e for e in members if e != m]
        if not others:
            return True
        for o in others:
            if m not in m_dict[o]:
                return False
            
    return all_match

    
def all_by_all_orthos(ortho_dict):
    full_groups = {}
    for k, v in ortho_dict.items():
        groups = []
        max_n = len(v)
        # special case where there is only one element to compare against
        if max_n == 1:
            if every_member_match(v, ortho_dict):
                groups.append(tuple(v))
        for i in range(2, max_n + 1):
            for g in combinations(v, i):
                if every_member_match(g, ortho_dict):
                    groups.append(g)
        full_groups[k] = groups
        
    all_sets = set()
    for k, group_list in full_groups.items():
        for g in group_list:
            full_group = sorted(set(list(g) + [k]))
            all_sets.add(tuple(full_group))
    
    return sorted(all_sets)


def aggregate_orthos_strict(orthos):
    """
    IN:
    [
        [('a', 'b'), ('a', 'c'), ('a', 'd')],
        [('b', 'c'), ('b', 'e'), ('b', 'f')],
        [('c', 'e'), ('c', 'f'), ('c', 'g')],
        [('z', 'x'), ('z', 'y'), ('z', 'w')]
    ]
    OUT:
    [
        ('a', 'b', 'c'), 
        ('a', 'd'), 
        ('b', 'c', 'e'), 
        ('b', 'c', 'f'), 
        ('c', 'g'), 
        ('w', 'z'), 
        ('x', 'z'), 
        ('y', 'z')
    ]
    
    """
    o_dict = make_ortho_dict(*orthos)
    aggregated = all_by_all_orthos(o_dict)
    
    return aggregated

def names_from_blastfile(blast_fn):
    file_pattern = r'(.+)-vs-(.+)\.t?blast[npx]'
    query_fn, subject_fn = re.findall(file_pattern, blast_fn)[0]

    return query_fn, subject_fn


def pair_reciprologs(query, subject, blast_type, qp, extra):
    # special case of BLASTing against self
    PARALOGS = query == subject

    q_lengths, s_lengths = {}, {}
    # we need to get sequence lengths for each file
    # for tie-breaking by length and qp
    for fa, target in zip([query, subject], [q_lengths, s_lengths]):
        target['query_match_threshold'] = qp
        for h, s in fasta_parse(fa):
            target[h] = len(s)

    # dictionaries are used as flag in top hits function
    # so need to be set here
    if not qp:
        qm_q, qm_s = {}, {}
    else:
        qm_q = q_lengths
        qm_s = s_lengths 

    # BLAST in both directions (unless PARALOGS)
    fw_names = unique_filenames(query, subject)
    rv_names = fw_names[::-1]

    fw_fn = '{}-vs-{}.{}'.format(*fw_names, blast_type)
    rv_fn = '{}-vs-{}.{}'.format(*rv_names, blast_type)

    # use existing BLAST output unless --overwrite is specified
    if not os.path.isfile(fw_fn) or OVERWRITE:
        forward_args = [BLAST, query, subject, blast_type, '-o', fw_fn] + extra
        subprocess.run(forward_args)
    else:
        print('[#] Using existing BLAST output \'{}\''.format(fw_fn))
    top_forward, fw_win_ledger = get_top_hits(
        fw_fn, PARALOGS, query_match=qm_q, seq_lengths=s_lengths)
    if PARALOGS:
        top_reverse = top_forward
        rv_win_ledger = defaultdict(list)
        ## don't need to run BLAST again; process the same file
        # top_reverse, rv_win_ledger = get_top_hits(
        #     fw_fn, PARALOGS, query_match=qm_q)
    else:
        if not os.path.isfile(rv_fn) or OVERWRITE:
            reverse_args = [BLAST, subject, query, blast_type, '-o', rv_fn] + extra
            subprocess.run(reverse_args)    
        else:
            print('[#] Using existing BLAST output \'{}\''.format(rv_fn))
        top_reverse, rv_win_ledger = get_top_hits(
            rv_fn, query_match=qm_s, seq_lengths=q_lengths)

    # Filter for reciprocal best hits
    reciprocal_hits = get_reciprocals(top_forward, top_reverse)

    # combine winners ledger

    win_ledger = {**fw_win_ledger, **rv_win_ledger}

    return reciprocal_hits, win_ledger

def remove_many_to_one(pairs):
    """
    Each element of >pairs< is a tuple: (hitA, hitB, (scoreX, scoreY))

    Takes a list of paired reciprocal hits (plus scores) and filters it
    such that each member of each pair only occur once, i.e. it removes
    any many-to-one hits, using the bitscores in the last element of
    each tuple.
    
    """
    uniques = {}
    to_remove = []
    for index, (a, b, scores) in enumerate(pairs):
        avg_score = sum(scores) / 2
        for e in [a, b]:
            if e not in uniques:
                uniques[e] = {'score': avg_score, 'index': index}
            elif uniques[e]['score'] >= avg_score:
                to_remove.append(index)
                continue
            else:
                to_remove.append(uniques[e]['index'])
                uniques[e] = {'score': avg_score, 'index': index}

    filtered = []
    for i, p in enumerate(pairs):
        if i in to_remove:
            names = p[0:2]
            print('Removed: {}'.format('\t'.join(names)), file=sys.stderr)
        else:
            filtered.append(p)

    return filtered

parser = argparse.ArgumentParser(
    description=(
        'Uses BLAST to find reciprocal best hits '
        'between two or more files.'),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    'input_files',
    metavar='file_1 file_2 ...',
    help='files to use to build reciprolog sets (space separated)',
    nargs='+'
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
    metavar='PERCENTAGE',
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
parser.add_argument(
    '--one_to_one',
    help=(
        'remove any many-to-one reciprolog relationships in each pairwise '
        'set, such that each member of each pairwise comparison is only '
        'present exactly one time in output'),
    action='store_true'
)
parser.add_argument(
    '--logging',
    help='output a log of BLAST best-hit choice criteria',
    action='store_true'
)

if len(sys.argv) == 1:
    sys.exit(parser.print_help())

t_start = time.time()

args, EXTRA_ARGS = parser.parse_known_args()

BLAST_TYPE = args.blast_type
PARALLEL = args.parallel_processes
# QUERY = args.file_1
# SUBJECT = args.file_2
INPUT_FILES = args.input_files
if len(INPUT_FILES) < 2:
    sys.exit('error: too few files specified (need >1)')
QUERY_PERCENTAGE = args.query_percentage_threshold
OVERWRITE = args.overwrite
ONE_TO_ONE = args.one_to_one
LOGGING = args.logging

BLAST = 'pblast.py'

# create a list with the flags/options to pass to the subsequence
# subprocess calls
call_options = {
    '-p': PARALLEL
}
optional = EXTRA_ARGS
for k, v in call_options.items():
    if v:
        optional.extend([str(k), str(v)])

# get reciprologs for each pairwise combo of files
pairwise_reciprolog_sets = []
for q, s in combinations(INPUT_FILES, 2):
    reciprolog_set, win_ledger = pair_reciprologs(
        q, s, BLAST_TYPE, QUERY_PERCENTAGE, optional)
    if LOGGING:
        ledger_file = '{}-{}.{}.log'.format(
            abbreviate(q), abbreviate(s), BLAST_TYPE)
        with open(ledger_file, 'w') as lf:
            for query, info in sorted(win_ledger.items()):
                winner = info['best']
                loser_tuples = info['losers']
                lf.write('>{}\t[{}]\n'.format(winner, query))
                for loser in sorted(loser_tuples):
                    lf.write('\t'.join(loser) + '\n')
    if ONE_TO_ONE is True:
        reciprolog_set = remove_many_to_one(reciprolog_set)
    # remove score info
    reciprolog_set = [tuple(x[0:2]) for x in reciprolog_set]

    pairwise_reciprolog_sets.append(reciprolog_set)


# reciprologs = aggregate_orthos_chained(pairwise_reciprolog_sets)
reciprologs = aggregate_orthos_strict(pairwise_reciprolog_sets)

basename = '-'.join([abbreviate(f) for f in INPUT_FILES])
out_file = '{}.{}.reciprologs'.format(basename, BLAST_TYPE)

with open(out_file, 'w') as out:
    for group in reciprologs:
        out.write('\t'.join(group) + '\n')

runtime = get_runtime(t_start)

print(
    '[#] Job finished in {}; {} reciprolog sets found: '
    '\'{}\''.format(runtime, len(reciprologs), out_file),
    file=sys.stderr)

sys.exit(0)

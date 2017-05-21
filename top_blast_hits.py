#!/usr/bin/env python3

"""
usage: top_blast_hits.py [-h] [-n NUMBER_OF_HITS] [-d] [-s] [-m] blastout_file

Reports the top n BLAST hits for each query in a blastout file

positional arguments:
  blastout_file

optional arguments:
  -h, --help            show this help message and exit
  -n NUMBER_OF_HITS, --number_of_hits NUMBER_OF_HITS
                        number of hits to report for each unique query
                        (default: 1)
  -d, --allow_duplicate_target_hits
                        allow multiple hits to the same subject to be included
                        (default: False)
  -s, --sort_by_name    sort output by name rather than bitscore (default:
                        False)
  -m, --memory_efficient
                        avoid reading entire dataset into memory at once, to
                        increase memory efficiency (default: False)

"""
import sys
from collections import defaultdict as dd
import argparse
import time

def get_runtime(start_time, p=5):
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

parser = argparse.ArgumentParser(description="Reports the top n BLAST hits "
                                             "for each query in a blastout file",
                                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("blastout_file")
parser.add_argument("-n", "--number_of_hits", type=int,
                    help="number of hits to report for each unique query",
                    default=1)
parser.add_argument("-d", "--allow_duplicate_target_hits",
                    action="store_true", default=False,
                    help="allow multiple hits to the same subject to be included")
parser.add_argument("-s", "--sort_by_name", action="store_true", default=False,
                    help="sort output by name rather than bitscore")
parser.add_argument("-m", "--memory_efficient", action="store_true", default=False,
                    help="avoid reading entire dataset into memory at once, "
                    "to increase memory efficiency")

args = parser.parse_args()

started = time.time()

# Figure out what arguments we're working with
BLASTFILE = args.blastout_file
N = args.number_of_hits
DUPES = args.allow_duplicate_target_hits
NAME_SORT = args.sort_by_name
MEM_REDUCE = args.memory_efficient

# Dictionary to store hits for each query
hitdict = dd(list)
total_hits = 0

# Do the data analysis at once or trimming as we go (~50% slower)
if not MEM_REDUCE:
    # fast
    with open(BLASTFILE) as blast:
        for hit in blast:
            bits = hit.strip().split('\t')
            query = bits[0]
            subject = bits[1]
            bitscore = float(bits[-1])
            # Attach the bitscore and subject to allow for sorting at the end
            index_tuple = (hit.strip(), bitscore, subject)
            hitdict[query].append(index_tuple)
            total_hits += 1

else:
    # slower
    with open(BLASTFILE) as blast:
        for hit in blast:
            bits = hit.strip().split('\t')
            query = bits[0]
            subject = bits[1]
            bitscore = float(bits[-1])
            # Attach the bitscore and subject to allow for sorting at the end
            index_tuple = (hit.strip(), bitscore, subject)
            hitdict[query].append(index_tuple)
            hitdict[query] = sorted(hitdict[query], key=lambda x: x[1], reverse=True)[:N]
            total_hits += 1

# The <file> argument directs output to standard error rather than
# standard out, to maintain clean redirection of the filtered query
# output to a file
print("[#] Indexed {} BLAST hits".format(total_hits), file=sys.stderr)

processed = 0

for query, hits in sorted(hitdict.items()):
    # Higher bitscores are better
    sorted_hits = sorted(hits, key=lambda x: x[1], reverse=True)
    filtered_hits = []
    if not DUPES:  # filter to unique subject hits only
        seen_subjects = set()
        for h in sorted_hits:
            subject = h[-1]
            if subject in seen_subjects:
                continue
            else:
                filtered_hits.append(h)
                seen_subjects.add(subject)
        final_hits = filtered_hits
    else:
        final_hits = sorted_hits
    if NAME_SORT:
        final_hits = sorted(final_hits[:N], key=lambda x: x[0])
    print("\n".join([x[0] for x in final_hits[:N]]))
    processed += len(final_hits[:N])

runtime = get_runtime(started)

sys.exit("[#] Processed {} hits in {}".format(processed, runtime))

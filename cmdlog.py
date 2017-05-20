#!/usr/bin/env python3

"""
A wrapper script that will write the command-line
args associated with any files generated to a log
file in the directory where the files were made.

"""

import sys
import os
from os import listdir
from os.path import isfile, join
import subprocess
import time
from datetime import datetime

def listFiles(mypath):
    """
    Return relative paths of all files in mypath

    """
    return [join(mypath, f) for f in listdir(mypath) if
            isfile(join(mypath, f))]

def read_log(log_file):
    """
    Reads a file history log and returns a dictionary
    of {filename: command} entries.

    Expects tab-separated lines of [time, filename, command]

    """
    entries = {}
    with open(log_file) as log:
        for l in log:
            l = l.strip()
            mod, name, cmd = l.split("\t")
            # cmd = cmd.lstrip("\"").rstrip("\"")
            entries[name] = [cmd, mod]
    return entries

def time_sort(t, fmt):
    """
    Turn a strftime-formatted string into a tuple
    of time info

    """
    parsed = datetime.strptime(t, fmt)
    return parsed

if len(sys.argv) == 1:
    sys.exit(
        'usage: cmdlog.py "[commands]"'
        '\n\nnote: all commands must be wrapped in double quotes')

ARGS = sys.argv[1]
ARG_LIST = ARGS.split()

# Guess where logfile should be put
if (">" or ">>") in ARG_LIST:
    # Get position after redirect in arg list
    redirect_index = max(ARG_LIST.index(e) for e in ARG_LIST if e in ">>")
    output = ARG_LIST[redirect_index + 1]
    output = os.path.abspath(output)
    out_dir = os.path.dirname(output)
elif ("cp" or "mv") in ARG_LIST:
    output = ARG_LIST[-1]
    out_dir = os.path.dirname(output)
else:
    out_dir = os.getcwd()

# Set logfile location within the inferred output directory
LOGFILE = out_dir + "/cmdlog_history.log"

# Get file list state prior to running
all_files = listFiles(out_dir)
pre_stats = [os.path.getmtime(f) for f in all_files]

# Run the desired external commands
subprocess.call(ARGS, shell=True)

# Get done time of external commands
TIME_FMT = "%Y-%m-%d@%H:%M:%S"
log_time = time.strftime(TIME_FMT)

# Get existing entries from logfile, if present
if LOGFILE in all_files:
    logged = read_log(LOGFILE)
else:
    logged = {}

# Get file list state after run is complete
post_stats = [os.path.getmtime(f) for f in all_files]
post_files = listFiles(out_dir)

# Find files whose states have changed since the external command
changed = [e[0] for e in zip(all_files, pre_stats, post_stats) if e[1] != e[2]]
new = [e for e in post_files if e not in all_files]
all_modded = list(set(changed + new))

if not all_modded:  # exit early, no need to log
    sys.exit(0)

# Replace files that have changed, add those that are new
for f in all_modded:
    name = os.path.basename(f)
    logged[name] = [ARGS, log_time]

# Write changed files to logfile
with open(LOGFILE, 'w') as log:
    for name, info in sorted(logged.items(), key=lambda x: time_sort(x[1][1], TIME_FMT)):
        cmd, mod_time = info
        if not cmd.startswith("\""):
            cmd = "\"{}\"".format(cmd)
        log.write("\t".join([mod_time, name, cmd]) + "\n")

sys.exit(0)

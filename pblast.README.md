# __pblast.py__

## __[tl;dr]__
`pblast.py` allows for significantly faster BLAST searches via parallelization, with some nice convenience functions to boot.

## __[details]__
BLASTing files can take a long time. `pblast.py` speeds up BLAST searches by first breaking the query file into chunks, and then BLASTing those chunks against the subject file in parallel.

For reasons that aren't entirely clear, this approach has significant speed gains over using the native `--num_threads` argument in the modern BLAST suite.

Additionally, `pblast.py` will auto-create the required database for a given BLAST run, and will format the input files to FASTA if necessary.

## __[example usage]__
To BLAST fileA against fileB using tblastx and 6 processes, simply do:

```pblast.py fileA fileB tblastx -p 6 > blast.output```

This will create the BLAST database, and convert filaA and fileB to FASTA files if they aren't already before running the search.

## __[background]__
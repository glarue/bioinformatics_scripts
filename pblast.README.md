# __palign__

## __[tl;dr]__
`palign` allows for significantly faster sequence alignment (via parallelization in the case of `BLAST`, and via the inherent speed of `DIAMOND`), with some nice convenience functions to boot.

## __[details]__
BLASTing files can take a long time. `palign` speeds up `BLAST+` searches by first breaking the query file into chunks, and then BLASTing those chunks 
against the subject file in parallel.

For reasons that aren't entirely clear, this approach has significant speed gains over using the native `--num_threads` argument in the modern `BLAST+` suite (note that this is _not_ true for DIAMOND; using `--threads` is recommended instead of `-p` in that case).

Additionally, `palign` will auto-create the required database for a given BLAST run, and will format the input files to FASTA if necessary.

Alternatively, `palign` can use the very fast aligner `DIAMOND` instead for the equivalent of `blastp` and `blastx` runs.

## __[example usage]__
To BLAST fileA against fileB using tblastx and 6 separate processes, simply do:

```palign fileA fileB tblastx -p 6 > blast.output```

This will create the BLAST database, and convert filaA and fileB to FASTA files if they aren't already before running the search.

## __[background]__

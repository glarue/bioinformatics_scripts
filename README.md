# General-use bioinformatics scripts

These scripts are relatively mature and well-documented, and should be useful
to a somewhat wider audience than just their author.

### `auto_dump.py`
Run NCBI's `fastq-dump` on one or more accession numbers faster than normal

### `cmdlog.py`
Keep a log of commands used and the files they produced

### `notify.py`
Automatically send an email to yourself at the end of a (usually long-running) command

### `pblast.py`
A parallel implementation of NCBI's BLAST+ suite (faster than their threading implementation)

### `random_lines.py`
Retrieve a random set of lines from a file (FASTA-format aware)

### `reciprologs.py`
Get the best reciprocal BLAST hits between two files

### `seq_from_fasta.py`
Retrieve a sequence (or sub-sequence) from a FASTA file

### `top_blast_hits.py`
Get only the top N hits for each query from a BLAST tabular output file

### `transcript_extractor.py`
Retrieve the coding sequences (defined by either CDS or exon entries) from a genome using an annotation file

### `translator.py`
Translate nucleotide sequence(s) into amino acid sequences
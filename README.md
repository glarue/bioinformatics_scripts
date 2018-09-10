# General-use bioinformatics scripts

These scripts are relatively mature and well-documented, and should be useful
to a somewhat wider audience than just their author.

### `cli_mail.py`
Send email over SSL from the command line

### `cmdlog.py`
Keep a log of commands used and the files they produce

### `fasterq_dump.py`
Run NCBI's `fastq-dump` on one or more accession numbers faster than normal

### `notify.py`
Automatically send an email to yourself at the end of a (usually long-running) command

### `pblast.py`
A parallel implementation of NCBI's BLAST+ suite (faster than their threading implementation for certain db/query size combinations)

### `random_lines.py`
Retrieve a random set of lines from a file (FASTA-format aware)

### `reciprologs.py`
Get the best reciprocal BLAST hits between two files

### `fasta_seq.py`
Retrieve one or more sequences (or sub-sequences) from a FASTA file

### `top_blast_hits.py`
Get only the top N hits for each query from a BLAST tabular output file

### `transcript_extractor.py`
Retrieve the coding sequences (defined by either CDS or exon entries) using a genome and annotation file

### `translator.py`
Translate nucleotide sequence(s) into amino acid sequences
# General-use bioinformatics scripts

These Pytyon scripts are relatively mature and well-documented, and may be useful to a somewhat wider audience than their author.

__Note:__ Some scripts require the [biogl](https://github.com/glarue/biogl) module to function properly; do 
`git clone --recursive https://github.com/glarue/bioinformatics_scripts.git` to pull the `biogl` directory as well

### `climail`
Send email over SSL from the command line

### `cmdlog`
Keep a log of commands used and the files they produce

### `fasterq_dump`
Run NCBI's `fastq-dump` on one or more accession numbers faster than normal

### `notify`
Automatically send an email to yourself at the end of a (usually long-running) command

### `palign`
A parallel implementation of NCBI's BLAST+ suite (faster than their threading implementation for certain db/query size combinations) or a convenience wrapper for DIAMOND

### `randl`
Retrieve a random set of lines from a file (FASTA-format aware)

### `reciprologs`
Get the best reciprocal BLAST/DIAMOND hits between two files

### `faseq`
Retrieve one or more sequences (or sub-sequences) from a FASTA file

### `tophits`
Get only the top N hits for each query from a BLAST tabular output file

### `cdseq`
Retrieve the coding sequences (defined by either CDS or exon entries) using a genome and annotation file

### `transeq`
Translate nucleotide sequence(s) into amino acid sequences

# General-use bioinformatics scripts

These Python scripts are relatively well-documented and battle tested, and may be useful to a somewhat wider audience than their author.

__Note:__ Some scripts require the [biogl](https://github.com/glarue/biogl) module to function properly; do 
`git clone --recursive https://github.com/glarue/bioinformatics_scripts.git` to pull the `biogl` directory as well

- [cdseq](#cdseq) | Extract transcript/coding sequences from a genome using an annotation file
- [climail](#climail) | Send an email over SSL from the command line
- [cmdlog](#cmdlog) | Log a shell commmand with optional note
- [faseq](#faseq) | Retrieve (sub)sequences from a FASTA file
- [fasterq_dump](#fasterq-dump) | Run NCBI's `fastq-dump` on one or more accession numbers faster than normal
- [notify](#notify) | Automatically send an email at the end of a (usually long-running) command
- [palign](#palign) | A parallel implementation of NCBI's BLAST+ suite (faster than their threaded implementation for certain db/query size combinations) or a convenience wrapper for DIAMOND
- [quickstats](#quickstats) | Calculate basic statistics on single-column input data
- [randl](#randl) | Retrieve a random set of lines from a file (FASTA-format aware)
- [reciprologs](#reciprologs) | Find the best reciprocal BLAST/DIAMOND hits between 2+ sets of sequences
- [revcomp](#revcomp) | Reverse-complement a nucleotide sequence
- [tophits](#tophits) | Get the top _n_ hits for each query sequence in a BLAST tabular output file
- [transeq](#transeq) | Translate a nucleotide sequence to amino acid sequence

### cdseq

```
usage: cdseq [-h] [-t] [-e] [-i] [-c] [-v] [-n] [--introns]
             [--truncate_introns int_length]
             genome annotation

Extract transcript/coding sequences from a genome using an annotation file

positional arguments:
  genome                genome file in FASTA format
  annotation            annotation file in GFF[3]/GTF format

optional arguments:
  -h, --help            show this help message and exit
  -t, --translate       translate the output sequence (default: False)
  -e, --exon            use exons instead of CDS entries to define coding
                        sequence (default: False)
  -i, --isoforms        allow multiple isoforms per gene, instead of only
                        longest (default: False)
  -c, --coord_based_isoforms
                        detect isoforms by overlapping coordindates in
                        addition to shared parent (useful for annotations
                        without gene entries) (default: False)
  -v, --verbose_headers
                        include coordinate info in output headers (default:
                        False)
  -n, --non_coding      include non-coding (intronic, UTR, etc.) sequence;
                        uses only the coordinates of the transcript itself
                        (incompatible with --introns) (default: False)
  --introns             include intron sequences in lowercase (incompatible
                        with -n) (default: False)
  --truncate_introns int_length
                        truncate intron sequences to length {int_length} (or
                        {int_length} -1 if odd) (default: None)
```
### climail

```
usage: climail [-h] -sadd SERVER_ADDRESS -fadd FROM_ADDRESS -pw PASSWORD -tadd
               TO_ADDRESS [-p PORT] [-s [SUBJECT]] [-b [BODY]] [--ID ID]
               [--no_host_tag] [--suffix_tag]

Send an email over SSL from the command line

optional arguments:
  -h, --help            show this help message and exit
  -p PORT, --port PORT  Port to use on email server (default: 587)
  -s [SUBJECT], --subject [SUBJECT]
                        A subject line for the email, wrapped in single quotes
                        ('') (default: None)
  -b [BODY], --body [BODY]
                        The body of the email, wrapped in single quotes ('')
                        (default: None)
  --ID ID               Arbitrary string to include in host tag (if present)
                        (default: None)
  --no_host_tag         disable host tag (time/server string) in subject
                        (default: False)
  --suffix_tag          place host tag at end of subject rather than beginning
                        (default: False)

required named arguments:
  -sadd SERVER_ADDRESS, --server_address SERVER_ADDRESS
                        Server address for outgoing email using the SSL
                        protocol (default: None)
  -fadd FROM_ADDRESS, --from_address FROM_ADDRESS
                        Full email address from which to send email (default:
                        None)
  -pw PASSWORD, --password PASSWORD
                        Account password for the 'from' address. Wrap in
                        single quotes ('') (default: None)
  -tadd TO_ADDRESS, --to_address TO_ADDRESS
                        Destination email address (default: None)
```
### cmdlog

```
usage: cmdlog [-h] shell commands [shell commands ...]

Log a shell commmand with optional note

positional arguments:
  shell commands  one or more commands to be sent to the shell (may require
                  single or double quotes)

optional arguments:
  -h, --help      show this help message and exit
```
### faseq

```
usage: faseq [-h] [-f INFO_FILE] [--flank FLANK] [-s SEPARATOR] [-u] [-e]
             [--full_header] [-t integer]
             FASTA_file [header strand start stop [label]
             [header strand start stop [label] ...]]

Retrieves sequences from a FASTA file using the following format: header
strand start stop [label]. If {label} is provided, will output in FASTA format
using {label} as header. Otherwise, outputs one sequence per line. If only
{header} is provided, will return the entire header/sequence pair.

positional arguments:
  FASTA_file            FASTA-formatted text file
  header strand start stop [label]
                        header strand start stop [label]

optional arguments:
  -h, --help            show this help message and exit
  -f INFO_FILE, --info_file INFO_FILE
                        file with information for sequences on separate lines.
  --flank FLANK         size of any desired flanking region around specified
                        sequence(s)
  -s SEPARATOR, --separator SEPARATOR
                        Character to use to separate flanking sequence (if
                        any) from main sequence, which defaults to {tab}.
  -u, --unformatted     leave original file formatting intact (e.g. don't join
                        multi-line entries into a single line)
  -e, --exclude_header  exclude the header line from the output (does not
                        apply to explicitly labeled queries)
  --full_header         match on full header string (including whitespace)
  -t integer, --truncate integer
                        truncate sequence to the nearest even integer <=
                        {integer} (equally from both ends)
```
### fasterq_dump

```
usage: fasterq_dump [-h] [-f FILE] [-k] [-t] [-w]
                    [accessions [accessions ...]]

A program to run fastq-dump sequentially on a list of accession numbers. Will
automatically detect if reads are single- or paired-end and will run fastq-
dump accordingly, adding proper ID suffixes as needed. Accession numbers may
be provided directly, or in a file using the -f option. If provided directly,
accessions may be comma or space separated, or hyphen-separated to specify a
range, e.g. SRR123455, SRR123456, SRR123457 or SRR123455-SRR123457

positional arguments:
  accessions            Space-separated list of SRA accession numbers (e.g.
                        SRR123456 SRR123457)

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  File with SRA accession numbers on separate lines
  -k, --keep_sra_files  Keep SRA files after dumping to FASTQ format
  -t, --trinity_compatible_ids
                        rename read IDs in paired-end data to ensure
                        downstream compatibility with Trinity
  -w, --overwrite       overwrite any matching local SRA files
```
### notify

```
usage: notify [-h] [-e EMAIL] [--add_email] [--view_config] [--ID ID]
              [external commands [external commands ...]]

Automatically sends an email to the specified address upon completion of the
specified command. Useful primarily for very long-running processes. In many
cases, the command being run (including the name of the other program) will
need to be placed in quotes. If no email address is provided with the -e flag,
a prompt will be displayed based upon the configuration file.

positional arguments:
  external commands     External commands to run, including external program
                        call (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -e EMAIL, --email EMAIL
                        the email address to notify (default: None)
  --add_email           add or change an email address in the config file
                        (default: False)
  --view_config         view the contents of the configuration file (default:
                        False)
  --ID ID               additional string to include in email subject
                        (default: None)
```
### palign

```
usage: palign [-h] [-p PARALLEL_PROCESSES] [-s] [-f OUTPUT_FORMAT]
              [-o OUTPUT_NAME] [-t THREADS] [-e E_VALUE] [-naf] [--clobber_db]
              query subject
              {diamondp,diamondx,blastn,blastp,blastx,tblastn,tblastx}

Align one file against another. Any arguments not listed here will be passed
to the chosen aligner unmodified.

positional arguments:
  query                 query file to be aligned
  subject               subject file to be aligned against
  {diamondp,diamondx,blastn,blastp,blastx,tblastn,tblastx}
                        type of alignment to run

optional arguments:
  -h, --help            show this help message and exit
  -p PARALLEL_PROCESSES, --parallel_processes PARALLEL_PROCESSES
                        run the alignment step using multiple parallel
                        processes (default: 1)
  -s, --single          disable parallel processing (default: False)
  -f OUTPUT_FORMAT, --output_format OUTPUT_FORMAT
                        integer output format for alignment results (default:
                        6)
  -o OUTPUT_NAME, --output_name OUTPUT_NAME
                        filename for results (otherwise, automatic based on
                        input) (default: None)
  -t THREADS, --threads THREADS
                        number of threads per process. Be careful when
                        combining this with multiple processes! (default: 1)
  -e E_VALUE, --e_value E_VALUE
                        e-value threshold to use for search (default: 1e-10)
  -naf, --no_auto_format
                        disable helper operations that auto-format
                        subject/query as needed and build database if not
                        already present (default: False)
  --clobber_db          create new database even if one already exists
                        (default: False)
```
### quickstats

```

```
### randl

```
usage: randl [-h] [-s SEED] file n

Retrieves a random sample of lines from a file.

positional arguments:
  file                  input file. If FASTA, will retrieve random
                        header/sequence pairs.
  n                     size of sample to retrieve.

optional arguments:
  -h, --help            show this help message and exit
  -s SEED, --seed SEED  seed value for the pseudo-random algorithm
```
### reciprologs

```
usage: reciprologs [-h] [-p PARALLEL_PROCESSES] [-q PERCENTAGE] [--chain]
                   [--overwrite] [--one_to_one] [--logging]
                   file_1 file_2 ... [file_1 file_2 ... ...]
                   {diamondp,diamondx,blastn,blastp,blastx,tblastn,tblastx}

Find reciprocal best hits between two or more files.

positional arguments:
  file_1 file_2 ...     files to use to build reciprolog sets (space
                        separated)
  {diamondp,diamondx,blastn,blastp,blastx,tblastn,tblastx}
                        type of alignment program to run

optional arguments:
  -h, --help            show this help message and exit
  -p PARALLEL_PROCESSES, --parallel_processes PARALLEL_PROCESSES
                        run the alignment step using multiple parallel
                        processes (default: 1)
  -q PERCENTAGE, --query_percentage_threshold PERCENTAGE
                        require a specified fraction of the query length to
                        match in order for a hit to qualify (lowest allowable
                        percentage (default: None)
  --chain               cluster reciprologs without requiring all-by-all
                        pairwise relationships, e.g. A-B, A-C, A-D --> A-B-C-D
                        (default: False)
  --overwrite           overwrite existing output files (instead of using to
                        bypass alignment) (default: False)
  --one_to_one          remove any many-to-one reciprolog relationships in
                        each pairwise set, such that each member of each
                        pairwise comparison is only present exactly one time
                        in output (default: False)
  --logging             output a log of best-hit choice criteria (default:
                        False)
```
### revcomp

```
usage: revcomp [-h] [-f FILE] [--no_mask] [--no_lower] [seq]

reverse-complement a nucleotide sequence

positional arguments:
  seq                   nucleotide sequence to be reverse-complemented
                        (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  file to reverse-complement (FASTA-compatible)
                        (default: None)
  --no_mask             don't mask non-ACTG characters with N (default: False)
  --no_lower            consider lowercase input invalid (will still reverse)
                        (default: False)
```
### tophits

```
usage: tophits [-h] [-n NUMBER_OF_HITS] [-e E_VALUE_CUTOFF] [-d] [-s] [-m]
               [-r]
               blast_file

Reports the top n BLAST hits for each query in a (tabular) blast output file

positional arguments:
  blast_file

optional arguments:
  -h, --help            show this help message and exit
  -n NUMBER_OF_HITS, --number_of_hits NUMBER_OF_HITS
                        number of hits to report for each unique query
                        (default: 1)
  -e E_VALUE_CUTOFF, --e_value_cutoff E_VALUE_CUTOFF
                        exclude hits with e-value greater than this value
                        (default: None)
  -d, --allow_duplicate_target_hits
                        allow multiple hits to the same subject to be included
                        (default: False)
  -s, --sort_by_name    sort output by name rather than bitscore (default:
                        False)
  -m, --memory_efficient
                        avoid reading entire dataset into memory at once, to
                        increase memory efficiency (default: False)
  -r, --redundant_ids   allow query and subject to share the same identifier
                        (default: False)
```
### transeq

```
usage: transeq [-h] [-v {short,long}] [-p {1,2}] [-r] [-s STOP_CHARACTER]
               [--no_mask]
               [sequence_input]

Translate nucleotide sequence into amino acid sequence

positional arguments:
  sequence_input        sequence string or file containing sequence(s)
                        (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -v {short,long}, --verbosity {short,long}
                        amino acid abbreviation length (e.g. Glu/Glutamine)
                        (default: single)
  -p {1,2}, --phase {1,2}
                        change the phase of translation (default: 0)
  -r, --reverse_complement
                        reverse-complement the sequence before translation
                        (default: False)
  -s STOP_CHARACTER, --stop_character STOP_CHARACTER
                        the string to use for stop codons (default: *)
  --no_mask             do not mask unrecognized codons with X; report
                        lowercased (default: False)
```

# __`fasterq_dump`__

## __[tl;dr]__

Downloads sequencing data in FASTQ format from the [NCBI Short Read Archive](https://www.ncbi.nlm.nih.gov/sra) using accession numbers.

## __[details]__

`fasterq_dump` is a wrapper script for NCBI's SRA-to-FASTQ conversion program `fastq-dump`, part of its [SRA-Tools package](http://ncbi.github.io/sra-tools/). `fasterq_dump` has the following advantages over vanilla `fastq-dump`:
* it is generally faster*
* it auto-detects read type – either single- or paired-end – and splits the output accordingly
* it formats the read IDs in paired-end data for compatability with Trinity (appending /1 or /2 to the ends of the IDs)
* it can be run on a list of multiple SRA accession numbers via either direct input on the command line or a separate text file
* it does not add files to the system SRA cache

\* Benchmarks:

* (PE) 5.4 M reads (1.1 Gbp): 4.81 min; ~1.12 M reads/min; ~229 Mbp/min
* (PE) 199.6 M reads (30.5 Gbp): 2.594 hr; ~1.28 M reads/min; ~196 Mbp/min
* (SE) 43.3 M reads (2.2 Gbp): 11.064 min; ~3.91 M reads/min; ~199 Mbp/min

* In many cases it may also run faster than another option, `parallel-fastq-dump.py` (~37% faster for an 11 GB PE dataset)

## __[example usage]__

To obtain RNA-seq reads for [this small ebola dataset](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run=SRR3359559), run `fasterq_dump` on the corresponding accession number (SRR3359559), which will result in the following output:

```console
$ fasterq_dump SRR3359559
--2017-07-04 13:38:55--  ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR335/SRR3359559/SRR3359559.sra
           => ‘SRR3359559.sra’
Resolving ftp-trace.ncbi.nlm.nih.gov (ftp-trace.ncbi.nlm.nih.gov)... 130.14.250.10, 2607:f220:41e:250::11
Connecting to ftp-trace.ncbi.nlm.nih.gov (ftp-trace.ncbi.nlm.nih.gov)|130.14.250.10|:21... connected.
Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
==> TYPE I ... done.  ==> CWD (1) /sra/sra-instant/reads/ByRun/sra/SRR/SRR335/SRR3359559 ... done.
==> SIZE SRR3359559.sra ... 667481
==> PASV ... done.    ==> RETR SRR3359559.sra ... done.
Length: 667481 (652K) (unauthoritative)

SRR3359559.sra      100%[===================>] 651.84K   811KB/s    in 0.8s    

2017-07-04 13:38:57 (811 KB/s) - ‘SRR3359559.sra’ saved [667481]

[#] Detected read type for SRR3359559: paired
[#] Running command 'fastq-dump --defline-seq '@$ac.$si:$sn[_$rn]/$ri' --split-files SRR3359559.sra'
Read 4481 spots for SRR3359559.sra
Written 4481 spots for SRR3359559.sra
[#] All commands finished. Exiting now.
```

To retrieve multiple runs from the same experiment, multiple accessions may be listed as command line arguments,

```
$ fasterq_dump SRR3359557-SRR3359559
```

which will result in a corresponding set of paired FASTQ files:

```
-rw-r--r--  1 glarue glarue  17M Jul  4 13:53 SRR3359557_1.fastq
-rw-r--r--  1 glarue glarue  17M Jul  4 13:53 SRR3359557_2.fastq
-rw-r--r--  1 glarue glarue 5.1M Jul  4 13:53 SRR3359558_1.fastq
-rw-r--r--  1 glarue glarue 5.1M Jul  4 13:53 SRR3359558_2.fastq
-rw-r--r--  1 glarue glarue 1.2M Jul  4 13:53 SRR3359559_1.fastq
-rw-r--r--  1 glarue glarue 1.2M Jul  4 13:53 SRR3359559_2.fastq

```

Run `fasterq_dump` without arguments for further usage information.


## __[background]__

Based on in-house testing, downloading the SRA file directly (via `wget`) and then using `fastq-dump` locally is often significantly faster than using the automated downloading capabilities of `fastq-dump`. In addition, users may often want to retrieve multiple datasets in one fell swoop; while this could be achieved using a relatively simple Bash loop, user familarity with Bash is variable and the syntax of `fasterq_dumpy.py` is likely simpler.

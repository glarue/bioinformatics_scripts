# __fasta_seq__

## __[tl;dr]__

Retrieves one or more sub-sequences from a FASTA file based on location, coordinate and strand information. Will reverse-complement sub-sequences on the negative strand.

## __[details]__

`faseq` allows for the extraction of either header/sequence pairs in 
their entirety or sub-sequences thereof from a FASTA file. In a large genome FASTA, one might wish to isolate the sequence for scaffold42. This could of course be done using `grep -A 1 'scaffold42'`, unless (as is unfortunately often the case) the formatting of the file did not place the entire sequence on a single line.

If, however, one were interested in a single exon's sequence, such an operation proves more difficult with basic command-line tools, as it involves parsing a sub-sequence from within the parent sequence using positional information.

In order to pull a header/sequence pair, `faseq` requires only the header name; to retrieve sub-sequences, the minimal set of information required in addition to the FASTA file is (in order): `header strand start stop`

## __[example usage]__

For illustrative purposes, here is a dummy FASTA file:

```
$ more dummy.fasta
>one
11111111111111111111111111111111111111111111111111111111111111111111111111111111
>two
22222222222222222222222222222222222222222222222222222222222222222222222222222222
22222222222222222222222222222222222222222222222222222222222222222222222222222222
>three
33333333333333333333333333333333333333333333333333333333333333333333333333333333
33333333333333333333333333333333333333333333333333333333333333333333333333333333
33333333333333333333333333333333333333333333333333333333333333333333333333333333
>four
44444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
>one-to-nine
123456789
```

Note the variable formatting. The first entry has a single 80 character sequence, where the second and third entries contain multiple lines of sequence per entry. The fourth is better behaved and has all of its sequence on a single line. The fifth sequence will be used to demonstrate sub-sequence retrieval.

### __whole entries__

To pull out the sequence for header `one`, the syntax is simply

```
$ faseq dummy.fasta one
>one
11111111111111111111111111111111111111111111111111111111111111111111111111111111
```

For headers that have sequences across multiple lines, `faseq` will output the sequence in a single line by default:

```
$ faseq dummy.fasta three
>three
333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
```

Use `-u` to maintain the original formatting during output:

```
$ faseq dummy.fasta three -u
>three
33333333333333333333333333333333333333333333333333333333333333333333333333333333
33333333333333333333333333333333333333333333333333333333333333333333333333333333
33333333333333333333333333333333333333333333333333333333333333333333333333333333
```

### __sub-sequences__

To retrieve a sub-sequence from an entry, one must supply the following: `header strand start stop`. Optionally, a final `label` value may be supplied to name the sub-sequence, which will result in FASTA-formatted output; without labels, the output is one sub-sequence per line.

The `strand` value is used to determine whether the sub-sequence should be reverse complemented. If no reverse complementing is desired, simply use `+` for the `strand` value.

Note that the coordinate values are _inclusive_, and 1-indexed (as is standard in genomic coordinates).

To get the sub-sequence `4567` from header `one-to-nine` in our dummy FASTA:

```
$ faseq dummy.fasta one-to-nine + 4 7
4567
```

To retrieve multiple sequences in one go, include the appropriate information for each sequence on separate lines in a file,

```
$ more seqinfo.txt          
two + 1 3
one-to-nine + 1 3
```

and run `faseq` on the file using the `-f` flag:

```
$ faseq dummy.fasta -f seqinfo.txt 
222
123
```

## __[background]__

If one works with genomic data long enough, one is bound to want to check one or more sequences from within a FASTA file. The Roy lab has `get.pl`, but this only fetches a single header/sequence pair. Manually parsing out sequences of various genomic features a handful of times was enough to motivate the creation of this script.

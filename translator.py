#!/usr/bin/env python3
"""
usage: translator.py [-h] [-v {short,long}] [-p {1,2}] [-r]
                     [-s STOP_CHARACTER] [--no_mask]
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
                        the string to use for stop codons (default: None)
  --no_mask             do not mask unrecognized codons with X; report
                        lowercased (default: False)
"""
import sys
import os
import argparse

def fasta_parse(fasta, separator=None):
    """
    Iterator which takes FASTA as input. Yields
    header/values list pairs. If separator arg,
    will output a string of values joined by
    separator character (e.g separator="" will
    return a string with no spaces or intervening
    characters).
    """
    header, seq = None, []
    with open(fasta) as fastaFile:
        for line in fastaFile:
            if line.startswith('>'):
                if header:
                    if separator != None:
                        seq = separator.join(str(e) for e in seq)
                    yield header, seq
                header = line.strip().lstrip('>')
                seq = []
            elif line.startswith('#'):
                continue
            else:
                seq.append(line.rstrip('\n'))
        if separator != None:
            seq = separator.join(str(e) for e in seq)
        yield header, seq

def translate_seq(
    string, 
    verbosity="single", 
    phase=0, 
    rev=False, 
    stop_char='*'):
    codonMap = {
        'AAA': ('K', 'Lys', 'Lysine'),
        'AAC': ('N', 'Asn', 'Asparagine'),
        'AAG': ('K', 'Lys', 'Lysine'),
        'AAT': ('N', 'Asn', 'Asparagine'),
        'ACA': ('T', 'Thr', 'Threonine'),
        'ACC': ('T', 'Thr', 'Threonine'),
        'ACG': ('T', 'Thr', 'Threonine'),
        'ACT': ('T', 'Thr', 'Threonine'),
        'AGA': ('R', 'Arg', 'Arginine'),
        'AGC': ('S', 'Ser', 'Serine'),
        'AGG': ('R', 'Arg', 'Arginine'),
        'AGT': ('S', 'Ser', 'Serine'),
        'ATA': ('I', 'Ile', 'Isoleucine'),
        'ATC': ('I', 'Ile', 'Isoleucine'),
        'ATG': ('M', 'Met', 'Methionine'),
        'ATT': ('I', 'Ile', 'Isoleucine'),
        'CAA': ('Q', 'Gln', 'Glutamine'),
        'CAC': ('H', 'His', 'Histidine'),
        'CAG': ('Q', 'Gln', 'Glutamine'),
        'CAT': ('H', 'His', 'Histidine'),
        'CCA': ('P', 'Pro', 'Proline'),
        'CCC': ('P', 'Pro', 'Proline'),
        'CCG': ('P', 'Pro', 'Proline'),
        'CCT': ('P', 'Pro', 'Proline'),
        'CGA': ('R', 'Arg', 'Arginine'),
        'CGC': ('R', 'Arg', 'Arginine'),
        'CGG': ('R', 'Arg', 'Arginine'),
        'CGT': ('R', 'Arg', 'Arginine'),
        'CTA': ('L', 'Leu', 'Leucine'),
        'CTC': ('L', 'Leu', 'Leucine'),
        'CTG': ('L', 'Leu', 'Leucine'),
        'CTT': ('L', 'Leu', 'Leucine'),
        'GAA': ('E', 'Glu', 'Glutamate'),
        'GAC': ('D', 'Asp', 'Aspartate'),
        'GAG': ('E', 'Glu', 'Glutamate'),
        'GAT': ('D', 'Asp', 'Aspartate'),
        'GCA': ('A', 'Ala', 'Alanine'),
        'GCC': ('A', 'Ala', 'Alanine'),
        'GCG': ('A', 'Ala', 'Alanine'),
        'GCT': ('A', 'Ala', 'Alanine'),
        'GGA': ('G', 'Gly', 'Glycine'),
        'GGC': ('G', 'Gly', 'Glycine'),
        'GGG': ('G', 'Gly', 'Glycine'),
        'GGT': ('G', 'Gly', 'Glycine'),
        'GTA': ('V', 'Val', 'Valine'),
        'GTC': ('V', 'Val', 'Valine'),
        'GTG': ('V', 'Val', 'Valine'),
        'GTT': ('V', 'Val', 'Valine'),
        'TAC': ('Y', 'Tyr', 'Tyrosine'),
        'TAT': ('Y', 'Tyr', 'Tyrosine'),
        'TCA': ('S', 'Ser', 'Serine'),
        'TCC': ('S', 'Ser', 'Serine'),
        'TCG': ('S', 'Ser', 'Serine'),
        'TCT': ('S', 'Ser', 'Serine'),
        'TGC': ('C', 'Cys', 'Cysteine'),
        'TGG': ('W', 'Trp', 'Tryptophan'),
        'TGT': ('C', 'Cys', 'Cysteine'),
        'TTA': ('L', 'Leu', 'Leucine'),
        'TTC': ('F', 'Phe', 'Phenylalanine'),
        'TTG': ('L', 'Leu', 'Leucine'),
        'TTT': ('F', 'Phe', 'Phenylalanine'),
        'TAA': (stop_char, stop_char, 'STOP'),
        'TAG': (stop_char, stop_char, 'STOP'),
        'TGA': (stop_char, stop_char, 'STOP'),
        '---': ('-', '-', 'GAP')
    }

    def get_codons(s, p=0):
        codons = [s[i + p:i + p + 3] for i in range(0, len(s), 3)]
        return [c for c in codons if len(c) == 3]
    string = string.strip().replace(" ", "").upper()  # remove spaces if present
    if rev is True:
        string = reverse_complement(string)
    codons = get_codons(string, p=phase)
    if not codons:
        return '?', '?'
    stop_codons = ("TAA", "TAG", "TGA")
    verbosityD = {"single": 0, "short": 1, "long": 2}
    amino_acids = []
    v = verbosityD[verbosity]
    if v == 0:
        joinChar = ""
    else:
        joinChar = "-"
    for c in codons:
        try:
            amino_acids.append(codonMap[c.upper()][v])
        except KeyError:
            if NO_MASK:
                amino_acids.append(c.lower())
            else:
                amino_acids.append('X')
    return joinChar.join(amino_acids)

def parse_seq(s):
    if s.startswith(">"):
        print(s.strip())
        return
    aas = translate_seq(
        s, 
        verbosity=v_level, 
        phase=phase_choice,
        rev=rev_choice,
        stop_char=STOP_CHAR)
    print(aas, flush=True)

def reverse_complement(seq):
    """
    Returns reverse complement of seq, with
    any non-ACTG characters replaced with Ns

    """
    transform = {'A': 'T',
                 'T': 'A',
                 'C': 'G',
                 'G': 'C',
                 'N': 'N'}
    try:
        comp = [transform[e] for e in seq]
    except KeyError:  # non-ATCGN characters in seq
        seq = [e if e in "ACTGN" else "N" for e in seq]
        comp = [transform[e] for e in seq]
    rev_comp = comp[::-1]
    rev_comp_string = ''.join(rev_comp)
    return rev_comp_string

def parse_file(seqFile):
    """
    Do translation if input is a file.

    """
    isFasta = False
    with open(seqFile) as f:
        if any(l.startswith(">") for l in f):
            isFasta = True

    if isFasta:
        for h, s in fasta_parse(seqFile, separator=""):
            aas = translate_seq(
                s, 
                verbosity=v_level, 
                phase=phase_choice, 
                rev=rev_choice, 
                stop_char=STOP_CHAR)
            print(">" + h + "\t")
            print(aas, flush=True)
    else:
        with open(seqFile) as f:
            for l in f:
                aas = translate_seq(
                    l, 
                    verbosity=v_level, 
                    phase=phase_choice, 
                    rev=rev_choice, 
                    stop_char=STOP_CHAR)
                print(aas, flush=True)

parser = argparse.ArgumentParser(
    description='Translate nucleotide sequence into amino acid sequence',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
    'sequence_input',
    nargs='?',
    help='sequence string or file containing sequence(s)'
)
parser.add_argument(
    '-v',
    '--verbosity',
    help='amino acid abbreviation length (e.g. Glu/Glutamine)',
    choices=['short', 'long'],
    default='single'
)
parser.add_argument(
    '-p',
    '--phase',
    help='change the phase of translation',
    type=int,
    choices=[1, 2],
    default=0
)
parser.add_argument(
    '-r',
    '--reverse_complement',
    help='reverse-complement the sequence before translation',
    action='store_true'
)
parser.add_argument(
    '-s',
    '--stop_character',
    help='the string to use for stop codons',
    type=str,
    default='*'
)
parser.add_argument(
    '--no_mask',
    action='store_true',
    help='do not mask unrecognized codons with X; report lowercased'
)

args = parser.parse_args()

# this is convoluted, but allows for piping in of sequence
seq_input = args.sequence_input
if not seq_input:
    if not sys.stdin.isatty():
        seq_input = sys.stdin
    else:
        parser.print_help()
        sys.exit(0)
        
v_level = args.verbosity
phase_choice = args.phase
rev_choice = args.reverse_complement
NO_MASK = args.no_mask
STOP_CHAR = args.stop_character

try:
    if os.path.isfile(seq_input):
        parse_file(seq_input)  
    else:
        parse_seq(seq_input)
except TypeError:
    for line in seq_input.readlines():
        parse_seq(line)


sys.exit(0)

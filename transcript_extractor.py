#!/usr/bin/env python3
"""
usage: transcript_extractor.py [-h] [-t] [-e] [-i] [-c] [-v] genome annotation

Extract transcript/coding sequences from an annotation + genome file
combination

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

"""

import sys
import argparse
from collections import defaultdict
from operator import itemgetter


class GFFLineInfo(object):
    """
    Takes a gff3/gtf/gff annotation line, and returns available metadata
    about it.

    """
    def __init__(self, line, line_number):
        self.bits = self.__split_on_tabs(line)
        try:
            self.region = self.bits[0]
            self.start = min(map(int, self.bits[3:5]))
            self.stop = max(map(int, self.bits[3:5]))
            self.strand = self.bits[6]
            self.infostring = self.bits[8]
            self.feat_type = self.get_type()
            self.parent = self.get_parent()
            self.name = self.get_ID()
            self.line_number = line_number
            self.protein_coding = self.is_protein_coding()
        except TypeError:
            raise

    @staticmethod
    def __split_on_tabs(l, n=9):
        """
        Checks for a valid line (does not start with #).

        Splits valid line on tabs, and returns a list of the bits
        if the list length is <<n>>. Setting <<n>> to None will return
        the split line regardless of length.

        """
        if l.startswith('#'):
            return None
        l = l.strip()
        columns = l.split("\t")
        if n and len(columns) < n:
            return None
        return columns

    @staticmethod
    def __field_match(infostring, tags, delimiter):
        info_bits = infostring.split(delimiter)
        try:
            match = next(
                e for e in info_bits
                if any(p.lower() in e.lower() for p in tags))
        except StopIteration:  # no matches found
            return None
        if "=" in match:
            substring = match.split("=")[1]
        else:
            substring = match.split()[1]
        return substring.strip("\"")

    def get_type(self, delimiter=';'):
        """
        Classifies annotation lines into type categories,
        taking into account edge cases like 'processed transcript'
        and 'pseudogenic transcript'.
        
        """
        og_type = self.bits[2].lower()
        if og_type == 'mrna':
            og_type = 'transcript'
        if og_type in ('gene', 'transcript', 'exon', 'cds'):
            return og_type

        disqualifying = ['utr', 'start', 'stop']
        if any(kw in og_type for kw in disqualifying):
            return None

        # Not an obvious type, so search for features of transcripts
        # and genes in infostring to try to infer type
        gene_tags = ["gene_id", "geneId"]
        transcript_tags = ["transcriptId", "transcript_ID"]
        # Transcripts first because genes shouldn't have transcript IDs,
        # but transcripts may have gene IDs
        for ftype, tags in zip(
                ['transcript', 'gene'], [transcript_tags, gene_tags]
        ):
            match = self.__field_match(self.infostring, tags, delimiter)
            if match:
                return ftype

    def is_protein_coding(self):
        """
        Infer whether a annotation feature is protein coding or
        not by a couple of rough heuristics.
        
        """
        og_type = self.bits[2].lower()
        # Clear cases where it should be protein-coding by default
        if og_type in ('gene', 'transcript', 'mrna'):
            return True
        infostring = self.infostring
        # This amounts to a blacklist of all things that don't explicitly
        # state that they're protein coding
        if 'protein_coding' in infostring:
            return True
        else:
            return False

    def get_ID(self, delimiter=";"):
        """
        Finds the ID of a given annotation file line.

        """
        # Constrain feature types to simplify indexing
        feat_type_converter = {"cds": "exon", "mrna": "transcript"}
        feat_type = self.feat_type
        if feat_type in feat_type_converter:
            feat_type = feat_type_converter[feat_type]
        # all get lowercased in the comparison
        child_tags = ["ID="]
        # if is no 'ID=', should reference self via others
        gene_tags = ["ID=", "gene_id", "geneId"]
        transcript_tags = ["ID=", "transcriptId", "transcript_ID"]
        tag_selector = {
            "gene": gene_tags,
            "transcript": transcript_tags,
            "exon": child_tags
        }
        try:
            tags = tag_selector[feat_type]
        except KeyError:
            return None
        infostring = self.infostring
        match = self.__field_match(infostring, tags, delimiter)
        return match

    def get_parent(self, delimiter=";"):
        """
        Retrieves parent information from an annotation line.
        
        """
        feat_type_converter = {"cds": "exon", "mrna": "transcript"}
        feat_type = self.feat_type
        if feat_type in feat_type_converter:
            feat_type = feat_type_converter[feat_type]
        child_tags = [
            "Parent=", "proteinId", "protein_ID", "transcriptId",
            "transcript_ID"
        ]
        transcript_tags = ["Parent=", "geneId", "gene_ID"]
        gene_tags = ["Parent="]
        tag_selector = {
            "gene": gene_tags,
            "transcript": transcript_tags,
            "exon": child_tags
        }
        try:
            tags = tag_selector[feat_type]
        except KeyError:
            return None
        infostring = self.infostring
        match = self.__field_match(infostring, tags, delimiter)
        if not match and feat_type == "transcript":
            match = self.get_ID()
        return match


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


def get_subseq(region_seq, start, stop):
    # Correct for 1-based indexing in start and stop
    start -= 1
    # Pull sequence and reverse if necessary
    seq = region_seq[start:stop]

    return seq


def fasta_parse(fasta, delimiter=">", separator="", trim_header=True):
    """
    Iterator which takes FASTA as input. Yields
    header/value pairs. Separator will be
    used to join the return value; use separator=
    None to return a list.

    If trim_header, parser will return the
    FASTA header up to the first space character.
    Otherwise, it will return the full, unaltered
    header string.

    """
    header, seq = None, []
    with open(fasta) as f:
        for line in f:
            if line.startswith(delimiter):
                if header:  # associate accumulated seq with header
                    if separator is not None:
                        seq = separator.join(str(e) for e in seq)
                    yield header, seq
                # Assign a new header
                header = line.strip().lstrip(delimiter)
                if trim_header:
                    header = header.split()[0]
                # Clear seq for new round of collection
                seq = []
            elif line.startswith('#'):
                continue
            else:
                if line.strip():  # don't collect blank lines
                    seq.append(line.rstrip('\n'))
        if separator is not None:  # make string
            seq = separator.join(str(e) for e in seq)
        yield header, seq


def translate_seq(string, verbosity="single", phase=0):
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
        'TAG': ('*', '*', 'STOP'),
        'TGA': ('*', '*', 'STOP'),
        'TAA': ('*', '*', 'STOP'),
    }

    def _get_codons(s, p=0):
        for i in range(0, len(s), 3):
            yield s[i + p:i + p + 3]

    verbosityD = {"single": 0, "short": 1, "long": 2}
    string = string.replace(" ", "")  # remove spaces if present
    codons = _get_codons(string, p=phase)
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
            amino_acids.append(c.lower())
    return joinChar.join(amino_acids)

def get_transcripts(gff, child_type):
    transcripts = defaultdict(dict)
    child_type_found = False
    regions_with_content = []
    with open(gff) as annot:
        for ln, line in enumerate(annot):
            feat = None  ###!!! testing
            try:
                feat = GFFLineInfo(line, ln)
            except TypeError:
                continue
            if feat.feat_type == 'transcript':
                info_dict = {
                    'name': feat.name,
                    'strand': feat.strand,
                    'parent': feat.parent,
                    'coords': (feat.start, feat.stop)
                }
                if feat.name not in transcripts[feat.region]:
                    transcripts[feat.region][feat.name] = {
                        'info': info_dict,
                        'children': []}
                else:  # made by child
                    transcripts[feat.region][feat.name]['info'] = info_dict
            elif feat.feat_type == child_type:
                if not child_type_found:
                    child_type_found = True
                parent = feat.parent
                start = feat.start
                stop = feat.stop
                region = feat.region
                if parent not in transcripts[region]:
                    transcripts[region][parent] = {'children': []}
                entry = transcripts[region][parent]['children']
                entry.append((start, stop))
                if region not in regions_with_content:
                    regions_with_content.append(region)
    transcripts = {
        k: v for k, v in transcripts.items() if k in regions_with_content}
    
    return transcripts, child_type_found

def overlap_check(a, b):
    """
    a and b are both (start, stop) coord pairs
    
    """
    lowest = min([a, b], key=lambda x: min(x))
    highest = max([a, b], key=lambda x: min(x))
    if min(highest) <= max(lowest):
        return True
    else:
        return False

# def longest_isoforms_by_coord(transcript_dict):
#     """
#     Identifies longest isoforms, and returns a dictionary.
    
#     """
#     # identify longest isoforms
#     # sort by length, then use overlap() function to determine if
#     # subsequent transcripts are isoforms of longest version and skip
#     # if they are.
#     longest_isoforms = {}
#     for region, transcripts in transcript_dict.items():
#         if region not in longest_isoforms:
#             longest_isoforms[region] = {}
#         seen_coords = set()
#         # sort by longest transcripts first
#         for name, meta in sorted(transcripts.items(), 
#         key=lambda x: coding_length(x[1]['children']), reverse=True):
#             coords = meta['info']['coords']
#             if not any(overlap_check(coords, c) for c in seen_coords):
#                 length = coding_length(meta['children'])
#                 meta['info']['length'] = length
#                 longest_isoforms[region][name] = meta
#             seen_coords.add(coords)
            
#     return longest_isoforms

def longest_isoforms(transcript_dict, use_coords=False):
    """
    Identifies longest isoforms, and returns a dictionary.
    
    """
    # identify longest isoforms
    # sort by length, then use either gene name or overlap() function 
    # to determine if subsequent transcripts are isoforms of longest 
    # version and skip if they are
    longest_isoforms = {}
    seen_genes = set()
    for region, transcripts in transcript_dict.items():
        if region not in longest_isoforms:
            longest_isoforms[region] = {}
        seen_coords = set()
        # sort by longest transcripts first
        for name, meta in sorted(transcripts.items(), 
        key=lambda x: coding_length(x[1]['children']), reverse=True):
            try:
                gene = meta['info']['parent']
            except:
                print(name, meta, file=sys.stderr)
            if gene not in seen_genes:
                if use_coords:
                    # skip those overlapping longer transcripts
                    # even if they have unique gene name
                    coords = meta['info']['coords']
                    if any(overlap_check(coords, c) for c in seen_coords):
                        seen_coords.add(coords)
                        continue
                seen_genes.add(gene)
                length = coding_length(meta['children'])
                meta['info']['length'] = length
                longest_isoforms[region][name] = meta
    
    return longest_isoforms

def finalize_transcripts(transcript_dict):
    """
    Reformats transcripts, and returns a dictionary.
    
    """
    # identify longest isoforms
    # sort by length, then use either gene name or overlap() function 
    # to determine if subsequent transcripts are isoforms of longest 
    # version and skip if they are
    finalized = defaultdict(dict)
    for region, transcripts in transcript_dict.items():
        for name, meta in sorted(transcripts.items()):
            try:
                gene = meta['info']['parent']
            except:
                print(name, meta, file=sys.stderr)
            length = coding_length(meta['children'])
            meta['info']['length'] = length
            finalized[region][name] = meta
    
    return finalized

# def longest_isoforms_by_gene(transcript_dict):
#     """
#     Identifies longest isoforms, and returns a dictionary.
    
#     """
#     # identify longest isoforms
#     # sort by length, then use overlap() function to determine if
#     # subsequent transcripts are isoforms of longest version and skip
#     # if they are.
#     longest_isoforms = {}
#     seen_genes = set()
#     for region, transcripts in transcript_dict.items():
#         if region not in longest_isoforms:
#             longest_isoforms[region] = {}
#         # sort by longest transcripts first
#         for name, meta in sorted(transcripts.items(), 
#         key=lambda x: coding_length(x[1]['children']), reverse=True):
#             gene = meta['info']['parent']
#             if gene not in seen_genes:
#                 seen_genes.add(gene)
#                 length = coding_length(meta['children'])
#                 meta['info']['length'] = length
#                 longest_isoforms[region][name] = meta
    
#     return longest_isoforms


## old way, works
# def longest_isoforms_by_gene(transcript_dict):
#     """
#     Identifies longest isoforms, and returns a dictionary.
    
#     """
#     # identify longest isoforms
#     longest_isoforms = {}
#     for region, transcripts in transcript_dict.items():
#         for name, meta in transcripts.items():
#             gene = meta['info']['parent']
#             length = coding_length(meta['children'])
#             if gene in longest_isoforms:
#                 if length <= longest_isoforms[gene]['length']:
#                     continue
#             longest_isoforms[gene] = {
#                 'transcript': name,
#                 'length': length}

#     # filter based on longest found transcripts
#     filtered = {}
#     for region, transcripts in transcript_dict.items():
#         if region not in filtered:
#             filtered[region] = {}
#         for name, meta in transcripts.items():
#             gene = meta['info']['parent']
#             if name != longest_isoforms[gene]['transcript']:
#                 continue
#             else:
#                 filtered[region][name] = meta

#     return filtered

def coding_length(coords):
    return sum([abs(stop-start) for start, stop in coords])

def get_coding_seq(seq, coord_list):
    full_seq = ''
    for coord in sorted(coord_list):
        full_seq += get_subseq(seq, *coord)
    
    return full_seq

def format_output(region_seq, t_dict, verbose=False):
    t_info = t_dict['info']
    t_name = t_info['name']
    strand = t_info['strand']
    span = ':'.join(map(str, t_info['coords']))
    length = str(t_info['length'])
    seq_coords = t_dict['children']
    seq = get_coding_seq(region_seq, seq_coords)
    if strand == '-':
        seq = reverse_complement(seq)
    if TRANSLATE is True:
        seq = translate_seq(seq)
    gene = t_info['parent']
    header_bits = [t_name, gene, strand, span, length]
    if verbose is True:
        coord_string = ','.join(['-'.join(list(map(str, c))) for c in seq_coords])
        header_bits.append(coord_string)
    header = '\t'.join(header_bits)

    return '>{}\n{}'.format(header, seq)

parser = argparse.ArgumentParser(
    description='Extract transcript/coding sequences from '
    'an annotation + genome file combination',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    'genome',
    help='genome file in FASTA format')
parser.add_argument(
    'annotation',
    help='annotation file in GFF[3]/GTF format')
parser.add_argument(
    '-t',
    '--translate',
    action='store_true',
    help='translate the output sequence')
parser.add_argument(
    '-e',
    '--exon',
    action='store_true',
    help='use exons instead of CDS entries to define coding sequence',
)
parser.add_argument(
    '-i',
    '--isoforms',
    help='allow multiple isoforms per gene, instead of only longest',
    action='store_true'
)
parser.add_argument(
    '-c',
    '--coord_based_isoforms',
    help=('detect isoforms by overlapping coordindates in addition to shared '
          'parent (useful for annotations without gene entries)'),
    action='store_true'
)
parser.add_argument(
    '-v',
    '--verbose_headers',
    action='store_true',
    help='include coordinate info in output headers'
)

if len(sys.argv) == 1:
    sys.exit(parser.print_help())

args = parser.parse_args()

ANNOTATION = args.annotation
GENOME = args.genome
TRANSLATE = args.translate
EXON_DEF = args.exon
VERBOSE = args.verbose_headers
ISOFORMS = args.isoforms
FILTER_BY_COORDS = args.coord_based_isoforms

if EXON_DEF:
    child_type = 'exon'
else:
    child_type = 'cds'

child_types = ('exon', 'cds')

transcripts, child_found = get_transcripts(ANNOTATION, child_type)

# ensure there are features we can use in file
if not child_found:
    alt_child = next(c for c in child_types if c != child_type)
    print(
        '[!] No {} entries found in annotation; using {} instead'.format(
            child_type, alt_child), file=sys.stderr)
    child_type = alt_child
    transcripts, child_found = get_transcripts(ANNOTATION, child_type)
    if not child_found:
        sys.exit('[!] No {} entries found in annotation.'.format(child_type))

if ISOFORMS:  # don't filter out extra isoformsy
    transcripts = finalize_transcripts(transcripts)
else:
    transcripts = longest_isoforms(transcripts, FILTER_BY_COORDS)

seq_count = 0
total_regions = len(transcripts)

for region, region_seq in fasta_parse(GENOME):
    if total_regions == 0:  # don't keep looping if we're done
        break
    try:
        for t_name, t_dict in sorted(transcripts[region].items()):
            if not t_dict['children']:
                continue
            print(format_output(region_seq, t_dict, verbose=VERBOSE), flush=True)
            seq_count += 1
        total_regions -= 1
    except KeyError:
        continue
    

print('[#] Extracted {} coding sequences'.format(seq_count), file=sys.stderr)

sys.exit(0)

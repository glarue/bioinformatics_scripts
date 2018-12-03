#!/usr/bin/python3
# the above uses specific Python version; allows script name in top
##!/usr/bin/env python3  
# the above sources Python from $PATH

"""
usage: cdseq.py [-h] [-t] [-e] [-i] [-c] [-v] [-n]
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
  -n, --non_coding      include non-coding (intronic, UTR, etc.) sequence; uses
                        only the coordinates of the transcript itself
                        (default: False)

"""

import sys
import argparse
from collections import defaultdict
from operator import itemgetter
import time
try:
    from biogl import flex_open as open_func
except ModuleNotFoundError:
    open_func = open

class GFFLineInfo(object):
    """
    Takes a gff3/gtf/gff annotation line, and returns available metadata
    about it.

    """
    def __init__(self, line, line_number):
        self.bits = self.__split_on_tabs(line)
        try:
            self.region = self.bits[0]
            try:
                self.start = min(map(int, self.bits[3:5]))
                self.stop = max(map(int, self.bits[3:5]))
            except ValueError:
                self.start = None
                self.stop = None
            self.strand = self.bits[6]
            if self.strand not in ('+', '-'):
                self.strand = '+'
            if self.bits[7] in ['0', '1', '2']:
                self.phase = int(self.bits[7])
            else:
                self.phase = None
            self.infostring = self.bits[8]
            self.feat_type = self.bits[2].lower()
            self.parent = self.get_parent()
            self.grandparent = None
            # try to get grandparent for child types
            if self.feat_type in ('cds', 'exon'):
                self.grandparent = self.__field_match(self.infostring, ['gene_id', 'geneID'], delimiter=';')
            self.name = self.get_ID()
            self.line_number = line_number
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
    def __field_match(infostring, tags, delimiter, tag_order=False):
        if tag_order:
            # check for first match of tags in order
            try:
                tags = [next(t for t in tags if t.lower() in infostring.lower())]
            except StopIteration:
                return None
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
            return og_type

        # Not an obvious type, so search for features of transcripts
        # and genes in infostring to try to infer type

        # check for explicit mention of transcript in ID
        try:
            id_string = next(
                (f for f in delimiter.split(self.infostring)
                if f.startswith("ID")))

            if any(tag in id_string for tag in ('transcript', 'mrna')):
                return 'transcript'
        except StopIteration:
            pass

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
            else:
                return og_type


    def get_ID(self, delimiter=";"):
        """
        Finds the ID of a given annotation file line.

        """
        # first, do it the easy way
        prefix = None
        infostring = self.infostring
        match = self.__field_match(infostring, ["ID="], delimiter)
        if match:
            return match

        # Constrain feature types to simplify indexing
        feat_type = self.feat_type
        if feat_type == "mrna":
            feat_type = "transcript"
        # all get lowercased in the comparison
        # if is no 'ID=', should reference self via others
        gene_tags = ["ID=", "gene_id", "geneId"]
        transcript_tags = ["ID=", "transcriptId", "transcript_ID"]
        tag_selector = {
            "gene": gene_tags,
            "transcript": transcript_tags
        }
        try:
            tags = tag_selector[feat_type]
        except KeyError:
            # get any ID available, prepended with the feature type
            # to keep different features of the same transcript unique
            prefix = self.feat_type
            tags = ['transcriptID', 'transcript_ID', 'gene_ID', 'geneID']

        match = self.__field_match(
            infostring, tags, delimiter, tag_order=True)

        # if nothing matches, return infostring if there's only
        # one tag in it (common for gtf parent features)
        if match is None and infostring.count(";") < 2:
            match = infostring.split(";")[0]

        if prefix:
            match = '{}_{}'.format(prefix, match)

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
            "Parent=", "transcript_ID",
            "transcriptId", "proteinId", "protein_ID"
        ]
        transcript_tags = ["Parent=", "gene_ID", "geneId"]
        gene_tags = ["Parent="]
        tag_selector = {
            "gene": gene_tags,
            "transcript": transcript_tags,
            "exon": child_tags
        }
        try:
            tags = tag_selector[feat_type]
        except KeyError:
            tags = list(set(child_tags + transcript_tags))
        infostring = self.infostring
        match = self.__field_match(infostring, tags, delimiter, tag_order=True)
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
    with open_func(fasta) as f:
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


def translate_seq(
    string, 
    verbosity="single", 
    phase=0,
    codonMap={
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
    }):

    def _get_codons(s, p=0):
        for i in range(0, len(s), 3):
            codon = s[i + p:i + p + 3]
            if codon:  # don't return blank strings
                yield codon

    verbosityD = {"single": 0, "short": 1, "long": 2}
    string = string.replace(" ", "").upper()  # remove spaces if present
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
            amino_acids.append('X')
    return joinChar.join(amino_acids)


def get_transcripts(gff, child_type):
    transcripts = defaultdict(dict)
    feature_info = defaultdict(dict)
    grandparents = defaultdict(dict)
    child_type_found = False
    regions_with_content = set()
    orphans = 0
    with open_func(gff) as annot:
        for ln, line in enumerate(annot):
            feat = None
            try:
                feat = GFFLineInfo(line, ln)
            except TypeError:
                continue
            if feat.feat_type != child_type:
                info_dict = {
                    'name': feat.name,
                    'strand': feat.strand,
                    # use transcript names as genes if no 
                    # parent field in transcript
                    'parent': feat.parent if feat.parent else feat.name,
                    'coords': (feat.start, feat.stop),
                    'region': feat.region,
                    'inferred': False
                }
                if feat.name not in feature_info[feat.region]:
                    feature_info[feat.region][feat.name] = info_dict
                if feat.name in transcripts[feat.region]:  # made by child; update
                    # don't override existing values with Nones
                    updater = {k: v for k, v in info_dict.items() if v is not None}
                    transcripts[feat.region][feat.name]['info'].update(updater)
                else:
                    transcripts[feat.region][feat.name] = {
                        'info': info_dict,
                        'children': set()}
            elif feat.feat_type == child_type:
                if not child_type_found:
                    child_type_found = True
                parent = feat.parent
                if not parent:
                    orphans += 1
                    continue
                start = feat.start
                stop = feat.stop
                region = feat.region
                if parent not in transcripts[region]:
                    try:
                        parent_info = feature_info[region][parent]
                    except KeyError:
                        parent_info = {
                            'strand': feat.strand,
                            'region': region,
                            'name': parent,
                            'inferred': True,
                            'parent': feat.grandparent,
                            'coords': None
                        }
                    transcripts[region][parent] = {
                        'info': parent_info,
                        'children': set()}
                if feat.grandparent and parent not in grandparents[region]:
                        grandparents[region][parent] = feat.grandparent
                    # parent_record = transcripts[region][parent]
                    # existing_gp = parent_record['info']['parent']
                    # if existing_gp != feat.grandparent:
                    #     parent_record['info']['parent'] = feat.grandparent
                tr = transcripts[region][parent]
                tr['children'].add((start, stop))
                regions_with_content.add(region)

    if orphans:
        print('[#] Skipped {} orphan {} features'.format(
              orphans, child_type), file=sys.stderr)

    final_transcripts = {}
    for region, trs in transcripts.items():
        if region not in regions_with_content:
            continue
        elif region not in final_transcripts:
            final_transcripts[region] = defaultdict(dict)
        for name, tr in trs.items():
            if not tr['children']:
                continue
            tr['info']['length'] = coding_length(tr['children'])
            if tr['info']['inferred']:  # need to infer coords from children
                min_coord = min([e[0] for e in tr['children']])
                max_coord = max([e[1] for e in tr['children']])
                tr['info']['coords'] = (min_coord, max_coord)
            # use child grandparent info for transcripts without preexisting
            # parent info (if available)
            if name in grandparents[region] and tr['info']['parent'] == name:
                tr['info']['parent'] = grandparents[region][name]
            final_transcripts[region][name] = tr

    return final_transcripts, child_type_found


def overlap_check(a, b):
    """
    a and b are both (start, stop) coord pairs
    
    """
    lowest = min([a, b], key=lambda x: min(x))
    highest = max([a, b], key=lambda x: min(x))

    return min(highest) <= max(lowest)


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
        # coords only have meaning within a given region
        seen_coords = set()
        # sort by longest transcripts first
        
        # if there are transcripts with the same coding length,
        # sorting may be different each time within those equal-length
        # transcripts, which may change the output if a given transcript's 
        # coord span is larger/smaller
        for name, meta in sorted(transcripts.items(), 
        key=lambda x: (-x[1]['info']['length'], x[1]['info']['name'])):
            gene = meta['info']['parent']
            if gene not in seen_genes:
                if use_coords:
                    # skip those overlapping longer transcripts
                    # even if they have unique gene name
                    coords = meta['info']['coords']
                    if any(overlap_check(coords, c) for c in seen_coords):
                        continue
                    seen_coords.add(coords)
                if gene is not None:
                    seen_genes.add(gene)
                # length = coding_length(meta['children'])
                # meta['info']['length'] = length
                longest_isoforms[region][name] = meta
    
    return longest_isoforms


def finalize_transcripts(transcript_dict):
    """
    Reformats transcripts, and returns a dictionary.
    
    """
    finalized = defaultdict(dict)
    for region, transcripts in transcript_dict.items():
        for name, meta in sorted(transcripts.items()):
            # try:
            #     gene = meta['info']['parent']
            # except:
            #     print(name, meta, file=sys.stderr)
            length = coding_length(meta['children'])
            meta['info']['length'] = length
            finalized[region][name] = meta
    
    return finalized


def coding_length(coords):
    return sum([abs(stop-start) + 1 for start, stop in coords])


def get_coding_seq(seq, coord_list):
    seq_list = []
    for coord in sorted(coord_list):
        seq_list.append(get_subseq(seq, *coord))
    full_seq = ''.join(seq_list)

    return full_seq


def format_output(region_seq, t_dict, non_coding=False, verbose=False):
    t_info = t_dict['info']
    t_name = t_info['name']
    strand = t_info['strand']
    region = t_info['region']
    span = ':'.join(map(str, t_info['coords']))
    # length = str(t_info['length'])
    if non_coding is True:
        seq_coords = [t_info['coords']]
    else:
        seq_coords = t_dict['children']
    seq = get_coding_seq(region_seq, seq_coords)
    if strand == '-':
        seq = reverse_complement(seq)
    if TRANSLATE is True:
        seq = translate_seq(seq)
    length = str(len(seq))
    gene = t_info['parent']
    if gene is None: 
        gene = 'unknown'
    header_bits = [t_name, gene, region, strand, span, length]
    if verbose is True:
        coord_string = ','.join(['-'.join(list(map(str, c))) for c in sorted(seq_coords)])
        header_bits.append(coord_string)
    header = '\t'.join(header_bits)

    return '>{}\n{}'.format(header, seq)
    

def get_runtime(start_time, p=3):
    """
    Takes a start time and optional decimal precision p,
    returns a string of the total run-time until current
    time with appropriate units.

    """
    total = time.time() - start_time  # start with seconds
    divided = total/60.0
    if divided < 2:
        run_time = total
        units = "seconds"
    elif divided < 60:
        run_time = divided
        units = "minutes"
    else:
        run_time = divided/60.0
        units = "hours"
    rounded = round(run_time, p)

    return "{} {}".format(rounded, units)


parser = argparse.ArgumentParser(
    description='Extract transcript/coding sequences from a '
    'genome using an annotation file',
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
parser.add_argument(
    '-n',
    '--non_coding',
    action='store_true',
    help=(
        'include non-coding (intronic, UTR, etc.) sequence; uses only the '
        'coordinates of the transcript itself')
)

if len(sys.argv) == 1:
    sys.exit(parser.print_help())

start = time.time()

args = parser.parse_args()

ANNOTATION = args.annotation
GENOME = args.genome
TRANSLATE = args.translate
EXON_DEF = args.exon
VERBOSE = args.verbose_headers
ISOFORMS = args.isoforms
FILTER_BY_COORDS = args.coord_based_isoforms
NON_CODING = args.non_coding

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

if not ISOFORMS:
    transcripts = longest_isoforms(transcripts, FILTER_BY_COORDS)

seq_count = 0
total_regions = len(transcripts)

for region, region_seq in fasta_parse(GENOME):
    if region not in transcripts:
        continue
    region_dict = transcripts[region]
    for t_name, t_dict in sorted(region_dict.items()):
        if not t_dict['children']:
            continue
        print(format_output(
            region_seq, t_dict, non_coding=NON_CODING, verbose=VERBOSE), 
            flush=True)
        seq_count += 1
    total_regions -= 1
    if total_regions == 0:  # don't keep looping if we're done
        break

runtime = get_runtime(start)

print('[#] Extracted {} coding sequences in {}'.format(
        seq_count, runtime), file=sys.stderr)

sys.exit(0)

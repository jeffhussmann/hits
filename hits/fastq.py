''' Utilities for dealing with fastq files. '''

import array
import functools
import gzip
import heapq
import os
import s3fs
import tempfile
from pathlib import Path
from itertools import chain
from collections import namedtuple

import numpy as np

from . import utilities
from .fastq_cython import *
from .utilities import identity, base_order, group_by

# SANGER_OFFSET = 33 is imported from fastq_cython
SOLEXA_OFFSET = 64
MAX_QUAL = 93
MAX_EXPECTED_QUAL = 42

MIN_SOLEXA_ORD = SOLEXA_OFFSET - 5
MAX_SANGER_ORD = SANGER_OFFSET + 41

def decode_sanger(qual):
    ''' Converts a string of sanger-encoded quals to a list of integers. '''
    return [ord(q) - SANGER_OFFSET for q in qual]

def decode_sanger_to_array(qual):
    ''' Converts a string of sanger-encoded quals to a (python, not numpy) array of unsigned chars
    suitable for assignment to a pysam.AlignedSegment.query_qualities '''
    q_array = array.array('B', qual.encode())
    subtract_sanger_offset(q_array)
    return q_array

def decode_solexa(qual):
    ''' Converts a string of solexa-encoded quals to a list of integers. '''
    return [ord(q) - SOLEXA_OFFSET for q in qual]

def encode_sanger(ints):
    ''' Converts a list of integer quals to a sanger-encoded string. '''
    return ''.join(chr(i + SANGER_OFFSET) for i in ints)

def encode_solexa(ints):
    ''' Converts a list of integer quals to a solexa-encoded string. '''
    return ''.join(chr(i + SOLEXA_OFFSET) for i in ints)

_solexa_to_sanger_table = {}
for q in range(ord('!') - SOLEXA_OFFSET, MAX_EXPECTED_QUAL + 1):
    # Old solexa encoding was -10 log(p / (1 - p)), which could be negative.
    # Character encodings of negative values cause problems, and we don't really
    # care about fine distinctions in low quality scores, so just set to a
    # minimum of zero.
    # Some files assign ! to N's, so this range needs to go down to ord('!').
    _solexa_to_sanger_table[chr(q + SOLEXA_OFFSET)] = chr(max(q, 0) + SANGER_OFFSET)

def solexa_to_sanger(qual):
    return ''.join(_solexa_to_sanger_table[c] for c in qual)

# If a qname for a paired mapping ends in '/1', '/2', or '/3', bowtie2 chops off
# the last two characters of the qname. If qual strings of trimmed portions of
# reads are to be put in qnames, '/' needs to be downgraded to
# chr(ord('/') - 1). 
# '_' is used as a field separator in annotations, so similarly needs to be
# downgraded to chr(ord('_') - 1).
_chars_to_sanitize = '/_'
_sanitized_chars = ''.join(chr(ord(c) - 1) for c in _chars_to_sanitize)
try:
    _sanitize_table = str.maketrans(_chars_to_sanitize, _sanitized_chars)
    period_to_N = str.maketrans('.', 'N')
except AttributeError:
    import string
    _sanitize_table = string.maketrans(_chars_to_sanitize, _sanitized_chars)
    period_to_N = string.maketrans('.', 'N')

def sanitize_qual(qual):
    sanitized = qual.translate(_sanitize_table)
    return sanitized

def quality_and_complexity(reads_iter, max_read_length, alignments=False, min_q=0):
    stats = {
        'q': np.zeros((max_read_length, MAX_EXPECTED_QUAL + 1), int),
        'c': np.zeros((max_read_length, 256), int),
        'c_above_min_q': np.zeros((max_read_length, 256), int),
        'average_q': np.zeros((max_read_length, 256), int),
    }

    if isinstance(reads_iter, (str, Path)):
        reads_iter = reads(reads_iter)

    for read in reads_iter:
        if alignments:
            process_Alignment(read.query_sequence.encode(),
                              read.query_qualities,
                              stats['q'],
                              stats['average_q'],
                              stats['c'],
                              stats['c_above_min_q'],
                              min_q,
                             )
        else:
            process_read(read.seq.encode(), read.qual.encode(),
                         stats['q'],
                         stats['average_q'],
                         stats['c'],
                         stats['c_above_min_q'],
                         min_q,
                        )
        
    # To avoid a lookup at every single base, base-specific arrays are 2*max_read_length x 256.
    # This pulls out only the columns corresponding to possible base
    # identities. 
    for k in ['c', 'c_above_min_q', 'average_q']:
        stats[k] = np.vstack([stats[k].T[ord(b)] for b in base_order]).T

    stats['average_q'] = stats['average_q'] / np.maximum(1, stats['c'])
    
    return stats

def quality_and_complexity_paired(read_pairs, max_read_length, min_q=0):
    stats = {}
    for which in ['R1', 'R2']:
        stats[which] = {
            'q': np.zeros((max_read_length, MAX_EXPECTED_QUAL + 1), int),
            'c': np.zeros((max_read_length, 256), int),
            'c_above_min_q': np.zeros((max_read_length, 256), int),
            'average_q': np.zeros((max_read_length, 256), int),
        }
    
    joint_average_q_distribution = np.zeros((MAX_EXPECTED_QUAL + 1, MAX_EXPECTED_QUAL + 1), int)
    
    for R1, R2 in read_pairs:
        process_read(R1.seq.encode(),
                     R1.qual.encode(),
                     stats['R1']['q'],
                     stats['R1']['average_q'],
                     stats['R1']['c'],
                     stats['R1']['c_above_min_q'],
                     min_q,
                    )
        process_read(R2.seq.encode(),
                     R2.qual.encode(),
                     stats['R2']['q'],
                     stats['R2']['average_q'],
                     stats['R2']['c'],
                     stats['R2']['c_above_min_q'],
                     min_q,
                    )
        
    # See comment in quality_and_complexity above. 
    for which in ['R1', 'R2']:
        for k in ['c', 'c_above_min_q', 'average_q']:
            stats[which][k] = np.vstack([stats[which][k].T[ord(b)] for b in base_order]).T
    
    return stats

def get_line_groups(line_source):
    if isinstance(line_source, Path):
        line_source = str(line_source)

    if isinstance(line_source, str):
        s3_prefix = 's3://'
        if line_source.startswith(s3_prefix):
            s3 = s3fs.S3FileSystem()
            def opener(s3_path):
                if s3_path.endswith('.gz'):
                    with s3.open(s3_path[len(s3_prefix):], 'rb') as s3_fh:
                        with gzip.open(s3_fh, mode='rt') as gz_fh:
                            yield from gz_fh
                else:
                    with s3.open(s3_path[len(s3_prefix):]) as s3_fh:
                        yield from s3_fh

        elif line_source.endswith('.gz'):
            opener = functools.partial(gzip.open, mode='rt')
        else:
            opener = open
    else:
        opener = iter

    lines = opener(line_source)
    groups = zip(*[lines]*4)
    return groups

class Read(object):
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual
        
    def __str__(self):
        return f'@{self.name}\n{self.seq}\n+\n{self.qual}\n'
    
    def reverse_complement(self):
        return Read(self.name,
                    utilities.reverse_complement(self.seq),
                    self.qual[::-1],
                   )
    
    def __getitem__(self, sl):
        return Read(self.name, self.seq[sl], self.qual[sl])
    
    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.seq)

    def __add__(self, other):
        return Read(self.name, self.seq + other.seq, self.qual + other.qual)

    @property
    def query_sequence(self):
        return self.seq

    @property
    def query_name(self):
        return self.name

    @property
    def qname(self):
        return self.name
    
    @utilities.memoized_property
    def query_qualities(self):
        return array.array('B', decode_sanger(self.qual))

    @utilities.memoized_property
    def Q30_fraction(self):
        return np.mean(np.array(self.query_qualities) >= 30)

    @utilities.memoized_property
    def Q93_fraction(self):
        return np.mean(np.array(self.query_qualities) == 93)

    @utilities.memoized_property
    def mean_Q(self):
        return np.mean(self.query_qualities)
    
def line_group_to_read(line_group, name_standardizer=identity, qual_convertor=identity):
    name_line, seq_line, _, qual_line = line_group

    name = name_standardizer(name_line.rstrip().lstrip('@'))

    seq = seq_line.strip()
    seq = seq.translate(period_to_N)

    qual = qual_convertor(qual_line.strip())

    return Read(name, seq, qual)

def get_read_length(file_name):
    read = next(reads(file_name))
    return len(read)

def reads(file_name, standardize_names=False, ensure_sanger_encoding=False, up_to_space=False):
    ''' Yields Read's from a file name or line iterator.
        If standardize_names == True, infers read name structure and
        standardizes read names.
        If ensure_sanger_encoding == True, detects the quality score encoding
        and converts to sanger if necessary.
    '''
    if isinstance(file_name, list):
        args = [standardize_names, ensure_sanger_encoding, up_to_space]
        files = [reads(fn, *args) for fn in file_name]
        return chain.from_iterable(files)

    line_groups = get_line_groups(file_name)

    if standardize_names:
        name_standardizer, line_groups = detect_structure(line_groups)
    elif up_to_space:
        name_standardizer = lambda n: n.split()[0]
    else:
        name_standardizer = identity

    if ensure_sanger_encoding:
        qual_convertor, line_groups = detect_encoding(line_groups)
    else:
        qual_convertor = identity
    
    reads_ = (
        line_group_to_read(line_group, name_standardizer, qual_convertor)
        for line_group in line_groups
    )

    return reads_

def reverse_complement_reads(file_name, **kwargs):
    for read in reads(file_name, **kwargs):
        yield read.reverse_complement()

def detect_structure(line_groups):
    ''' Look at the first read to figure out the read name structure. '''
    try:
        first_group = next(line_groups)
        first_read = line_group_to_read(first_group)
        name_standardizer = get_read_name_standardizer(first_read.name)
        line_groups = chain([first_group], line_groups)
    except StopIteration:
        name_standardizer = identity
        # Note: line_groups is now an empty iterator.
    
    return name_standardizer, line_groups

def detect_encoding(line_groups):
    groups_examined = []
    try:
        for line_group in line_groups:
            groups_examined.append(line_group)
            read = line_group_to_read(line_group)
            ords = [ord(q) for q in read.qual]
            if min(ords) < MIN_SOLEXA_ORD:
                encoding = 'SANGER'
                break
            if max(ords) > MAX_SANGER_ORD:
                encoding = 'SOLEXA'
                break
    except StopIteration:
        encoding = 'SANGER'
        # Note: line_groups is now an empty iterator.

    if encoding == 'SOLEXA':
        qual_convertor = solexa_to_sanger
    else:
        qual_convertor = identity

    line_groups = chain(groups_examined, line_groups)
    
    return qual_convertor, line_groups

def unambiguous_sanger_Q40(length):
    ''' Return a qual str of specified length with first character
        equal to the highest sanger-encoding-specific value and
        the remaining character equal to sanger Q40.
    '''
    first_char = chr(MIN_SOLEXA_ORD - 1) 
    remaining_chars = encode_sanger([40]*(length - 1))
    return first_char + remaining_chars

def read_pairs(R1_file_name, R2_file_name, **kwargs):
    R1_reads = reads(R1_file_name, **kwargs)
    R2_reads = reads(R2_file_name, **kwargs)
    return zip(R1_reads, R2_reads)

def read_pairs_interleaved(lines, **kwargs):
    interleaved_reads = reads(lines, **kwargs)
    grouped = group_by(interleaved_reads, key=lambda r: get_pair_name(r.name))
    for pair_name, group in grouped:
        if len(group) != 2:
            raise ValueError(group)
        R1, R2 = group
        R1_renamed = Read(pair_name, R1.seq, R1.qual)
        R2_renamed = Read(pair_name, R2.seq, R2.qual)
        yield R1_renamed, R2_renamed

quartet_order = ['I1', 'I2', 'R1', 'R2']
Quartet = namedtuple('Quartet', quartet_order)

def read_quartets(fns, **kwargs):
    if isinstance(fns, dict):
        fns = [fns[which] for which in quartet_order]

    all_reads = [reads(fn, **kwargs) for fn in fns]

    return (Quartet(*reads) for reads in zip(*all_reads))

def get_read_name_parser(read_name):
    if read_name.startswith('test') or read_name.startswith('simulated'):
        # Simulated data sometimes needs read names to contain information
        # and can't use standard Illumina-formatted read names.
        parser = None

    elif read_name.startswith('SRR'):
        if len(read_name.split('.')) == 2:
            parser = parse_SRA_read_name
        elif len(read_name.split('.')) == 3:
            parser = parse_paired_SRA_read_name

    elif read_name.startswith('ERR') or read_name.startswith('DRR'):
        parser = parse_ERR_read_name

    else:
        num_words = len(read_name.split())

        if num_words > 1:
            location_info, member_info = read_name.split()[:2]
            fields = location_info.split(':')

            if len(fields) == 8 and any(b in fields[-1] for b in 'TCAG'):
                parser = parse_illumina_read_name_with_UMI
            else:
                parser = parse_new_illumina_read_name

        elif num_words == 1:
            if '#' in read_name:
                parser = parse_old_illumina_read_name
            elif '/' in read_name:
                parser = parse_unindexed_old_illumina_read_name
            else:
                parser = parse_standardized_name

        else:
            raise ValueError(f'read name format not recognized - {read_name}')

    return parser

def get_read_name_standardizer(read_name):
    ''' Looks at structure of read_name to determine the appropriate read name
        standardizer.
    '''
    parser = get_read_name_parser(read_name)

    if parser == parse_SRA_read_name or parser == parse_ERR_read_name:
        standardize = templates['SRA'].format
        def standardizer(read_name):
            accession, number = parser(read_name)
            standardized = standardize(accession, number)
            return standardized
            
    elif parser == parse_paired_SRA_read_name:
        standardize = templates['paired_SRA'].format
        def standardizer(read_name):
            accession, number, member = parser(read_name)
            standardized = standardize(accession, number, member)
            return standardized

    elif parser == parse_illumina_read_name_with_UMI:
        standardize = templates['UMI'].format
        def standardizer(read_name):
            lane, tile, x, y, member, index, UMI = parser(read_name)
            standardized = standardize(lane, tile, x, y, UMI)
            return standardized

    elif parser is not None:
        standardize = templates['default'].format
        def standardizer(read_name):
            lane, tile, x, y, member, index = parser(read_name)
            standardized = standardize(lane, tile, x, y)
            return standardized

    else:
        standardizer = identity

    return standardizer

def standardize(template, *args):
    return templates[template].format(*args)

templates = {
    'default': '{0:0>2.2s}:{1:0>5.5s}:{2:0>6.6s}:{3:0>6.6s}',
    'SRA': '{0:0>9.9s}:{1:0>10.10s}',
    'paired_SRA': '{0:0>9.9s}:{1:0>10.10s}:{2:0>1.1s}',
}
templates['UMI'] = templates['default'] + ':{4}'

def parse_illumina_read_name_with_UMI(read_name):
    location_info, member_info = read_name.split()[:2]
    fields = location_info.split(':')

    if len(fields) == 8 and any(b in fields[-1] for b in 'TCAG'):
        UMI = fields[-1]
        fields = fields[:-1]
    else:
        UMI = ''

    lane, tile, x, y = fields[-4:]

    member, _, _, index = member_info.split(':')
    
    return lane, tile, x, y, member, index, UMI

def parse_new_illumina_read_name(read_name):
    location_info, member_info = read_name.split()[:2]

    lane, tile, x, y = location_info.split(':')[-4:]

    member, _, _, index = member_info.split(':')
    
    return lane, tile, x, y, member, index

def parse_old_illumina_read_name(read_name):
    location_info, member_info = read_name.split('#')
    lane, tile, x, y = location_info.split(':')[-4:]
    index, member = member_info.split('/')
    return lane, tile, x, y, member, index

def parse_unindexed_old_illumina_read_name(read_name):
    location_info, member = read_name.split('/')
    lane, tile, x, y = location_info.split(':')[-4:]
    return lane, tile, x, y, member, ''

def parse_standardized_name(read_name):
    lane, tile, x, y, member = read_name.split(':')
    return lane, tile, x, y, member, ''

def parse_SRA_read_name(read_name):
    accession, number = read_name.split()[0].split('.')
    # Remove the leading 'SRR'
    accession = accession[3:]
    return accession, number

def parse_paired_SRA_read_name(read_name):
    accession, number, member = read_name.split('.')
    # Remove the leading 'SRR'
    accession = accession[3:]
    return accession, number, member

def parse_ERR_read_name(read_name):
    accession, number = read_name.split()[0].split('.')
    # Remove the leading 'SRR'
    accession = accession[3:]
    return accession, number

def coordinates_from_standardized(standardized):
    coordinates = standardized.split(':')[:-1]
    return coordinates

def get_pair_name(standardized):
    ''' Returns the part of a standardized name that is common to both members of
        the pair it came from.
    '''
    coordinates = coordinates_from_standardized(standardized)
    pair_name = ':'.join(coordinates)
    return pair_name

class ExternalSorter:
    def __init__(self, read_iter, sorted_fn, chunk_size=int(1e5), sort_key=lambda read: read.query_name):
        self.read_iter = read_iter
        self.sorted_fn = sorted_fn
        self.chunk_size = chunk_size
        self.sort_key = sort_key

        self.chunk_fns = []

    def sort_and_write_chunk(self, chunk):
        sorted_chunk = sorted(chunk, key=self.sort_key)
        
        with tempfile.NamedTemporaryFile(prefix=f'{self.sorted_fn}.', mode='w', delete=False) as chunk_file:
            self.chunk_fns.append(chunk_file.name)

            for read in sorted_chunk:
                chunk_file.write(str(read))

    def merge_sorted_chunks(self):
        sorted_chunks = (reads(fn) for fn in self.chunk_fns)
        yield from heapq.merge(*sorted_chunks, key=self.sort_key)

        for chunk_fn in self.chunk_fns:
            os.remove(chunk_fn)

    def sort(self):
        for chunk in utilities.chunks(self.read_iter, self.chunk_size):
            self.sort_and_write_chunk(chunk)

        with gzip.open(self.sorted_fn, 'wt', compresslevel=1) as fh:
            for read in self.merge_sorted_chunks():
                fh.write(str(read))
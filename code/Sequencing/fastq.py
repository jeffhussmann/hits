''' Utilities for dealing with fastq files. '''

from itertools import izip, chain
from collections import namedtuple
from .fastq_cython import *
from .utilities import identity
import numpy as np

# Note that SANGER_OFFSET and SOLEXA_OFFSET are imported from fastq_cython

def decode_sanger(qual):
    ''' Converts a string of sanger-encoded quals to a list of integers. '''
    return [ord(q) - SANGER_OFFSET for q in qual]

def decode_solexa(qual):
    ''' Converts a string of solexa-encoded quals to a list of integers. '''
    return [ord(q) - SOLEXA_OFFSET for q in qual]

def encode_sanger(ints):
    ''' Converts a list of integer quals to a sanger-encoded string. '''
    return ''.join(chr(i + SANGER_OFFSET) for i in ints)

def solexa_to_sanger(qual):
    ''' Converts a string of solexa-encoded quals to sanger encoding. '''
    ints = decode_solexa(qual)
    # Old solexa encoding was -10 log(p / (1 - p)), which could be negative.
    # Character encodings of negative values cause problems, and we don't really
    # care about fine distinctions in low quality scores, so just set to a
    # minimum of zero.
    nonnegative = np.maximum(ints, 0)
    return encode_sanger(nonnegative)

def get_line_groups(line_source):
    if type(line_source) == str:
        # line_source is a file name.
        lines = open(line_source)
    else:
        # line_source is an open file.
        lines = iter(line_source)
    return izip(*[lines]*4)

Read = namedtuple('Read', ['name', 'seq', 'qual'])

def line_group_to_read(line_group, name_standardizer=identity, qual_convertor=identity):
    name_line, seq_line, _, qual_line = line_group
    name = name_standardizer(name_line.rstrip().lstrip('@'))
    seq = seq_line.strip()
    qual = qual_convertor(qual_line.strip())
    read = Read(name, seq, qual)
    return read

def reads(file_name, standardize_names=False, ensure_sanger_encoding=False):
    ''' Yields Read's from a file name or line iterator.
        If standardize_names == True, infers read name structure and
        standardizes read names.
        If ensure_sanger_encoding == True, detects the quality score encoding
        and converts to sanger if necessary.
    '''
    line_groups = get_line_groups(file_name)

    if standardize_names:
        name_standardizer, line_groups = detect_structure(line_groups)
    else:
        name_standardizer = identity

    if ensure_sanger_encoding:
        qual_convertor, line_groups = detect_encoding(line_groups)
    else:
        qual_convertor = identity
    
    reads = (line_group_to_read(line_group, name_standardizer, qual_convertor)
             for line_group in line_groups)

    return reads

def detect_structure(line_groups):
    ''' Look at the first read to figure out the read name structure. '''
    try:
        first_group = line_groups.next()
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
            if min(ords) < SOLEXA_OFFSET - 5:
                encoding = 'SANGER'
                break
            if max(ords) > SANGER_OFFSET + 41:
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

def read_pairs(R1_file_name, R2_file_name, **kwargs):
    R1_reads = reads(R1_file_name, **kwargs)
    R2_reads = reads(R2_file_name, **kwargs)
    return izip(R1_reads, R2_reads)

make_record = '@{0}\n{1}\n+\n{2}\n'.format

def get_read_name_parser(read_name):
    if read_name.startswith('test'):
        # Simulated data sometimes needs read names to contain information
        # and can't use standard Illumina-formatted read names.
        parser = None
    elif read_name.startswith('SRR'):
        if len(read_name.split('.')) == 2:
            parser = parse_SRA_read_name
        elif len(read_name.split('.')) == 3:
            parser = parse_paired_SRA_read_name
    else:
        num_words = len(read_name.split())
        if num_words == 2:
            parser = parse_new_illumina_read_name
        elif num_words == 1:
            if '#' in read_name:
                parser = parse_old_illumina_read_name
            else:
                parser = parse_standardized_name
        else:
            raise ValueError('read name format not recognized - {}'.format(read_name))
    return parser

def get_read_name_standardizer(read_name):
    ''' Looks at structure of read_name to determine the appropriate read name
        standardizer.
    '''
    parser = get_read_name_parser(read_name)
    if parser == parse_SRA_read_name:
        def standardizer(read_name):
            accession, number = parser(read_name)
            standardized = _standardize_SRA(accession, number)
            return standardized
    elif parser == parse_paired_SRA_read_name:
        def standardizer(read_name):
            accession, number, member = parser(read_name)
            standardized = _standardize_paired_SRA(accession, number, member)
            return standardized
    elif parser:
        def standardizer(read_name):
            lane, tile, x, y, member, index = parser(read_name)
            standardized = _standardize(lane, tile, x, y, member)
            return standardized
    else:
        standardizer = identity

    return standardizer

_standardize = '{0:0>2.2s}:{1:0>5.5s}:{2:0>6.6s}:{3:0>6.6s}:{4:0>1.1s}'.format
_standardize_SRA = '{0:0>9.9s}:{1:0>10.10s}'.format
_standardize_paired_SRA = '{0:0>9.9s}:{1:0>10.10s}:{2:0>1.1s}'.format

def parse_new_illumina_read_name(read_name):
    location_info, member_info = read_name.split()
    lane, tile, x, y = location_info.split(':')[-4:]
    member, _, _, index = member_info.split(':')
    return lane, tile, x, y, member, index

def parse_old_illumina_read_name(read_name):
    location_info, member_info = read_name.split('#')
    lane, tile, x, y = location_info.split(':')[-4:]
    index, member = member_info.split('/')
    return lane, tile, x, y, member, index

def parse_standardized_name(read_name):
    lane, tile, x, y, member = read_name.split(':')
    return lane, tile, x, y, member, ''

def parse_SRA_read_name(read_name):
    accession, number = read_name.split('.')
    # Remove the leading 'SRR'
    accession = accession[3:]
    return accession, number

def parse_paired_SRA_read_name(read_name):
    accession, number, member = read_name.split('.')
    # Remove the leading 'SRR'
    accession = accession[3:]
    return accession, number, member

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

''' Utilities for dealing with sam files. '''

import utilities
import re
import subprocess
from collections import Counter
from itertools import izip
import os
import shutil
import external_sort
import pysam
import fastq

BAM_CMATCH = 0 # M
BAM_CINS = 1 # I
BAM_CDEL = 2 # D
BAM_CREF_SKIP = 3 # N
BAM_CSOFT_CLIP = 4 # S
BAM_CHARD_CLIP = 5 # H
BAM_CPAD = 6 # P
BAM_CEQUAL = 7 # =
BAM_CDIFF = 8 # X

_unmapped_template = '{0}\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n'.format

def unmapped_aligned_read(qname):
    aligned_read = pysam.AlignedRead()
    aligned_read.qname = qname
    aligned_read.flag = 0x10
    aligned_read.rname = -1
    aligned_read.pos = -1
    aligned_read.mapq = 0
    aligned_read.cigar = None
    aligned_read.rnext = -1
    aligned_read.pnext = -1
    aligned_read.tlen = 0
    aligned_read.seq = '*'
    aligned_read.qual = '*'
    return aligned_read

def splice_in_name(line, new_name):
    return '\t'.join([new_name] + line.split('\t')[1:])

def count_header_lines(sam_file_name):
    ''' Returns the total number of header lines in sam_file_name. '''
    count = 0
    with open(sam_file_name) as sam_file:
        for line in sam_file:
            if line.startswith('@'):
                count += 1
            else:
                break
    return count

def get_header_lines(sam_file_name):
    ''' Returns a list of all header lines in sam_file_name. '''
    header_lines = []
    with open(sam_file_name) as sam_file:
        for line in sam_file:
            if not line.startswith('@'):
                break
            else:
                header_lines.append(line)
    return header_lines

def extract_seq_lengths(lines):
    seq_lengths = {}
    for line in lines:
        if line.startswith('@SQ'):
            split_line = line.split('\t')
            
            SN_tag = split_line[1]
            seq_name = ':'.join(SN_tag.split(':')[1:])
            
            LN_tag = split_line[2]
            seq_length = int(LN_tag.split(':')[1])
            
            seq_lengths[seq_name] = seq_length
    
    return seq_lengths

def get_sq_lines(sam_file_name):
    ''' Returns a list of the SQ header lines in sam_file_name. '''
    sq_lines = []
    with open(sam_file_name) as sam_file:
        for line in sam_file:
            if not line.startswith('@'):
                break
            if line.startswith('@SQ'):
                sq_lines.append(line)
    return sq_lines

def open_to_reads(sam_file_name):
    ''' Returns an open file that has been advanced to the first read line in
        sam_file_name (i.e. past all header lines.)
    '''
    header_line_count = count_header_lines(sam_file_name)
    sam_file = open(sam_file_name)
    for i in range(header_line_count):
        sam_file.readline()
    return sam_file

def mapped_reads(sam_file_name):
    ''' Returns an iterator over the lines of sam_file_name that correspond to
        mapped reads.
    '''
    def is_mapped(line):
        parsed = parse_line(line)
        return parsed['mapped']
    
    lines = open_to_reads(sam_file_name)
    mapped_lines = (line for line in lines if is_mapped(line))
    return mapped_lines
    
def parse_line(line):
    ''' Returns a dictionary of the information in a SAM line.
    '''
    fields = line.rstrip().split('\t')
    
    flags = int(fields[1])
    parsed_line = {'mapped':    not (flags & 0x4),
                   'QNAME':     fields[0],
                   'strand':    '-' if (flags & 0x10) else '+',
                   'RNAME':     fields[2],
                   'POS':       int(fields[3]) - 1, # SAM is 1-indexed
                   'MAPQ':      int(fields[4]),
                   'CIGAR':     fields[5],
                   'TLEN':      int(fields[8]),
                   'SEQ':       fields[9],
                   'QUAL':      fields[10],  
                  }
    
    # There are 11 mandatory fields, followed by optional tags.
    for field in fields[11:]:
        name, data_type = field.split(':')[:2]
        value = ':'.join(field.split(':')[2:])
        
        # BWA-specific XA field
        if name == 'XA':
            entries = value.rstrip(';').split(';')
            parsed_line['XA'] = []
            for entry in entries:
                ref_seq_name, strand_and_position, cigar, nm = entry.split(',')
                # first character of strand_and_position is [+-], rest is position
                strand = strand_and_position[0]
                position = int(strand_and_position[1:]) - 1 # SAM is 1-indexed
                xa_dict = {'QNAME':     fields[0],
                           'RNAME':     ref_seq_name,
                           'strand':    strand,
                           'POS':       position,
                           'CIGAR':     cigar,
                           'NM':        int(nm),
                          }
                parsed_line['XA'].append(xa_dict)
        else:
            if data_type == 'i':
                value = int(value)
            elif data_type == 'f':
                value = float(value)
            parsed_line[name] = value

    return parsed_line

cigar_block = re.compile(r'(\d+)([MIDNSHP=X])')

md_number = re.compile(r'[0-9]+')
md_text = re.compile(r'[A-Z]+')

def cigar_string_to_blocks(cigar_string):
    """ Decomposes a CIGAR string into a list of its operations. """
    return [(int(l), k) for l, k in re.findall(cigar_block, cigar_string)]

def total_reference_nucs(cigar_blocks):
    """ Returns the number of nucleotides from the reference involved in a list of CIGAR operations. """
    return sum(l for l, k in cigar_blocks if k in ['M', '=', 'X', 'D'])

def contains_indel(parsed_line):
    cigar_blocks = cigar_string_to_blocks(parsed_line['CIGAR'])
    kinds = [k for l, k in cigar_blocks]
    return ('I' in kinds or 'D' in kinds)

def contains_indel_pysam(read):
    kinds = [k for k, l in read.cigar]
    return (1 in kinds or 2 in kinds)

def contains_splicing_pysam(read):
    kinds = [k for k, l in read.cigar]
    return (3 in kinds)

def contains_soft_clipping(parsed_line):
    cigar_blocks = cigar_string_to_blocks(parsed_line['CIGAR'])
    kinds = [k for k, l in cigar_blocks]
    return ('S' in kinds)

def contains_soft_clipping_pysam(read):
    kinds = [k for k, l in read.cigar]
    return (BAM_CSOFT_CLIP in kinds)

def cigar_blocks_to_string(cigar_blocks):
    """ Builds a CIGAR string out of a corresponding list of operations. """
    return ''.join([str(l) + k for l, k in cigar_blocks])

def alignment_to_cigar_blocks(ref_aligned, read_aligned):
    """ Builds a list of CIGAR operations from an alignment. """
    expanded_sequence = []
    for ref_char, read_char in zip(ref_aligned, read_aligned):
        if ref_char == '-':
            expanded_sequence.append('I')
        elif read_char == '-':
            expanded_sequence.append('D')
        elif ref_char == read_char:
            #expanded_sequence.append('=')
            expanded_sequence.append('M')
        else:
            #expanded_sequence.append('X')
            expanded_sequence.append('M')
    sequence, counts = utilities.decompose_homopolymer_sequence(expanded_sequence)
    return [[count, char] for char, count in zip(sequence, counts)]

def truncate_cigar_blocks(cigar_blocks, truncated_length):
    ''' Given pysam-style cigar_blocks, truncates the blocks to explain
        truncated_length read bases.
    '''
    bases_so_far = 0
    truncated_blocks = []

    consuming_operation = set([0, 1, 7, 8])
    for operation, length in cigar_blocks:
        if bases_so_far == truncated_length:
            break
        if operation in consuming_operation:
            length_to_use = min(truncated_length - bases_so_far, length)
            bases_so_far += length_to_use
        else:
            length_to_use = length
        truncated_blocks.append((operation, length_to_use))
    return truncated_blocks

def alignment_to_cigar_string(ref_aligned, read_aligned):
    """ Builds a CIGAR string from an alignment. """
    cigar_blocks = alignment_to_cigar_blocks(ref_aligned, read_aligned)
    return cigar_blocks_to_string(cigar_blocks)

def md_string_to_ops_string(md_string):
    ''' Converts an MD string into a list of operations for supplying reference
        characters, either '=' if equal to the read, or any other char if equal
        to that char.
    '''
    md_string = md_string.translate(None, '^')
    
    match_lengths = map(int, re.findall(md_number, md_string))
    text_blocks = re.findall(md_text, md_string)
    
    ops_string = '='*match_lengths[0]
    for text_block, match_length in zip(text_blocks, match_lengths[1:]):
        ops_string += text_block
        ops_string += '='*match_length
    
    return ops_string

def produce_alignment(parsed_sam_line, from_pysam=False):
    ''' Returns a list of (ref_char, read_char, qual_char, ref_pos, read_pos)
        tuples.
    '''
    if from_pysam:
        read_seq = parsed_sam_line.seq
        read_quals = parsed_sam_line.qual
        MD_string = dict(parsed_sam_line.tags)['MD']
        CIGAR_string = parsed_sam_line.cigarstring
        POS = parsed_sam_line.pos
    else:
        read_seq = parsed_sam_line['SEQ']
        read_quals = parsed_sam_line['QUAL']
        MD_string = parsed_sam_line['MD']
        CIGAR_string = parsed_sam_line['CIGAR']
        POS = parsed_sam_line['POS']
    
    columns = []
    read_seq_iter = iter(read_seq)
    read_quals_iter = iter(fastq.decode_sanger(read_quals))

    ref_ops = iter(md_string_to_ops_string(MD_string))
    cigar_blocks = cigar_string_to_blocks(CIGAR_string)

    ref_pos = POS
    read_pos = 0
    for (length, kind) in cigar_blocks:
        if kind == 'M' or kind == '=' or kind == 'X':
            for i in range(length):
                read_char = read_seq_iter.next()
                ref_op = ref_ops.next()
                if ref_op == '=':
                    ref_char = read_char
                else:
                    ref_char = ref_op
                qual = read_quals_iter.next()
                column = (ref_char, read_char, qual, ref_pos, read_pos)
                columns.append(column)
                ref_pos += 1
                read_pos += 1
        elif kind == 'D':
            # Deletion results in gap in read
            for i in range(length):
                read_char = '-'
                ref_char = ref_ops.next()
                qual = 0 
                column = (ref_char, read_char, qual, ref_pos, read_pos)
                columns.append(column)
                ref_pos += 1
        elif kind == 'I' or kind == 'S':
            # Insertion or soft-clipping results in gap in ref
            for i in range(length):
                read_char = read_seq_iter.next()
                ref_char = '-'
                qual = read_quals_iter.next()
                column = (ref_char, read_char, qual, ref_pos, read_pos)
                columns.append(column)
                read_pos += 1
    
    return columns

def alignment_to_MD_string(ref_aligned, read_aligned):
    ''' Produce an MD string from an alignment. '''
    # Mark all blocks of matching with numbers, all deletion bases with '^*0', and all mismatch bases.
    MD_list = []
    current_match_length = 0
    for ref_char, read_char in izip(ref_aligned, read_aligned):
        if ref_char == read_char:
            current_match_length += 1
        elif ref_char != '-':
            if current_match_length > 0:
                MD_list.append(current_match_length)
                current_match_length = 0

            if read_char == '-':
                MD_list.append('^{0}'.format(ref_char))
                MD_list.append(0)
            else:
                MD_list.append(ref_char)
    
    if current_match_length > 0:
        MD_list.append(current_match_length)

    # Remove all zeros that aren't a deletion followed by a mismatch
    reduced_MD_list = []
    for i in range(len(MD_list)):
        if isinstance(MD_list[i], int):
            if MD_list[i] > 0:
                reduced_MD_list.append(MD_list[i])
            elif 0 < i < len(MD_list) - 1:
                if isinstance(MD_list[i - 1], str) and isinstance(MD_list[i + 1], str) and MD_list[i - 1][0] == '^' and MD_list[i + 1][0] != '^':
                    reduced_MD_list.append(MD_list[i])
        else:
            reduced_MD_list.append(MD_list[i])
    
    # Collapse all deletions.
    collapsed_MD_list = [reduced_MD_list[0]]
    for i in range(1, len(reduced_MD_list)):
        if isinstance(collapsed_MD_list[-1], str) and collapsed_MD_list[-1][0] == '^' and \
           isinstance(reduced_MD_list[i], str) and reduced_MD_list[i][0] == '^':
            
            collapsed_MD_list[-1] += reduced_MD_list[i][1]
        else:
            collapsed_MD_list.append(reduced_MD_list[i])

    # The standard calls for a number to start and to end, zero if necessary.
    if isinstance(collapsed_MD_list[0], str):
        collapsed_MD_list.insert(0, 0)
    if isinstance(collapsed_MD_list[-1], str):
        collapsed_MD_list.append(0)
    
    MD_string = ''.join(map(str, collapsed_MD_list))
    return ''.join(map(str, collapsed_MD_list))

def line_groups(sam_file_name, key):
    ''' Yields (key value, list of consecutive lines from sam_file_name
        that are all transormed to key value).
    '''
    sam_file = open_to_reads(sam_file_name)
    groups = utilities.group_by(sam_file, key)
    return groups

def sort(input_file_name, output_file_name):
    sq_lines = get_sq_lines(input_file_name)
    input_file = open_to_reads(input_file_name)
    with open(output_file_name, 'w') as output_file:
        for sq_line in sq_lines:
            output_file.write(sq_line)
        external_sort.external_sort(input_file, output_file)
    
def sort_bam(input_file_name, output_file_name, by_name=False, num_threads=1):
    root, ext = os.path.splitext(output_file_name)
    samtools_command = ['samtools', 'sort']
    if by_name:
        samtools_command.append('-n')
    samtools_command.extend(['-@', str(num_threads),
                             '-T', root,
                             '-o', output_file_name,
                             input_file_name,
                            ])
    samtools_process = subprocess.Popen(samtools_command)
    samtools_process.communicate()

def merge_sorted_bam_files(input_file_names, merged_file_name):
    if len(input_file_names) == 1:
        shutil.copy(input_file_names[0], merged_file_name)
    else:
        merge_command = ['samtools', 'merge', '-f', merged_file_name] + input_file_names
        subprocess.check_call(merge_command)
    index_bam(merged_file_name)

def make_sorted_bam(sam_file_name, bam_file_name):
    bam_command = ['samtools', 'view', '-ubh', sam_file_name]
    sort_command = ['samtools', 'sort', '-T', bam_file_name, '-o', bam_file_name, '-']
    bam_process = subprocess.Popen(bam_command,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                  )
    sort_process = subprocess.Popen(sort_command,
                                    stdin=bam_process.stdout,
                                    stderr=subprocess.PIPE,
                                   )
    sort_process.communicate()
    _, err_output = bam_process.communicate()
    if bam_process.returncode:
        raise RuntimeError(err_output)

def bam_to_sam(bam_file_name, sam_file_name):
    view_command = ['samtools', 'view', '-o', sam_file_name, bam_file_name]
    subprocess.check_call(view_command)

def index_bam(bam_file_name):
    samtools_command = ['samtools', 'index', bam_file_name]
    subprocess.check_call(samtools_command)

def make_sorted_indexed_bam(sam_file_name, bam_file_name):
    make_sorted_bam(sam_file_name, bam_file_name)
    index_bam(bam_file_name)

def get_length_counts(bam_file_name, only_primary=True):
    bam_file = pysam.Samfile(bam_file_name)
    if only_primary:
        qlen_counts = Counter(ar.qlen for ar in bam_file if not ar.is_unmapped and not ar.is_secondary)
    else:
        qlen_counts = Counter(ar.qlen for ar in bam_file)
    
    return qlen_counts

def get_mapq_counts(bam_file_name):
    bam_file = pysam.Samfile(bam_file_name)
    mapq_counts = Counter(ar.mapq for ar in bam_file)
    return mapq_counts

def sam_to_fastq(sam_file_name):
    sam_file = pysam.Samfile(sam_file_name)
    for mapping in sam_file:
        if mapping.is_unmapped:
            seq = mapping.seq
            qual = mapping.qual
        elif mapping.is_reverse:
            seq = utilities.reverse_complement(mapping.seq)
            qual = mapping.qual[::-1]
        else:
            seq = mapping.seq
            qual = mapping.qual

        read = fastq.Read(mapping.qname, seq, qual)
        yield read

bam_to_fastq = sam_to_fastq

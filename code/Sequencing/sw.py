import numpy as np
from Sequencing import utilities, fastq, adapters, annotation
from sw_cython import *

def print_local_alignment(query, target, alignment_path):
    def index_to_char(seq, index):
        if index == GAP:
            return '-'
        else:
            if index >= len(seq):
                raise ValueError(index, len(seq), seq)
            return seq[index]
    
    first_query_index, first_target_index = alignment_path[0]
    last_query_index, last_target_index = alignment_path[-1]

    left_query = query[:first_query_index]
    left_target = target[:first_target_index]
    right_query = query[last_query_index + 1:]
    right_target = target[last_target_index + 1:]

    left_length = max(first_query_index, first_target_index)
    left_query = '{:>{length}s} '.format(left_query, length=left_length)
    left_target = '{:>{length}s} '.format(left_target, length=left_length)
    
    query_string = [left_query]
    target_string = [left_target]

    for q, t in alignment_path:
        query_string.append(index_to_char(query, q))
        target_string.append(index_to_char(target, t))

    query_string.extend([' ', right_query])
    target_string.extend([' ', right_target])

    print 'query ', ''.join(query_string)
    print 'target', ''.join(target_string)

def print_barcode_alignment(query, target, alignment_path):
    def index_to_char(seq, index):
        if index == -1:
            return '-'
        else:
            if index >= len(seq):
                raise ValueError(index, len(seq), seq)
            return seq[index]
    
    for _, last_query_index in alignment_path[::-1]:
        if last_query_index != GAP:
            break
    right_query = query[last_query_index + 1:]

    for last_target_index, _ in alignment_path[::-1]:
        if last_target_index != GAP:
            break
    right_target = target[last_target_index + 1:]

    query_string = []
    target_string = []

    for q, t in alignment_path:
        query_string.append(index_to_char(query, q))
        target_string.append(index_to_char(target, t))

    query_string.append(right_query)
    target_string.append(right_target)

    print 'query ', ''.join(query_string)
    print 'target', ''.join(target_string)

trimmed_annotation_fields = [
    ('original_name', 's'),
    ('insert_length', '04d'),
    ('adapter_seq', 's'),
    ('adapter_qual', 's'),
]

TrimmedAnnotation = annotation.Annotation_factory(trimmed_annotation_fields)

def trim_pairs(read_pairs, index_sequence, primer_type):
    before_R1, before_R2 = adapters.build_before_adapters(index_sequence, primer_type)

    for R1, R2, in read_pairs:
        status, insert_length, alignment = infer_insert_length(R1, R2, before_R1, before_R2)
        if status == 'illegal' or status == 'bad':
            trim_at = len(R1.seq)
            insert_length = -2
        else:
            trim_at = max(0, insert_length)

        payload_slice = slice(None, trim_at)
        adapter_slice = slice(trim_at, None)

        R1_annotation = TrimmedAnnotation(original_name=R1.name,
                                          insert_length=insert_length,
                                          adapter_seq=R1.seq[adapter_slice],
                                          adapter_qual=R1.qual[adapter_slice],
                                         )
        R1_trimmed_record = fastq.make_record(R1_annotation.identifier,
                                              R1.seq[payload_slice],
                                              R1.qual[payload_slice],
                                             )

        R2_annotation = TrimmedAnnotation(original_name=R2.name,
                                          insert_length=insert_length,
                                          adapter_seq=R2.seq[adapter_slice],
                                          adapter_qual=R2.qual[adapter_slice],
                                         )
        R2_trimmed_record = fastq.make_record(R2_annotation.identifier,
                                              R2.seq[payload_slice],
                                              R2.qual[payload_slice],
                                             )

        yield R1_trimmed_record, R2_trimmed_record

def infer_insert_length(R1, R2, before_R1, before_R2):
    ''' Infer the length of the insert represented by R1 and R2 by performing
        a semi-local alignment of R1 and the reverse complement of R2 with
        the expected adapter sequences prepended to each read.
    '''
    extended_R1 = before_R1 + R1.seq
    extended_R2 = utilities.reverse_complement(before_R2 + R2.seq)
    alignment = overlap_alignment(extended_R1, extended_R2, 2, -1, -5)
    
    R1_start = len(before_R1)
    R2_start = len(R2.seq) - 1
    R1_start_in_R2 = alignment['query_mappings'][len(before_R1)]
    R2_start_in_R1 = alignment['target_mappings'][len(R2.seq) - 1]
    
    # Since R1 is the query and R2 is the target, bases in R1 that aren't in
    # R2 are called insertions, and bases in R2 that aren't in R1 are called
    # deletions.
    # An indel in the insert is non-physical.
    if R2_start_in_R1 != SOFT_CLIPPED:
        illegal_insertion = any(R1_start <= i <= R2_start_in_R1 for i in alignment['insertions'])
    else:
        illegal_insertion = any(R1_start <= i for i in alignment['insertions'])

    if R1_start_in_R2 != SOFT_CLIPPED:
        illegal_deletion = any(R1_start_in_R2 <= d <= R2_start for d in alignment['deletions'])
    else:
        illegal_deletion = any(d <= R2_start for d in alignment['deletions'])
    
    if illegal_insertion or illegal_deletion:
        return 'illegal', 500, -1

    if R1_start_in_R2 != SOFT_CLIPPED and R2_start_in_R1 != SOFT_CLIPPED:
        length_from_R1 = R2_start_in_R1 - R1_start + 1
        length_from_R2 = R2_start - R1_start_in_R2 + 1
    else:
        # overlap_alignment forces the alignment to start with either the
        # beginning of R1 or R2 and end with either the end of R1 or R2. 
        # Making it to this else brach means that either the first base of R1 or
        # the first base of R2 or both wasn't aligned. This either means that
        # the insert is longer than the read length or a pathological alignment
        # has been produced in which only adapter bases are involved in the 
        # alignment. Flag the second case as illegal.

        first_R1_index, first_R2_index = alignment['path'][0]
        length_from_R1 = (first_R1_index - R1_start + 1) + (len(R2.seq) - 1)

        last_R1_index, last_R2_index = alignment['path'][-1]
        length_from_R2 = (R2_start - last_R2_index + 1) + (len(R1.seq) - 1)
        
        if first_R1_index == 0 or last_R2_index == 0:
            return 'illegal', 500, -1 

    if length_from_R1 < -1 or length_from_R2 < -1:
        # Negative insert lengths are non-physical. Even though I don't
        # understand it, -1 is relatively common so is tolerated.
        return 'illegal', 500, -1

    if length_from_R1 != length_from_R2:
        print length_from_R1, length_from_R2
        print_diagnostic(R1, R2, before_R1, before_R2, alignment)
        # This shouldn't be possible without an illegal indel.
        raise ValueError
    
    insert_length = length_from_R1

    if 2 * len(alignment['path']) - alignment['score'] > 20:
        status = 'bad'
    else:
        status = 'good'

    return status, insert_length, alignment

def print_diagnostic(R1, R2, before_R1, before_R2, alignment):
    extended_R1 = before_R1.lower() + R1.seq
    extended_R2 = utilities.reverse_complement(before_R2.lower() + R2.seq)
    print R1.name
    print R1.qual
    print R2.qual
    print alignment['score'], len(alignment['path']) * 2, alignment['score'] - len(alignment['path']) * 2
    print alignment['path']
    print_local_alignment(extended_R1, extended_R2, alignment['path'])
    print alignment['insertions']
    print alignment['deletions']
    print sorted(alignment['mismatches'])
    for q, t in sorted(alignment['mismatches']):
        print '\t', extended_R1[q], extended_R2[t]
    
if __name__ == '__main__':
    R1_fn = '/home/jah/projects/mutations/experiments/shiroguchi/E_coli_transcriptome_1/data/smaller_R1.fastq'
    R2_fn = '/home/jah/projects/mutations/experiments/shiroguchi/E_coli_transcriptome_1/data/smaller_R2.fastq'
    
    R1_trimmed_fn = '/home/jah/projects/mutations/experiments/shiroguchi/E_coli_transcriptome_1/results/smaller_R1_trimmed.fastq'
    R2_trimmed_fn = '/home/jah/projects/mutations/experiments/shiroguchi/E_coli_transcriptome_1/results/smaller_R2_trimmed.fastq'
    primer_type = 'PE'
    index_sequence = ''
    
    read_pairs = fastq.read_pairs(R1_fn, R2_fn)
    
    with open(R1_trimmed_fn, 'w') as R1_trimmed_fh, open(R2_trimmed_fn, 'w') as R2_trimmed_fh:
        trimmed_pairs = trim_pairs(read_pairs, index_sequence, primer_type)
        for R1_trimmed_record, R2_trimmed_record in trimmed_pairs:
            R1_trimmed_fh.write(R1_trimmed_record)
            R2_trimmed_fh.write(R2_trimmed_record)

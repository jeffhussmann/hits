import numpy as np
from Sequencing import utilities, fastq, adapters, annotation
from sw_cython import *

empty_alignment = {'score': -1e6,
                   'path': [],
                   'query_mappings': [],
                   'target_mappings': [],
                   'insertions': set(),
                   'deletions': set(),
                   'mismatches': set(),
                  }

def first_query_index(alignment_path):
    for q, t in alignment_path:
        if q != GAP:
            return q

    return -1

def first_target_index(alignment_path):
    for q, t in alignment_path:
        if t != GAP:
            return t

    return -1

def last_query_index(alignment_path):
    for q, t in alignment_path[::-1]:
        if q != GAP:
            return q

    return -1

def last_target_index(alignment_path):
    for q, t in alignment_path[::-1]:
        if t != GAP:
            return t

    return -1

def last_n_query_pairs(alignment_path, n):
    ''' Returns (q, t) pairs for the last n elements of alignment_path for which
    q is not a gap.
    '''
    pairs = []

    for q, t in alignment_path[::-1]:
        if q != GAP:
            pairs.append((q, t))
        if len(pairs) == n:
            break
    
    pairs = pairs[::-1]

    return pairs

def make_char_pairs(index_pairs, query, target):
    char_pairs = []
    for q, t in index_pairs:
        if q == GAP:
            q_char = '-'
        else:
            q_char = query[q]
        if t == GAP:
            t_char = '-'
        else:
            t_char = target[t]
        char_pairs.append((q_char, t_char))
    return char_pairs

def print_local_alignment(query, target, alignment_path):
    def index_to_char(seq, index):
        if index == GAP:
            return '-'
        else:
            if index >= len(seq):
                raise ValueError(index, len(seq), seq)
            return seq[index]
    
    first_q = first_query_index(alignment_path)
    first_t = first_target_index(alignment_path)
    last_q = last_query_index(alignment_path)
    last_t = last_target_index(alignment_path)

    left_query = query[:first_q]
    left_target = target[:first_t]
    right_query = query[last_q + 1:]
    right_target = target[last_t + 1:]

    left_length = max(first_q, first_t)
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

trimmed_annotation_fields = [
    ('original_name', 's'),
    ('insert_length', '04d'),
    ('adapter_seq', 's'),
    ('adapter_qual', 's'),
]

TrimmedAnnotation = annotation.Annotation_factory(trimmed_annotation_fields)
NO_DETECTED_OVERLAP = -2

def generate_alignments(query,
                        target,
                        alignment_type,
                        match_bonus=2,
                        mismatch_penalty=-2,
                        indel_penalty=-5,
                        max_alignments=-1,
                        min_score=60,
                       ):
    if alignment_type == 'local':
        force_query_start = False
        force_target_start = False
        force_either_start = False
        force_edge_end = False
    elif alignment_type == 'barcode':
        force_query_start = True
        force_target_start = True
        force_either_start = False
        force_edge_end = True
    elif alignment_type == 'overlap':
        force_query_start = False
        force_target_start = False
        force_either_start = True
        force_edge_end = True
    elif alignment_type == 'unpaired_adapter':
        force_query_start = True
        force_target_start = False
        force_either_start = False
        force_edge_end = True

    matrices = generate_matrices(query,
                                 target,
                                 match_bonus,
                                 mismatch_penalty,
                                 indel_penalty,
                                 force_query_start,
                                 force_target_start,
                                 force_either_start,
                                )
    cells_seen = set()
    if force_edge_end:
        possible_ends = propose_edge_ends(matrices['score'], cells_seen, min_score)
    else:
        possible_ends = propose_all_ends(matrices['score'], cells_seen, min_score)

    alignments = []
    for possible_end in possible_ends:
        #print possible_end, matrices['score'][possible_end],
        alignment = backtrack(query,
                              target,
                              matrices,
                              cells_seen,
                              possible_end,
                              force_query_start,
                              force_target_start,
                              force_either_start,
                             )
        if alignment != None:
            alignments.append(alignment)
            #print alignment['path']
            #sw.print_local_alignment(query, target, alignment['path'])
            if len(alignments) == max_alignments:
                break
        else:
            #print 'collision'
            pass

    return alignments

def propose_edge_ends(score_matrix, cells_seen, min_score):
    num_rows, num_cols = score_matrix.shape
    right_edge = [(i, num_cols - 1) for i in range(num_rows)]
    bottom_edge = [(num_rows - 1, i) for i in range(num_cols)]
    edge_cells = right_edge + bottom_edge
    sorted_edge_cells = sorted(edge_cells,
                               key=lambda cell: score_matrix[cell],
                               reverse=True,
                              )
    for cell in sorted_edge_cells:
        if score_matrix[cell] <= min_score:
            break

        if cell in cells_seen:
            continue

        yield cell

def propose_all_ends(score_matrix, cells_seen, min_score):
    sorted_indices = score_matrix.ravel().argsort()[::-1]
    for index in sorted_indices:
        cell = np.unravel_index(index, score_matrix.shape)

        if score_matrix[cell] < min_score:
            break

        if cell in cells_seen:
            continue

        yield cell

def backtrack(query,
              target,
              matrices,
              cells_seen,
              end_cell,
              force_query_start,
              force_target_start,
              force_either_start,
             ):
    query_mappings = np.ones(len(query), int) * SOFT_CLIPPED
    target_mappings = np.ones(len(target), int) * SOFT_CLIPPED

    path = []
    insertions = set()
    deletions = set()
    mismatches = set()
    
    unconstrained_start = not(force_query_start or force_target_start or force_either_start)

    row, col = end_cell

    while True:
        if (row, col) in cells_seen:
            #print (row, col), 'was already seen'
            return None
        cells_seen.add((row, col))

        next_col = col + matrices['col_direction'][row, col]
        next_row = row + matrices['row_direction'][row, col]
        if next_col == col:
            target_index = GAP
            insertions.add(row - 1)
        else:
            target_index = col - 1
        if next_row == row:
            query_index = GAP
            deletions.add(col - 1)
        else:
            query_index = row - 1
        
        if target_index != GAP:
            target_mappings[target_index] = query_index
        if query_index != GAP:
            query_mappings[query_index] = target_index
        if target_index != GAP and query_index != GAP and query[query_index] != target[target_index]:
            mismatches.add((query_index, target_index))

        path.append((query_index, target_index))

        row = next_row
        col = next_col

        if unconstrained_start:
            if matrices['score'][row, col] <= 0:
                break
        elif force_query_start and force_target_start:
            if row == 0 and col == 0:
                break
        elif force_either_start:
            if row == 0 or col == 0:
                break
        elif force_query_start:
            if row == 0:
                break
        elif force_target_start:
            if col == 0:
                break

    path = path[::-1]

    alignment = {'score': matrices['score'][end_cell],
                 'path': path,
                 'query_mappings': query_mappings,
                 'target_mappings': target_mappings,
                 'insertions': insertions,
                 'deletions': deletions,
                 'mismatches': mismatches,
                 'XM': len(mismatches),
                 'XO': len(insertions) + len(deletions),
                }

    return alignment

def infer_insert_length(R1, R2, before_R1, before_R2):
    ''' Infer the length of the insert represented by R1 and R2 by performing
        a semi-local alignment of R1 and the reverse complement of R2 with
        the expected adapter sequences prepended to each read.
    '''
    extended_R1 = before_R1 + R1.seq
    extended_R2 = utilities.reverse_complement(before_R2 + R2.seq)
    alignment,  = generate_alignments(extended_R1,
                                      extended_R2, 
                                      'overlap',
                                      2,
                                      -1,
                                      -5,
                                      1,
                                      0,
                                     )

    #print_diagnostic(R1, R2, before_R1, before_R2, alignment)
    
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
        # overlap alignment forces the alignment to start with either the
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

    insert_length = length_from_R1

    if 2 * len(alignment['path']) - alignment['score'] > .2 * len(alignment['path']):
        status = 'bad'
    else:
        status = 'good'
    
    if status == 'good' and (length_from_R1 != length_from_R2):
        print 'length from R1', length_from_R1
        print 'length from R2', length_from_R2
        print_diagnostic(R1, R2, before_R1, before_R2, alignment)
        # This shouldn't be possible without an illegal indel.
        raise ValueError

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

import numpy as np
import sys
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

    return None

def first_target_index(alignment_path):
    for q, t in alignment_path:
        if t != GAP:
            return t

    return None

def last_query_index(alignment_path):
    for q, t in alignment_path[::-1]:
        if q != GAP:
            return q

    return None

def last_target_index(alignment_path):
    for q, t in alignment_path[::-1]:
        if t != GAP:
            return t

    return None

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

def print_local_alignment(query, target, alignment_path, fh=sys.stdout):
    if alignment_path == []:
        fh.write('{0}\n\n{1}\n'.format(query, target))
        return

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

    fh.write('query\t{0}\n'.format(''.join(query_string)))
    fh.write('target\t{0}\n'.format(''.join(target_string)))

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
                        mismatch_penalty=-1,
                        indel_penalty=-5,
                        max_alignments=1,
                        min_score=None,
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
        possible_ends = propose_edge_ends(matrices['scores'], cells_seen, min_score)
    else:
        possible_ends = propose_all_ends(matrices['scores'], cells_seen, min_score)

    alignments = []
    for end_row, end_col in possible_ends:
        #print possible_end, matrices['score'][possible_end],
        alignment = backtrack_cython(query,
                                     target,
                                     matrices,
                                     cells_seen,
                                     end_row,
                                     end_col,
                                     force_query_start,
                                     force_target_start,
                                     force_either_start,
                                    )
        if alignment != None:
            alignments.append(alignment)
            #print alignment['path']
            #print_local_alignment(query, target, alignment['path'])
            if len(alignments) == max_alignments:
                break
        else:
            #print 'collision'
            pass

    return alignments

def propose_edge_ends(score_matrix,
                      cells_seen,
                      min_score=None,
                      max_alignments=1,
                     ):
    num_rows, num_cols = score_matrix.shape
    if max_alignments == 1:
        max_row = np.argmax(score_matrix[:, num_cols - 1])
        max_row_score = score_matrix[max_row, num_cols -1]
        max_col = np.argmax(score_matrix[num_rows - 1, :])
        max_col_score = score_matrix[num_rows - 1, max_col]
        if max_row_score > max_col_score:
            sorted_edge_cells = [(max_row, num_cols - 1)]
        else:
            sorted_edge_cells = [(num_rows - 1, max_col)]
    else:
        right_edge = [(i, num_cols - 1) for i in range(num_rows)]
        # Note: range(num_cols - 1) prevents including the corner twice
        bottom_edge = [(num_rows - 1, i) for i in range(num_cols - 1)]
        edge_cells = right_edge + bottom_edge
        sorted_edge_cells = sorted(edge_cells,
                                   key=lambda cell: score_matrix[cell],
                                   reverse=True,
                                  )
    for cell in sorted_edge_cells:
        if min_score != None and score_matrix[cell] < min_score:
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

def print_diagnostic(R1, R2, before_R1, before_R2, alignment, fh=sys.stdout):
    extended_R1 = before_R1.lower() + R1.seq
    extended_R2 = utilities.reverse_complement(before_R2.lower() + R2.seq)
    fh.write(R1.name + '\n')
    fh.write(R1.qual + '\n')
    fh.write(R2.qual + '\n')
    fh.write('{0}\t{1}\t{2}\n'.format(alignment['score'], len(alignment['path']) * 2, alignment['score'] - len(alignment['path']) * 2))
    fh.write(str(alignment['path']) + '\n')
    print_local_alignment(extended_R1, extended_R2, alignment['path'], fh=fh)
    fh.write(str(alignment['insertions']) + '\n')
    fh.write(str(alignment['deletions']) + '\n')
    fh.write(str(sorted(alignment['mismatches'])) + '\n')
    for q, t in sorted(alignment['mismatches']):
        fh.write('\t{0}\t{1}\n'.format(extended_R1[q], extended_R2[t]))

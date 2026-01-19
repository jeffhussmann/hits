import array
import copy
import logging
import sys
from collections import Counter, defaultdict

import numpy as np
import pysam

from . import fasta
from . import fastq
from . import genomes
from . import interval
from . import sam
from . import utilities
from .sw_cython import *

logger = logging.getLogger(__name__)

empty_alignment = {
    'score': -1e6,
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

def print_local_alignment(alignment, fh=sys.stdout):
    query = alignment['query']
    target = alignment['target']
    path = alignment['path']
    if path == []:
        fh.write(f'{query}\n\n{target}\n')
        return

    def index_to_char(seq, index):
        if index == GAP:
            return '-'
        else:
            if index >= len(seq):
                raise ValueError(index, len(seq), seq)
            return seq[index]
    
    first_q = first_query_index(path)
    first_t = first_target_index(path)
    last_q = last_query_index(path)
    last_t = last_target_index(path)

    left_query = query[:first_q]
    left_target = target[:first_t]
    right_query = query[last_q + 1:]
    right_target = target[last_t + 1:]

    left_length = max(first_q, first_t)
    left_query = '{:>{length}s} '.format(left_query, length=left_length)
    left_target = '{:>{length}s} '.format(left_target, length=left_length)
    
    query_string = [left_query]
    target_string = [left_target]

    for q, t in path:
        query_string.append(index_to_char(query, q))
        target_string.append(index_to_char(target, t))

    query_string.extend([' ', right_query])
    target_string.extend([' ', right_target])

    fh.write('query\t{0}\n'.format(''.join(query_string)))
    fh.write('target\t{0}\n'.format(''.join(target_string)))

def generate_alignments(query,
                        target,
                        alignment_type='local',
                        match_bonus=2,
                        mismatch_penalty=-1,
                        indel_penalty=-5,
                        max_alignments=1,
                        min_score=-np.inf,
                        N_matches=True,
                        insertion_penalty=None,
                        deletion_penalty=None,
                       ):

    force_query_start = False
    force_target_start = False
    force_either_start = False

    force_query_end = False
    force_target_end = False
    force_edge_end = False

    if alignment_type == 'local':
        # no forces need to be set to True
        pass

    elif alignment_type == 'edge':
        force_edge_end = True

    elif alignment_type == 'barcode':
        force_query_start = True
        force_target_start = True
        force_edge_end = True

    elif alignment_type == 'overlap':
        force_either_start = True
        force_edge_end = True
    
    elif alignment_type == 'query_end':
        force_query_end = True

    elif alignment_type == 'query_start':
        force_query_start = True

    elif alignment_type == 'unpaired_adapter':
        force_query_start = True
        force_edge_end = True

    elif alignment_type == 'IVT':
        force_query_start = True

    elif alignment_type == 'global':
        force_query_start = True
        force_target_start = True
        force_edge_end = True

    elif alignment_type == 'whole_query':
        force_query_start = True
        force_query_end = True

    elif alignment_type == 'fixed_start':
        force_query_start = True
        force_target_start = True

    elif alignment_type == 'fixed_end':
        force_query_end = True
        force_target_end = True

    query_bytes = query.encode()
    target_bytes = target.encode()

    if indel_penalty is not None:
        insertion_penalty = indel_penalty
        deletion_penalty = indel_penalty

    matrices = generate_matrices(query_bytes,
                                 target_bytes,
                                 match_bonus,
                                 mismatch_penalty,
                                 insertion_penalty,
                                 deletion_penalty,
                                 force_query_start,
                                 force_target_start,
                                 force_either_start,
                                 N_matches,
                                )
    
    cells_seen = set()

    if force_query_end and force_target_end:
        possible_ends = propose_bottom_right_corner(matrices['scores'])
    elif force_edge_end:
        possible_ends = propose_edge_ends(matrices['scores'], cells_seen, min_score, max_alignments=max_alignments)
    elif force_query_end:
        possible_ends = propose_query_edge_ends(matrices['scores'], cells_seen, min_score, max_alignments=max_alignments)
    else:
        possible_ends = propose_all_ends(matrices['scores'], cells_seen, min_score)

    alignments = []
    for end_row, end_col in possible_ends:
        alignment = backtrack_cython(query_bytes,
                                     target_bytes,
                                     matrices,
                                     cells_seen,
                                     end_row,
                                     end_col,
                                     force_query_start,
                                     force_target_start,
                                     force_either_start,
                                    )
        if alignment != None and len(alignment['path']) != 0:
            alignment['query'] = query
            alignment['target'] = target
            alignment['score_ratio'] = alignment['score'] / (len(alignment['path']) * match_bonus)
            alignments.append(alignment)
            if len(alignments) == max_alignments:
                break
        else:
            pass

    return alignments

def propose_bottom_right_corner(score_matrix):
    num_rows, num_cols = score_matrix.shape
    yield (num_rows - 1, num_cols - 1)

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

def propose_query_edge_ends(score_matrix,
                            cells_seen,
                            min_score=None,
                            max_alignments=1,
                           ):
    num_rows, num_cols = score_matrix.shape
    if max_alignments == 1:
        max_col = np.argmax(score_matrix[num_rows - 1, :])
        sorted_edge_cells = [(num_rows - 1, max_col)]
    else:
        bottom_edge = [(num_rows - 1, i) for i in range(num_cols)]
        sorted_edge_cells = sorted(bottom_edge,
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

def global_alignment(query, target, **kwargs):
    al, = generate_alignments(query, target, 'global', **kwargs)
    return al

def infer_insert_length(R1, R2, before_R1, before_R2,
                        indel_penalty=-5,
                        in_R1_prefix=None,
                        in_R2_prefix=None,
                        min_overlap_length=10,
                       ):
    ''' Infer the length of the insert represented by R1 and R2 by performing
        a semi-local alignment of R1 and the reverse complement of R2 with
        the expected adapter sequences prepended to each read.
    '''
    # Attempt a fast approach - check if expected sequencing primers occur at the same offset
    # in both reads and the implied insert in each read is close to identical.
    fast_insert_length = None

    if in_R1_prefix is not None and in_R2_prefix is not None:
        if in_R1_prefix in R1.seq and in_R2_prefix in R2.seq:
            R1_end = R1.seq.index(in_R1_prefix)
            R2_end = R2.seq.index(in_R2_prefix)
            if R1_end == R2_end:
                insert_length = R1_end
                R1_trimmed = R1.seq[:insert_length].encode()
                R2_trimmed = utilities.reverse_complement(R2.seq[:insert_length]).encode()
                
                if hamming_distance(R1_trimmed, R2_trimmed) / max(insert_length, 1) < 0.05:
                    fast_insert_length = insert_length

    if fast_insert_length is not None:
        results = {'length': fast_insert_length}
        return results

    # If not, fall back on dynamic programming.

    extended_R1 = before_R1 + R1.seq
    extended_R2 = utilities.reverse_complement(before_R2 + R2.seq)
    alignments = generate_alignments(extended_R1,
                                     extended_R2, 
                                     'overlap',
                                     2,
                                     -1,
                                     indel_penalty,
                                     1,
                                     0,
                                    )

    if len(alignments) == 0:
        results = {'failed': 'no alignments detected'}
        return results
    else:
        alignment = alignments[0]

    results = {
        'alignment': alignment,
    }

    if len(alignment['path']) < min_overlap_length:
        results['failed'] = 'overlap shorter than minimum required'
        return results

    R1_start = len(before_R1)
    R2_start = len(R2.seq) - 1
    R1_start_in_R2 = alignment['query_mappings'][R1_start]
    R2_start_in_R1 = alignment['target_mappings'][R2_start]
    R1_end_in_R2 = alignment['query_mappings'][len(extended_R1) - 1]
    R2_end_in_R1 = alignment['target_mappings'][0]
    
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
        if illegal_insertion:
            results['failed'] = 'illegal insertion'
        elif illegal_deletion:
            results['failed'] = 'illegal insertion'

        return results

    if len(alignment['path']) == 0:
        results['failed'] = 'no alignment'
        return results

    if R1_start_in_R2 != SOFT_CLIPPED and R2_start_in_R1 != SOFT_CLIPPED:
        length_from_R1 = R2_start_in_R1 - R1_start + 1
        length_from_R2 = R2_start - R1_start_in_R2 + 1
    elif R1_start_in_R2 != SOFT_CLIPPED and R1_end_in_R2 != SOFT_CLIPPED:
        # R1 entirely contained in R2
        length_from_R2 = len(R2) - R1_start_in_R2
        length_from_R1 = length_from_R2
    elif R2_start_in_R1 != SOFT_CLIPPED and R2_end_in_R1 != SOFT_CLIPPED:
        # R2 entirely contained in R1
        length_from_R1 = R2_start_in_R1 - R1_start + 1
        length_from_R2 = length_from_R1
    else:
        # overlap alignment forces the alignment to start with either the
        # beginning of R1 or R2 and end with either the end of R1 or R2. 
        # Making it to this else branch means that either the first base of R1 or
        # the first base of R2 or both wasn't aligned. This either means that
        # the insert is longer than the read length or a pathological alignment
        # has been produced in which only adapter bases are involved in the 
        # alignment. Flag the second case as illegal.

        try:
            first_R1_index, first_R2_index = alignment['path'][0]
        except IndexError:
            print(R1)
            print(R2)
            print(alignment)
            raise
        length_from_R1 = (first_R1_index - R1_start + 1) + (len(R2.seq) - 1)

        last_R1_index, last_R2_index = alignment['path'][-1]
        length_from_R2 = (R2_start - last_R2_index + 1) + (len(R1.seq) - 1)
        
        if first_R1_index == 0 or last_R2_index == 0:
            results['failed'] = 'only adapters'
            return results

    if length_from_R1 < -1 or length_from_R2 < -1:
        # Negative insert lengths are non-physical. Even though I don't
        # understand it, -1 is relatively common so is tolerated.
        results['failed'] = 'negative insert length'
        return results

    if alignment['score_ratio'] < 0.7:
        results['failed'] = 'bad score'
    else:
        if length_from_R1 != length_from_R2:
            # This shouldn't be possible without an illegal indel.
            results['failed'] = 'length mismatch'
        else:
            results['length'] = length_from_R1
    
    return results

def print_diagnostic(R1, R2, before_R1, before_R2, alignment, fh=sys.stdout):
    extended_R1 = before_R1.lower() + R1.seq
    extended_R2 = utilities.reverse_complement(before_R2.lower() + R2.seq)
    #fh.write(R1.name + '\n')
    #fh.write(R1.qual + '\n')
    #fh.write(R2.qual + '\n')
    #fh.write('{0}\t{1}\t{2}\n'.format(alignment['score'], len(alignment['path']) * .2, alignment['score'] - len(alignment['path']) * 2))
    #fh.write(str(alignment['path']) + '\n')
    print_local_alignment(alignment, fh=fh)
    #fh.write(str(alignment['insertions']) + '\n')
    #fh.write(str(alignment['deletions']) + '\n')
    #fh.write(str(sorted(alignment['mismatches'])) + '\n')
    #for q, t in sorted(alignment['mismatches']):
    #    fh.write('\t{0}\t{1}\n'.format(extended_R1[q], extended_R2[t]))

def align_read(read,
               targets,
               min_path_length,
               header,
               min_score_ratio=0.8,
               max_alignments_per_target=1,
               both_directions=True,
               read_interval=None,
               ref_intervals=None,
               **kwargs):

    if read_interval is None:
        read_interval = interval.Interval(0, len(read) - 1)

    if ref_intervals is None:
        ref_intervals = {}

    for target_name, target_seq in targets:
        if target_name not in ref_intervals:
            ref_intervals[target_name] = interval.Interval(0, len(target_seq) - 1)

    skipped_at_start = max(read_interval.start, 0)
    skipped_at_end = max(len(read) - 1 - read_interval.end, 0)

    # Register full seq and quals before possibly trimming read
    original_forward_seq = read.seq
    original_forward_qual = array.array('B', fastq.decode_sanger(read.qual))

    rc_read = read.reverse_complement()
    original_rc_seq = rc_read.seq
    original_rc_qual = array.array('B', fastq.decode_sanger(rc_read.qual))

    # Trim read to read interval
    read = read[read_interval.start:read_interval.end + 1]

    forward_seq = read.seq

    rc_read = read.reverse_complement()
    rc_seq = rc_read.seq

    alignments = []

    if both_directions:
        directions = [False, True]
    else:
        directions = [False]

    for target_name, target_seq in targets:
        ref_interval = ref_intervals[target_name]
        target_seq = target_seq[ref_interval.start:ref_interval.end + 1]

        target_alignments = []

        for is_reverse in directions:
            if is_reverse:
                seq = rc_seq

                original_seq = original_rc_seq
                original_qual = original_rc_qual

                extra_soft_clip_start = skipped_at_end
                extra_soft_clip_end = skipped_at_start

            else:
                seq = forward_seq

                original_seq = original_forward_seq
                original_qual = original_forward_qual

                extra_soft_clip_start = skipped_at_start
                extra_soft_clip_end = skipped_at_end

            for alignment in generate_alignments(seq, target_seq, max_alignments=max_alignments_per_target, **kwargs):
                path = alignment['path']

                if len(path) >= min_path_length and alignment['score_ratio'] >= min_score_ratio:
                    al = pysam.AlignedSegment(header)
                    al.seq = original_seq
                    al.query_qualities = original_qual
                    al.is_reverse = is_reverse

                    char_pairs = make_char_pairs(path, seq, target_seq)

                    cigar = sam.aligned_pairs_to_cigar(char_pairs)

                    clip_from_start = first_query_index(path) + extra_soft_clip_start
                    if clip_from_start > 0:
                        cigar = [(sam.BAM_CSOFT_CLIP, clip_from_start)] + cigar

                    clip_from_end = len(seq) - 1 - last_query_index(path) + extra_soft_clip_end
                    if clip_from_end > 0:
                        cigar = cigar + [(sam.BAM_CSOFT_CLIP, clip_from_end)]

                    al.cigar = cigar
                    
                    if al.query_length != al.infer_query_length():
                        raise ValueError(f'CIGAR implies different query length - {al.query_name}: {al.query_length}, {al.infer_query_length()}')

                    # NOTE: some evidence this might mess up if alignment ends in an indel.
                    md = sam.aligned_pairs_to_MD_string(char_pairs)
                    al.set_tag('MD', md)

                    al.set_tag('AS', alignment['score'])
                    al.reference_name = target_name
                    al.query_name = read.name
                    al.next_reference_id = -1

                    # Add back target seq that was ignored
                    al.reference_start = first_target_index(path) + ref_interval.start

                    target_alignments.append(al)

        target_alignments = sorted(target_alignments, key=lambda al: al.get_tag('AS'), reverse=True)
        alignments.extend(target_alignments[:max_alignments_per_target])

    return alignments

def align_reads(target_fasta_fn,
                reads,
                bam_fn,
                min_path_length=15,
                error_fn='/dev/null',
                alignment_type='overlap',
                yield_unaligned=False,
               ):
    ''' Aligns reads to targets in target_fasta_fn by Smith-Waterman, storing
    alignments in bam_fn and yielding unaligned reads.
    '''
    header = sam.header_from_fasta(target_fasta_fn)

    targets = fasta.to_dict(target_fasta_fn)
    targets = {r: targets[r] for r in header.references}

    alignment_sorter = sam.AlignmentSorter(bam_fn, header)

    statistics = Counter()
    
    def generator():
        with alignment_sorter:
            for read in reads:
                statistics['input'] += 1

                alignments = align_read(read, targets, alignment_type, min_path_length, header)
                
                if alignments:
                    statistics['aligned'] += 1

                    sorted_alignments = sorted(alignments, key=lambda m: m.get_tag('AS'), reverse=True)
                    grouped = utilities.group_by(sorted_alignments, key=lambda m: m.get_tag('AS'))
                    _, highest_group = next(grouped)
                    primary_already_assigned = False
                    for alignment in highest_group:
                        if len(highest_group) == 1:
                            alignment.mapping_quality = 2
                        else:
                            alignment.mapping_quality = 1

                        if not primary_already_assigned:
                            primary_already_assigned = True
                        else:
                            alignment.is_secondary = True

                        alignment_sorter.write(alignment)
                else:
                    statistics['unaligned'] += 1

                    yield read

            with open(error_fn, 'w') as error_fh:
                for key in ['input', 'aligned', 'unaligned']:
                    error_fh.write(f'{key}: {statistics[key]:,}\n')

    if yield_unaligned:
        return generator()
    else:
        for _ in generator():
            pass

def stitch_read_pair(R1, R2, before_R1='', before_R2='', indel_penalty=-5):
    if len(before_R1) >= 10 and len(before_R2) >= 10:
        in_R1_prefix = utilities.reverse_complement(before_R2)[:10]
        in_R2_prefix = utilities.reverse_complement(before_R1)[:10]
    else:
        in_R1_prefix = None
        in_R2_prefix = None

    results = infer_insert_length(R1, R2, before_R1, before_R2, indel_penalty,
                                  in_R1_prefix=in_R1_prefix,
                                  in_R2_prefix=in_R2_prefix,
                                 )

    if 'failed' in results:
        insert_length = len(R1) + len(R2)
    else:
        insert_length = results['length']

    R2_rc = R2.reverse_complement()

    R1_overlap_start = max(0, insert_length - len(R2))
    just_R1 = R1[:R1_overlap_start]
    overlap_R1 = R1[R1_overlap_start:insert_length]

    R2_overlap_start = max(0, len(R2) - insert_length)
    overlap_R2 = R2_rc[R2_overlap_start:R2_overlap_start + len(overlap_R1)]
    just_R2 = R2_rc[R2_overlap_start + len(overlap_R1):]

    overlap_seq = []
    overlap_qual = []
    for R1_s, R1_q, R2_s, R2_q in zip(overlap_R1.seq,
                                      overlap_R1.qual,
                                      overlap_R2.seq,
                                      overlap_R2.qual,
                                     ):
        if R1_q > R2_q:
            s, q = R1_s, R1_q
        else:
            s, q = R2_s, R2_q

        overlap_seq.append(s)
        overlap_qual.append(q)

    overlap = fastq.Read('', ''.join(overlap_seq), ''.join(overlap_qual))

    stitched = just_R1 + overlap + just_R2

    return stitched

class SeedAndExtender():
    def __init__(self, target_seq_bytes, max_seed_length, header, target_name, min_seed_length=1):
        if not isinstance(target_seq_bytes, bytes):
            target_seq_bytes = target_seq_bytes.encode()
            
        self.target_seq_bytes = target_seq_bytes
        self.min_seed_length = min_seed_length
        self.max_seed_length = max_seed_length
        self.header = header
        self.target_name = target_name
        
        seed_locations = defaultdict(set)
        
        for k in range(min_seed_length, max_seed_length + 1):
            for start in range(len(target_seq_bytes) - k + 1):
                k_mer = target_seq_bytes[start:start + k]
                seed_locations[k_mer].add(start)
                
        self.seed_locations = seed_locations

    def seed_and_extend(self, query_seq_bytes, true_query_start, true_query_end, query_name):
        ''' true_query_end should point one past the end of the seed '''
        if not isinstance(query_seq_bytes, bytes):
            query_seq_bytes = query_seq_bytes.encode()

        # If the requested seed is longer than max_seed_length, truncate it, then check if
        # the extended truncated seed at least covers the original seed.
        seed_length = true_query_end - true_query_start
        if seed_length > self.max_seed_length:
            initial_true_query_end = true_query_end
            true_query_end = true_query_start + self.max_seed_length
        else:
            initial_true_query_end = true_query_end

        als = []
        for is_reverse in [True, False]:
            if is_reverse:
                seed_query_start = len(query_seq_bytes) - true_query_end
                seed_query_end = len(query_seq_bytes) - true_query_start
                possibly_RCd_query_seq = utilities.reverse_complement(query_seq_bytes)

                def covers_original_seed(query_start, query_end):
                    return query_start <= (len(query_seq_bytes) - initial_true_query_end)

            else:
                seed_query_start = true_query_start
                seed_query_end = true_query_end
                possibly_RCd_query_seq = query_seq_bytes

                def covers_original_seed(query_start, query_end):
                    return query_end >= initial_true_query_end

            seed = possibly_RCd_query_seq[seed_query_start:seed_query_end]

            for target_start in self.seed_locations[seed]:
                query_start, query_end, target_start = extend_perfect_seed(possibly_RCd_query_seq,
                                                                           self.target_seq_bytes,
                                                                           seed_query_start,
                                                                           seed_query_end,
                                                                           target_start,
                                                                           target_start + len(seed),
                                                                          )
                if not covers_original_seed(query_start, query_end):
                    continue

                cigar = [(sam.BAM_CMATCH, query_end - query_start)]
                if query_start > 0:
                    cigar = [(sam.BAM_CSOFT_CLIP, query_start)] + cigar
                if query_end < len(query_seq_bytes):
                    cigar = cigar + [(sam.BAM_CSOFT_CLIP, len(query_seq_bytes) - query_end)]

                al = pysam.AlignedSegment(self.header)
                al.cigartuples = cigar
                al.is_reverse = is_reverse
                al.is_unmapped = False
                al.reference_name = self.target_name
                al.query_sequence = possibly_RCd_query_seq
                al.query_qualities = [41] * len(query_seq_bytes)
                al.query_name = query_name

                al.reference_start = target_start

                als.append(al)

        return als

def extend_alignment(initial_al, target_seq_bytes):
    query_seq_bytes = initial_al.query_sequence.encode()
    
    if not isinstance(target_seq_bytes, bytes):
        target_seq_bytes = target_seq_bytes.encode()

    new_query_start, new_query_end, new_target_start = extend_perfect_seed(query_seq_bytes,
                                                                           target_seq_bytes,
                                                                           initial_al.query_alignment_start,
                                                                           initial_al.query_alignment_end,
                                                                           initial_al.reference_start,
                                                                           initial_al.reference_end,
                                                                          )

    added_to_start = initial_al.query_alignment_start - new_query_start
    added_to_end = new_query_end - initial_al.query_alignment_end
    cigar = initial_al.cigar
    
    if added_to_start > 0:
        # Remove from starting soft clip...
        kind, length = cigar[0]
        if kind != sam.BAM_CSOFT_CLIP:
            raise ValueError(f'expected soft-clip, got {kind}')
        
        cigar[0] = (kind, length - added_to_start)

        # ... and add to subsequent match.
        kind, length = cigar[1]
        if kind != sam.BAM_CMATCH:
            raise ValueError(f'expected match, got {kind}')
            
        cigar[1] = (kind, length + added_to_start)

        if cigar[0][1] == 0:
            cigar = cigar[1:]
        
    if added_to_end > 0:
        # Remove from ending soft clip...
        kind, length = cigar[-1]
        if kind != sam.BAM_CSOFT_CLIP:
            raise ValueError(f'expected soft-clip, got {kind}')
        
        cigar[-1] = (kind, length - added_to_end)
        
        # ... and add to subsequent match.
        kind, length = cigar[-2]
        if kind != sam.BAM_CMATCH:
            raise ValueError(f'expected match, got {kind}')
            
        cigar[-2] = (kind, length + added_to_end)

        if cigar[-1][1] == 0:
            cigar = cigar[:-1]

    if added_to_start > 0 or added_to_end > 0:
        new_al = copy.deepcopy(initial_al)
        new_al.cigar = cigar
        new_al.reference_start = new_target_start
    else:
        new_al = initial_al
    
    return new_al

def extend_alignment_with_one_nt_deletion(initial_al,
                                          target_seq_bytes,
                                          extend_before=True,
                                          extend_after=True,
                                          minimum_gain=1,
                                         ):
    query_seq_bytes = initial_al.query_sequence.encode()
    
    if not isinstance(target_seq_bytes, bytes):
        target_seq_bytes = target_seq_bytes.encode()

    gained_before, gained_after = extend_perfect_seed_with_one_nt_deletion(query_seq_bytes, target_seq_bytes, initial_al.query_alignment_start, initial_al.query_alignment_end, initial_al.reference_start, initial_al.reference_end)
    cigar = initial_al.cigar

    extended = False
    
    if gained_before >= minimum_gain and extend_before:
        # Remove from starting soft clip...
        kind, length = cigar[0]
        if kind != sam.BAM_CSOFT_CLIP:
            raise ValueError(f'expected soft-clip, got {kind}')
        
        soft_clip_length = length - gained_before
        if soft_clip_length > 0:
            soft_clip = [(kind, soft_clip_length)]
        else:
            soft_clip = []
        
        # ... and add match and deletion.
        new_match = (sam.BAM_CMATCH, gained_before)
        deletion = (sam.BAM_CDEL, 1)

        cigar = soft_clip + [new_match, deletion] + cigar[1:]

        extended = True
            
    if gained_after >= minimum_gain and extend_after:
        # Remove from ending soft clip...
        kind, length = cigar[-1]
        if kind != sam.BAM_CSOFT_CLIP:
            raise ValueError(f'expected soft-clip, got {kind}')
        
        soft_clip_length = length - gained_after
        if soft_clip_length > 0:
            soft_clip = [(kind, soft_clip_length)]
        else:
            soft_clip = []
        
        # ... and add match and deletion
        new_match = (sam.BAM_CMATCH, gained_after)
        deletion = (sam.BAM_CDEL, 1)

        cigar = cigar[:-1] + [deletion, new_match] + soft_clip

        extended = True
        
    if extended:
        new_al = copy.deepcopy(initial_al)
        new_al.cigar = cigar
        if (gained_before > 0 and extend_before):
            # reference_start changes by 1 for deletion plus gained_before
            new_al.reference_start = initial_al.reference_start - 1 - gained_before

    else:
        new_al = initial_al
    
    return new_al

def extend_repeatedly(initial_al,
                      target_seq_bytes,
                      **extend_kwargs,
                     ):
    ''' Extend initial_al with additional sequence from target_seq_bytes
    separated by 1 nt deletions until doing so doesn't change the alignment.
    '''
    
    extended = initial_al # Confusing to call this 'extended' but helps make the while loop simpler.
    previous = None
    while extended is not previous:
        previous = extended
        extended = extend_alignment_with_one_nt_deletion(previous, target_seq_bytes, **extend_kwargs)

    return extended

def possibly_convert_primers_to_dictionary(primers):
    if not isinstance(primers, dict):
        if len(primers) != 2:
            raise ValueError

        primers = dict(zip(['left', 'right'], primers))

    return primers

def align_primers_to_sequence(primers, sequence_name, sequence):
    ''' Can this be replaced by align_primers_to_genome? ''' 

    primers = possibly_convert_primers_to_dictionary(primers)

    header = pysam.AlignmentHeader.from_references([sequence_name], [len(sequence)])

    mapper = SeedAndExtender(sequence.upper(), 10, header, sequence_name, min_seed_length=10)

    primer_alignments = {}

    for primer_name, primer in primers.items():
        primer = primer.upper()

        alignments = mapper.seed_and_extend(primer, len(primer) - 10, len(primer), primer_name)

        if not alignments:
            raise ValueError(f'{primer_name}: {primer} could not be located in {sequence_name}')

        primer_alignments[primer_name] = [max(alignments, key=lambda al: al.query_alignment_length)]
        
    return primer_alignments

def identify_concordant_primer_alignment_pair(primer_alignments, max_length=10000):
    ''' Find pairs of alignments that are on the same reference and point towards each other. '''

    def same_reference_name(first_al, second_al):
        return first_al.reference_name == second_al.reference_name

    def in_correct_orientation(first_al, second_al):
        ''' Check if first_al and second_al point towards each other on the same chromosome. '''
        
        if not same_reference_name(first_al, second_al):
            return None
        
        first_interval = interval.get_covered_on_ref(first_al)
        second_interval = interval.get_covered_on_ref(second_al)
        if interval.are_overlapping(first_interval, second_interval):
            return None
        
        left_al, right_al = sorted([first_al, second_al], key=lambda al: al.reference_start)

        if sam.get_strand(left_al) != '+' or sam.get_strand(right_al) != '-':
            return None
        else:
            return left_al, right_al

    def reference_extent(first_al, second_al):
        if not same_reference_name(first_al, second_al):
            return np.inf
        else:
            start = min(first_al.reference_start, second_al.reference_start)
            end = max(first_al.reference_end, second_al.reference_end)

            return end - start

    correct_orientation_pairs = []
    not_correct_orientation_pairs = []

    first_name, second_name = sorted(primer_alignments)

    for first_al in primer_alignments[first_name]:
        for second_al in primer_alignments[second_name]:
            if (oriented_pair := in_correct_orientation(first_al, second_al)) is not None:
                correct_orientation_pairs.append(oriented_pair)

            elif reference_extent(first_al, second_al) < max_length:
                not_correct_orientation_pairs.append((first_al, second_al))
                
    if len(correct_orientation_pairs) == 0:
        if len(not_correct_orientation_pairs) > 0:
            for first_al, second_al in not_correct_orientation_pairs:
                logger.warning(f'Found nearby primer alignments that don\'t point towards each other: {first_al.reference_name} {sam.get_strand(first_al)} {first_al.reference_start:,}-{first_al.reference_end:,}, {sam.get_strand(second_al)} {second_al.reference_start:,}-{second_al.reference_end:,}')

        raise ValueError(f'Could not identify primer binding sites that point towards each other.')

    # Rank pairs in the correct orientation by (shortest) amplicon length.

    correct_orientation_pairs = sorted(correct_orientation_pairs, key=lambda pair: reference_extent(*pair))

    if len(correct_orientation_pairs) > 1:
        logger.warning(f'{len(correct_orientation_pairs)} concordant primer alignments found.')
        for left_al, right_al in correct_orientation_pairs[:10]:
            logger.warning(f'Found {reference_extent(left_al, right_al):,} nt amplicon on {left_al.reference_name} from {left_al.reference_start:,} to {right_al.reference_end - 1:,} with cigars {left_al.cigarstring}, {right_al.cigarstring}')

        if len(correct_orientation_pairs) > 10:
            logger.warning(f'... and {len(correct_orientation_pairs) - 10} more')

    short_amplicons = [pair for pair in correct_orientation_pairs if reference_extent(*pair) <= 1000]
    if len(short_amplicons) > 1:
        logger.warning(f'{len(short_amplicons)} potential amplicons <= 1kb found. There is a risk the shortest amplicon may not be the intended target.')
        
    left_al, right_al = correct_orientation_pairs[0]

    if reference_extent(left_al, right_al) > max_length:
        logger.warning(f'Found {reference_extent(left_al, right_al):,} nt amplicon on {left_al.reference_name} from {left_al.reference_start:,} to {right_al.reference_end - 1:,}')
        raise ValueError(f'Could not identify an amplicon shorter than {max_length}.')

    return left_al, right_al

def amplify_sequence_with_primers(sequence, primers):
    primers = possibly_convert_primers_to_dictionary(primers)

    primer_alignments = align_primers_to_sequence(primers, 'sequence', sequence)

    left_al, right_al = identify_concordant_primer_alignment_pair(primer_alignments)

    after_left = left_al.reference_end
        
    right_start = right_al.reference_start

    left_primer = left_al.query_sequence
    # right_al is guaranteed to be rc'ed
    right_primer_rc = right_al.query_sequence
        
    return left_primer + sequence[after_left:right_start] + right_primer_rc

def find_all_matches(target_seq, query, case_sensitive=True):
    if not case_sensitive:
        target_seq = target_seq.upper()
        query = query.upper()

    all_matches = []

    for strand in ['+', '-']:
        if strand == '-':
            possibly_reversed = utilities.reverse_complement(query)
        else:
            possibly_reversed = query

        strand_matches = []

        while True:
            if len(strand_matches) == 0:
                start = 0
            else:
                _, previous_match = strand_matches[-1]
                start = previous_match + 1

            try:
                next_match = target_seq.index(possibly_reversed, start)
                strand_matches.append((strand, next_match))
            except ValueError:
                break

        all_matches.extend(strand_matches)

    return all_matches

def find_all_matches_in_genome(seq, genome, verbose=False):
    matches = []

    for ref_name, ref_seq in genome.items():
        if verbose:
            print(f'Searching {ref_name} ({len(ref_seq):,})')

        ref_matches = find_all_matches(ref_seq, seq)

        matches.extend([(ref_name, strand, position) for strand, position in ref_matches])

    return matches

def align_primers_to_genome(primers, genome, suffix_length, verbose=False):
    ''' Note: genome should be all upper '''

    header = genomes.get_header_from_dictionary(genome)

    def find_alignments_of_query_suffix_to_large_target(target_seq, target_name, query_seq, query_name, suffix_length):
        if not isinstance(target_seq, bytes) or not isinstance(query_seq, bytes):
            raise TypeError

        als = []

        query_suffix = query_seq[-suffix_length:]

        soft_clipped = len(query_seq) - len(query_suffix)

        matches = find_all_matches(target_seq, query_suffix)

        if verbose:
            print(f'{len(matches)} matches')

        for strand, position in matches:
            reverse = (strand == '-')

            if reverse:
                possibly_RCed_query_suffix = utilities.reverse_complement(query_suffix)
                possibly_RCed_query_seq = utilities.reverse_complement(query_seq)
            else:
                possibly_RCed_query_suffix = query_suffix
                possibly_RCed_query_seq = query_seq

            al = pysam.AlignedSegment(header)

            match_block = (sam.BAM_CMATCH, len(possibly_RCed_query_suffix))

            if soft_clipped > 0:
                soft_clipped_block = (sam.BAM_CSOFT_CLIP, soft_clipped)
                cigar_blocks = [soft_clipped_block, match_block]
            else:
                cigar_blocks = [match_block]

            if reverse:
                cigar_blocks = cigar_blocks[::-1]

            al.cigartuples = cigar_blocks

            al.is_reverse = reverse
            al.is_unmapped = False
            al.reference_name = target_name
            al.query_sequence = possibly_RCed_query_seq
            al.query_qualities = [41] * len(possibly_RCed_query_seq)
            al.query_name = query_name
            al.reference_start = position

            al = extend_alignment(al, target_seq)

            als.append(al)

        return als

    primer_alignments = {primer_name: [] for primer_name in primers}
    
    for ref_name, ref_seq in genome.items():
        if verbose:
            print(f'Searching {ref_name} ({len(ref_seq):,})')
 
        ref_seq_bytes = ref_seq.upper().encode()

        for primer_name, primer_seq in primers.items():
            primer_seq_bytes = primer_seq.upper().encode()
            
            als = find_alignments_of_query_suffix_to_large_target(ref_seq_bytes, ref_name, primer_seq_bytes, primer_name, suffix_length)
            
            primer_alignments[primer_name].extend(als)
            
    return primer_alignments

def find_primer_prefix_ignoring_Ns(read_seq, primer_prefix):
    mismatches = mismatches_at_offset(primer_prefix.encode(), read_seq.encode())

    indices_of_zeros = np.where(mismatches == 0)[0]

    if len(indices_of_zeros) > 0:
        start = indices_of_zeros[0]
    else:
        start = None

    return start
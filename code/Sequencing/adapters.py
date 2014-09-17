from adapters_cython import *
import numpy as np

tru_seq_R1_rc = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGA'
tru_seq_R2_rc = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

# Previous generation of adapters
paired_end_R1_rc = tru_seq_R1_rc
paired_end_R2_rc = 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG'

P5_rc = 'TCTCGGTGGTCGCCGTATCATT'
P7_rc = 'ATCTCGTATGCCGTCTTCTGCTTG'

A_tail = 'A' * 10

def build_adapters(index_sequence='', max_length=None):
    adapter_in_R1 = tru_seq_R2_rc + index_sequence + P7_rc + A_tail
    adapter_in_R2 = tru_seq_R1_rc + P5_rc + A_tail
    truncated_slice = slice(None, max_length)
    return adapter_in_R1[truncated_slice], adapter_in_R2[truncated_slice]

def build_adapter_ranges(index_sequence):
    def make_ranges(construct, names):
        cumulative_lengths = list(np.cumsum(map(len, construct)))
        bounds = zip([0] + cumulative_lengths, cumulative_lengths)
        ranges = zip(names, bounds)
        return ranges
    
    R1_construct = [tru_seq_R2_rc, index_sequence, P7_rc, A_tail]
    R1_names = ['TruSeq R2', 'I7', 'P7', 'A tail']

    chemistry_only_cycles = 7
    I5_length = 8
    R2_construct = [tru_seq_R1_rc[:-(I5_length + chemistry_only_cycles)],
                    tru_seq_R1_rc[-(I5_length + chemistry_only_cycles):-(chemistry_only_cycles)],
                    tru_seq_R1_rc[-chemistry_only_cycles:],
                    P5_rc,
                    A_tail,
                   ]
    R2_names = ['TruSeq R1',
                'I5',
                'Chemistry only',
                'P7',
                'A tail',
               ]
    
    R1_ranges = make_ranges(R1_construct, R1_names)
    R2_ranges = make_ranges(R2_construct, R2_names)
    return R1_ranges, R2_ranges

def get_barcode_ranges(index_sequence):
    start = len(tru_seq_R2_rc)
    end = start + len(index_sequence)
    barcode_in_R1_range = range(start, end)
    barcode_in_R2_range = []
    return barcode_in_R1_range, barcode_in_R2_range

def consistent_paired_position(R1_seq,
                               R2_seq,
                               adapter_in_R1,
                               adapter_in_R2,
                               min_comparison_length,
                               max_distance,
                               allow_prefix=True,
                              ):
    R1_positions = find_adapter_positions(R1_seq, adapter_in_R1, min_comparison_length, max_distance)
    R2_positions = find_adapter_positions(R2_seq, adapter_in_R2, min_comparison_length, max_distance)

    if allow_prefix:
        R1_prefix_position = find_adapter(adapter_in_R1, max_distance, R1_seq)
        if R1_prefix_position != len(R1_seq):
            R1_positions.append(R1_prefix_position)
        
        R2_prefix_position = find_adapter(adapter_in_R2, max_distance, R2_seq)
        if R2_prefix_position != len(R2_seq):
            R2_positions.append(R2_prefix_position)

    R1_positions = set(R1_positions)
    R2_positions = set(R2_positions)
    common_positions = R1_positions & R2_positions
    if common_positions:
        return min(common_positions)
    else:
        return None

def indel_position(R1_seq,
                   R2_seq,
                   adapter_in_R1,
                   adapter_in_R2,
                   min_comparison_length,
                   max_distance,
                  ):
    R1_positions = find_adapter_positions(R1_seq, adapter_in_R1, min_comparison_length, max_distance)
    R2_positions = find_adapter_positions(R2_seq, adapter_in_R2, min_comparison_length, max_distance)

    R1_positions = set(R1_positions)
    R2_positions = set(R2_positions)

    return R1_positions, R2_positions

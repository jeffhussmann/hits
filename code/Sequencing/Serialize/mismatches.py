import numpy as np
from Sequencing.utilities import base_order, base_to_index
from Circles import variants

def to_serialize_shape(type_counts):
    num_positions, num_quals, num_from, num_to = type_counts.shape
    num_types = num_from * num_to
    reordered = np.zeros((num_positions, num_quals, num_types), int)
    coords_order = variants.coords_order[:num_types]
    for p in range(num_positions):
        for q in range(num_quals):
            reordered[p, q] = [type_counts[p, q][cs] for cs in coords_order]
    return reordered

def to_standard_shape(reordered):
    num_positions, num_quals, num_types = reordered.shape
    num_from = int(np.sqrt(num_types))
    assert num_types == num_from**2
    num_to = num_from
    type_counts = np.zeros((num_positions, num_quals, num_from, num_to), int)
    coords_order = variants.coords_order[:num_types]
    for p in range(num_positions):
        for q in range(num_quals):
            for cs, count in zip(coords_order, reordered[p, q]):
                type_counts[p, q][cs] = count
    return type_counts

def write_file(type_counts, mismatch_file_name):
    num_positions, num_quals, num_from, num_to = type_counts.shape
    num_types = num_from * num_to
    with open(mismatch_file_name, 'w') as mismatch_file:
        shape_line = '#{0},{1},{2}\n'.format(num_positions, num_quals, num_types)
        mismatch_file.write(shape_line)
        reordered = to_serialize_shape(type_counts)
        for p in range(num_positions):
            mismatch_file.write('#{0}\n'.format(p))
            np.savetxt(mismatch_file, reordered[p], fmt='%i', delimiter='\t')

def read_file(mismatch_file_name):
    with open(mismatch_file_name) as mismatch_file:
        shape = mismatch_file.readline().strip('#\n')
        num_positions, num_quals, num_types = map(int, shape.split(','))
        reordered = np.zeros((num_positions, num_quals, num_types), int)
        for p in range(num_positions):
            delimiter_line = mismatch_file.readline()
            p_lines = (mismatch_file.readline() for q in range(num_quals))
            reordered[p] = np.loadtxt(p_lines, dtype=int)
        
    type_counts = to_standard_shape(reordered)
    return type_counts

def restrict(type_counts):
    ''' Leave out N and - from types. '''
    restricted = type_counts[:, :, :4, :4]
    return restricted

def combine_data(first_type_counts, second_type_counts):
    assert first_type_counts.shape == second_type_counts.shape
    return first_type_counts + second_type_counts

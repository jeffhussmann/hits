import numpy as np
from Circles import variants

def write_file(reordered, mismatch_file_name):
    num_positions, num_quals, num_types  = reordered.shape
    with open(mismatch_file_name, 'w') as mismatch_file:
        shape_line = '#{0},{1},{2}\n'.format(num_positions, num_quals, num_types)
        mismatch_file.write(shape_line)
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
        
    return reordered

def combine_data(first_type_counts, second_type_counts):
    assert first_type_counts.shape == second_type_counts.shape
    return first_type_counts + second_type_counts

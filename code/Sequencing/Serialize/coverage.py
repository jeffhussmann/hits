from itertools import izip
import numpy as np

def read_file(coverage_file_name):
    lines = open(coverage_file_name)
    line_pairs = izip(*[lines]*2)
    counts = {}
    for name_line, counts_line in line_pairs:
        name = name_line.strip().lstrip('>')
        these_counts = np.fromstring(counts_line.strip(), dtype=int, sep=' ')
        counts[name] = these_counts
    return counts

def combine_data(first_counts, second_counts):
    assert first_counts.viewkeys() == second_counts.viewkeys()

    counts = {name: first_counts[name] + second_counts[name] for name in first_counts}

    return counts

def write_file(counts, coverage_file_name):
    with open(coverage_file_name, 'w') as coverage_file:
        for name in sorted(counts):
            coverage_file.write('>{0}\n'.format(name))
            counts_string = ' '.join(str(c) for c in counts[name])
            coverage_file.write('{0}\n'.format(counts_string))

import numpy as np
from itertools import izip

def read_file(file_name):
    genes = {}

    lines = open(file_name)
    line_pairs = izip(*[lines]*2)

    for name_line, count_line in line_pairs:
        gene_name, length = name_line.strip().split()
        length = int(length)
        strand_counts = np.fromstring(count_line, dtype=int, sep=' ' )
        genes[gene_name] = {'CDS_length': length,
                            'expression': strand_counts,
                           }
    return genes

def write_file(data, file_name):
    with open(file_name, 'w') as output_file:
        for gene_name in sorted(data):
            info_line = '{0}\t{1}\n'.format(gene_name, data[gene_name]['CDS_length'])
            output_file.write(info_line)
            strand_counts = data[gene_name]['expression']
            counts_string = ' '.join(str(c) for c in strand_counts)
            output_file.write('{0}\n'.format(counts_string))

def combine_data(first_data, second_data):
    for gene in second_data:
        if gene not in first_data:
            first_data[gene] = second_data[gene]
        else:
            raise ValueError
    return first_data

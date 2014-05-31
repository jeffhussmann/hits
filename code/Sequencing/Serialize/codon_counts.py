import numpy as np
from itertools import izip

def write_file(codon_counts, file_name):
    with open(file_name, 'w') as fh:
        for name in sorted(codon_counts):
            counts_string = '\t'.join(str(count) for count in codon_counts[name])
            line = '{0}\t{1}\n'.format(name, counts_string)
            fh.write(line)

def read_file(fn):
    genes = {}
    for line in open(fn):
        name, values = line.strip().split('\t', 1)
        values = np.fromstring(values, dtype=int, sep='\t')
        # We are interested in counts in A-sites, for which no reliable
        # information can be obtained for the start codon or first codon
        # after this.
        # Responsibility of analysis to ignore first two values and last value.
        genes[name] = values
            
    return genes

def combine_data(first_data, second_data):
    for gene in second_data:
        if gene not in first_data:
            first_data[gene] = second_data[gene]
        else:
            first_data[gene] += second_data[gene]
    return first_data

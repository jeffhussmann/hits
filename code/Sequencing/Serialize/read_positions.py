''' Format of file: 
    The first line starts with '#' followed by a comma-separated list of
    fragment lengths for which positions are reported.
    The second line starts with a '#' followed by the length of the left buffer
    and length of the right buffer.
    The rest of the file consists of set of (number of relevant lengths + 2)
    lines. Each set starts with a line consisting of a gene name and an extent
    length, then (number of relevant lengths) lines of positions for each
    relevant length, then a line of positions for all lengths.
'''
from itertools import izip
import numpy as np
import positions

def possibly_int(string):
    try:
        value = int(string)
    except ValueError:
        value = string
    return value

def read_file(file_name):
    genes = {}

    lines = open(file_name)

    relevant_lengths_line = lines.next()
    fields = relevant_lengths_line.strip().lstrip('#').split(',')
    relevant_lengths = map(possibly_int, fields)

    buffer_lengths_line = lines.next()
    buffer_lengths_string = buffer_lengths_line.strip().lstrip('#')
    left_buffer, right_buffer = map(int, buffer_lengths_string.split(','))

    # Grab sets of len(relevant_lengths) + 1 lines, grouped into 1 and 
    # len(relevant_lengths). 
    info_lines = lines
    count_lines_iter = izip(*[lines]*len(relevant_lengths))
    grouped_lines = izip(info_lines, count_lines_iter)

    for info_line, count_lines in grouped_lines:
        gene_name = info_line.strip()
        
        position_counts = [positions.PositionCounts.from_string(line, left_buffer, right_buffer)
                           for line in count_lines]

        genes[gene_name] = dict(zip(relevant_lengths, position_counts))

    return genes

def write_file(data, file_name):
    with open(file_name, 'w') as output_file:
        random_gene = data.itervalues().next()

        # Note: 'all' or 'codons' sorts after any number
        relevant_lengths = sorted(random_gene.keys())
        relevant_lengths_string = ','.join(map(str, relevant_lengths))
        relevant_lengths_line = '#{0}\n'.format(relevant_lengths_string)
        output_file.write(relevant_lengths_line)

        random_length = relevant_lengths[0]
        random_counts = random_gene[random_length]
        buffer_lengths_line = '#{0},{1}\n'.format(random_counts.left_buffer,
                                                  random_counts.right_buffer,
                                                 )
        output_file.write(buffer_lengths_line)

        for gene_name in sorted(data):
            output_file.write('{0}\n'.format(gene_name))
            for length in relevant_lengths:
                counts = data[gene_name][length]
                output_file.write(str(counts))

def combine_data(first_data, second_data):
    for gene in second_data:
        if gene not in first_data:
            first_data[gene] = second_data[gene]
        else:
            for length in first_data[gene]:
                first_data[gene][length] += second_data[gene][length] 
    return first_data

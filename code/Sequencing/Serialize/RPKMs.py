def write_file(RPKMs, RPKMs_fn):
    with open(RPKMs_fn, 'w') as RPKMs_fh:
        for gene_name in sorted(RPKMs):
            line = '{0}\t{1:0.2f}\n'.format(gene_name, RPKMs[gene_name])
            RPKMs_fh.write(line)

def read_file(RPKMs_fn):
    genes = {}
    for line in open(RPKMs_fn):
        name, value = line.strip().split()
        genes[name] = float(value)

    return genes

def combine_data(first_data, second_data):
    for gene in second_data:
        if gene not in first_data:
            first_data[gene] = second_data[gene]
        else:
            raise ValueError
    return first_data

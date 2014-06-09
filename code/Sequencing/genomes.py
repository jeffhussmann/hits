import pysam
import glob
import os
from collections import namedtuple
import Bio.SeqIO

def get_all_fasta_file_names(genome_dir):
    fasta_file_names = [fn for fn in glob.glob('{0}/*.fa*'.format(genome_dir))
                        if not fn.endswith('.fai')]
    return fasta_file_names

def make_fais(genome_dir):
    fasta_file_names = get_all_fasta_file_names(genome_dir)
    map(pysam.faidx, fasta_file_names)

fai_entry_fields = ['file_name',
                    'length',
                    'offset',
                    'bases_per_line',
                    'bytes_per_line',
                   ]
fai_entry = namedtuple('fai_entry', fai_entry_fields)

def parse_fai(fai_file_name):
    fasta_file_name, _ = os.path.splitext(fai_file_name) 
    def parse_line(line):
        fields = line.strip().split()
        seq_name = fields[0]
        values = map(int, fields[1:])
        return seq_name, fai_entry(fasta_file_name, *values)
    entries = dict(parse_line(line) for line in open(fai_file_name))
    return entries

def get_genome_index(genome_dir):
    fasta_file_names = get_all_fasta_file_names(genome_dir)
    entries = {}
    for fasta_file_name in fasta_file_names:
        fai_file_name = fasta_file_name + '.fai'
        entries.update(parse_fai(fai_file_name))
    return entries

def build_region_fetcher(genome_dir, load_references=False):
    ''' Returns a function for fetching regions from the genome in genome_dir.
        If load_references == True, loads entire reference sequences into memory
        the first time they are fetched from.
        If the returned function is given a negative start or an end that is
        longer than the seq_name's sequence, the region returned will be
        truncated.
    '''
    genome_index = get_genome_index(genome_dir)
    fasta_files = {fasta_file_name: pysam.Fastafile(fasta_file_name)
                   for fasta_file_name in get_all_fasta_file_names(genome_dir)}
    seq_name_to_file = {seq_name: fasta_files[genome_index[seq_name].file_name]
                        for seq_name in genome_index}

    if load_references:
        references = {}
        def region_fetcher(seq_name, start, end):
            if start < 0:
                start = 0
            if seq_name not in references:
                references[seq_name] = seq_name_to_file[seq_name].fetch(seq_name)
            return references[seq_name][start:end]
    else:   
        def region_fetcher(seq_name, start, end):
            if start < 0:
                start = 0
            region = seq_name_to_file[seq_name].fetch(seq_name, start, end)
            return region

    return region_fetcher

def load_entire_genome(genome_dir):
    seqs = {}
    fasta_file_names = get_all_fasta_file_names(genome_dir)
    for fasta_file_name in fasta_file_names:
        for record in Bio.SeqIO.parse(fasta_file_name, 'fasta'):
            seqs[record.id] = str(record.seq)
    return seqs

def max_RNAME_length(genome_index):
    length = max(len(name) for name in genome_index)
    return length

def max_POS_length(genome_index):
    length = max(len(str(genome_index[name].length)) for name in genome_index)
    return length

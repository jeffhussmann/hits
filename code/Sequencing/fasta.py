from itertools import izip
from collections import namedtuple

Read = namedtuple('Read', ['name', 'seq'])

def line_groups(file_name):
    if type(file_name) == str:
        lines = open(file_name)
    else:
        lines = iter(file_name)
    return izip(*[lines]*2)

def reads(file_name):
    ''' Yields the name and sequence lines from a fasta file. '''
    for name_line, seq_line in line_groups(file_name):
        name = name_line.rstrip().lstrip('>')
        seq = seq_line.strip()
        read = Read(name, seq)
        yield read

make_record = '>{0}\n{1}\n'.format

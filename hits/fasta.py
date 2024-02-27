from collections import OrderedDict
from pathlib import Path

import Bio.SeqIO
import pandas as pd
import pysam

from . import utilities

class Record(object):
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def __str__(self):
        return make_record(self.name, self.seq)
    
    def reverse_complement(self):
        return Read(self.name,
                    utilities.reverse_complement(self.seq),
                   )
    
    def __getitem__(self, sl):
        return Read(self.name, self.seq[sl])
    
    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.seq)

    def __add__(self, other):
        return Read(self.name, self.seq + other.seq)


make_record = '>{0}\n{1}\n'.format

def records(file_name):
    ''' Yields the name and sequence lines from a fasta file. '''
    file_name = Path(file_name)

    if file_name.suffix == '.ab1':
        format = 'abi'
    else:
        format = 'fasta'

    for record in Bio.SeqIO.parse(file_name, format):
        yield Record(record.name, str(record.seq).upper())

Read = Record # for backwards compatibility
reads = records

def to_dict(file_name, upper_case=False):
    return OrderedDict((r.name, r.seq.upper() if upper_case else r.seq) for r in reads(file_name))

def write_dict(seq_dict, fn):
    with open(fn, 'w') as fh:
        for name, seq in seq_dict.items():
            record = Record(name, seq)
            fh.write(str(record))

    pysam.faidx(str(fn))

def load_fai(fasta_fn):
    fasta_fn = Path(fasta_fn)
    fai_fn = fasta_fn.with_suffix(fasta_fn.suffix + '.fai')

    if not fai_fn.exists():
        pysam.faidx(str(fasta_fn))

    column_names = [
        'NAME',
        'LENGTH',
        'OFFSET',
        'LINEBASES',
        'LINEWIDTH',
    ]

    fai = pd.read_csv(fai_fn, sep='\t', index_col=0, header=None, names=column_names)

    return fai

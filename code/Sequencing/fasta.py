from itertools import izip
from collections import namedtuple
import Bio.SeqIO

Read = namedtuple('Read', ['name', 'seq'])
make_record = '>{0}\n{1}\n'.format
def Read_to_record(self):
    return make_record(*self)
Read.__str__ = Read_to_record

def reads(file_name):
    ''' Yields the name and sequence lines from a fasta file. '''
    for record in Bio.SeqIO.parse(file_name, 'fasta'):
        read = Read(record.name, str(record.seq).upper())
        yield read

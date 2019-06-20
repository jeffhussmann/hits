import struct

MAGIC_STRING = b'BUS\x00'
PAD = 0

def read_header(fh):
    magic_string = fh.read(4)
    if magic_string != MAGIC_STRING:
        raise ValueError(magic_string)
        
    version, = struct.unpack('I', fh.read(4))
    barcode_length, = struct.unpack('I', fh.read(4))
    umi_length, = struct.unpack('I', fh.read(4))
    
    header_length, = struct.unpack('I', fh.read(4))
    header = fh.read(header_length)

    return barcode_length, umi_length

def write_header(fh, barcode_length, umi_length):
    version = 1
    fh.write(MAGIC_STRING)
    fh.write(struct.pack('I', version))
    fh.write(struct.pack('I', barcode_length))
    fh.write(struct.pack('I', umi_length))

    header = b'written by hits'
    fh.write(struct.pack('I', len(header)))
    fh.write(header)

i_to_b = 'ACGT'
b_to_i = {b: i for i, b in enumerate(i_to_b)}

def twobit_to_seq(n, length):
    seq = ''.join(i_to_b[(n >> (2 * (length - 1 - i))) & 3] for i in range(length))
    return seq

def seq_to_twobit(seq):
    b = seq[0]
    n = b_to_i[b]
    
    for b in seq[1:]:
        n <<= 2
        i = b_to_i[b]
        n |= i
        
    return n

RecordStruct = struct.Struct('LLiIII')

class Record():
    def __init__(self, barcode, umi, eq_class, count, flags):
        self.barcode = barcode
        self.umi = umi
        self.eq_class = eq_class
        self.count = count
        self.flags = flags

    @classmethod
    def from_file_handle(cls, fh, barcode_length, umi_length):
        bs = fh.read(RecordStruct.size)
        if len(bs) != RecordStruct.size:
            if len(bs) > 0:
                raise ValueError('malformed file')
            else:
                raise EOFError

        barcode, umi, eq_class, count, flags, pad = RecordStruct.unpack(bs)
        return cls(twobit_to_seq(barcode, barcode_length),
                   twobit_to_seq(umi, umi_length),
                   eq_class,
                   count,
                   flags,
                  )

    def pack(self):
        bs = RecordStruct.pack(seq_to_twobit(self.barcode),
                               seq_to_twobit(self.umi),
                               self.eq_class,
                               self.count,
                               self.flags,
                               PAD,
                              )
        return bs

def records(fn):
    with open(fn, 'rb') as fh:
        barcode_length, umi_length = read_header(fh)
        try:
            while True:
                record = Record.from_file_handle(fh, barcode_length, umi_length)
                yield record
        except EOFError:
            pass

def write_file(fn, records, barcode_length, umi_length):
    with open(fn, 'wb') as fh:
        write_header(fh, barcode_length, umi_length)
        for record in records:
            fh.write(record.pack())
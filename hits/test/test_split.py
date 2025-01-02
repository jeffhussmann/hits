from collections import defaultdict

import pytest

from hits import parallel
from hits.parallel import split_file as split_file
from hits import fastq as fastq

def test_piece_of_list():
    ''' Tests whether joining all which_pieces from piece_of_list recreates
        the original list.
    '''
    list_sizes = range(1000, 1100)
    for size in list_sizes:
        original_list = list(range(size))
        num_pieces = 96
        pieces = []
        pieces_joined = []
        for which_piece in range(num_pieces):
            piece = parallel.piece_of_list(original_list,
                                           num_pieces,
                                           which_piece,
                                          )
            pieces.append(piece)
            pieces_joined.extend(piece)

        assert (set(original_list) == set(pieces_joined))

        piece_lengths = [len(p) for p in pieces]
        assert (max(piece_lengths) - min(piece_lengths) <= 1)

@pytest.mark.skip
def test_split_fastq():
    ''' Tests whether every read in a fastq file shows up exactly once in
    piece of a split.
    '''
    fn = 'data.fastq'

    whole = defaultdict(list)

    for read in fastq.reads(fn):
        whole[read.name].append(read)

    num_pieces_list = [1, 10, 100]
    from_pieces = {n: defaultdict(list) for n in num_pieces_list}

    for num_pieces in num_pieces_list:
        for which_piece in range(num_pieces):
            piece = split_file.piece(fn, num_pieces, which_piece, 'fastq')

            for read in fastq.reads(piece):
                from_pieces[num_pieces][read.name].append(read)

        assert (whole == from_pieces[num_pieces])
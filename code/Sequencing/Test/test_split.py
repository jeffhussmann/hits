import unittest
import Sequencing.Parallel

class TestSplit(unittest.TestCase):
    def test_piece_of_list(self):
        ''' Tests whether joining all which_pieces from piece_of_list recreates
            the original list.
        '''
        piece_of_list = Sequencing.Parallel.piece_of_list
        list_sizes = range(1000, 1100)
        for size in list_sizes:
            original_list = range(size)
            num_pieces = 96
            pieces = []
            pieces_joined = []
            for which_piece in range(num_pieces):
                piece = piece_of_list(original_list,
                                      num_pieces,
                                      which_piece,
                                     )
                pieces.append(piece)
                pieces_joined.extend(piece)

            self.assertEqual(set(original_list), set(pieces_joined),
                             msg='Failed set equality for list size {0}'.format(size),
                            )

            piece_lengths = map(len, pieces)
            self.assertTrue(max(piece_lengths) - min(piece_lengths) <= 1,
                            msg='Pieces too variable in size for list size {0}'.format(size),
                           )


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSplit)
    unittest.TextTestRunner(verbosity=2).run(suite)

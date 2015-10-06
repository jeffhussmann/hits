import numpy as np

def piece_of_list(full_list, num_pieces, which_piece, interleaved=False):
    if which_piece == -1:
        # Sentinel value that means 'the whole thing'.
        piece = full_list
    elif interleaved:
        piece = [x for i, x in enumerate(full_list) if i % num_pieces == which_piece]
    else:
        base_size = len(full_list) // num_pieces
        num_with_extra = len(full_list) - base_size * num_pieces
        sizes = np.ones(num_pieces, int) * base_size
        sizes[:num_with_extra] += 1
        bounds = np.append([0], sizes.cumsum())
        piece = full_list[bounds[which_piece]:bounds[which_piece + 1]]

    return piece

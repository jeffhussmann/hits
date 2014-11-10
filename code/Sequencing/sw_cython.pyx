import numpy as np
cimport numpy as np
cimport cython

DTYPEINT = np.int
ctypedef np.int_t DTYPEINT_t

SOFT_CLIPPED = -2
GAP = -1

@cython.boundscheck(False)
def generate_matrices(char* query,
                      char* target,
                      int match_bonus,
                      int mismatch_penalty,
                      int indel_penalty,
                      force_query_start,
                      force_target_start,
                      force_either_start,
                     ):
    cdef unsigned int row, col, next_col, next_row
    cdef int match_or_mismatch, diagonal, from_left, from_above, new_score, unconstrained_start
    shape = (len(query) + 1, len(target) + 1)
    cdef np.ndarray[DTYPEINT_t, ndim=2] score_matrix = np.zeros(shape, DTYPEINT)
    cdef np.ndarray[DTYPEINT_t, ndim=2] row_direction_matrix = np.zeros(shape, DTYPEINT)
    cdef np.ndarray[DTYPEINT_t, ndim=2] col_direction_matrix = np.zeros(shape, DTYPEINT)

    # If the alignment is constrained to include the start of the query,
    # indel penalties need to be applied to cells in the first row.
    if force_query_start: 
        for row in range(1, len(query) + 1):
            score_matrix[row, 0] = score_matrix[row - 1, 0] + indel_penalty
            row_direction_matrix[row, 0] = -1

    # If the alignment is constrained to include the start of the target,
    # indel penalties need to be applied to cells in the first column.
    if force_target_start:
        for col in range(1, len(target) + 1):
            score_matrix[0, col] = score_matrix[0, col - 1] + indel_penalty
            col_direction_matrix[0, col] = -1

    unconstrained_start = not (force_query_start or force_target_start or force_either_start)

    for row in range(1, len(query) + 1):
        for col in range(1, len(target) + 1):
            if query[row - 1] == target[col - 1]:
                match_or_mismatch = match_bonus
            else:
                match_or_mismatch = mismatch_penalty
            diagonal = score_matrix[row - 1, col - 1] + match_or_mismatch
            from_left = score_matrix[row, col - 1] + indel_penalty
            from_above = score_matrix[row - 1, col] + indel_penalty
            new_score = max(diagonal, from_left, from_above)
            if unconstrained_start:
                new_score = max(0, new_score)
            score_matrix[row, col] = new_score
            if new_score > max_score:
                max_score = new_score
                max_row = row
                max_col = col
            if unconstrained_start and new_score == 0:
                pass
            elif new_score == diagonal:
                col_direction_matrix[row, col] = -1
                row_direction_matrix[row, col] = -1
            elif new_score == from_left:
                col_direction_matrix[row, col] = -1
            elif new_score == from_above:
                row_direction_matrix[row, col] = -1

    matrices = {'score': score_matrix,
                'row_direction': row_direction_matrix,
                'col_direction': col_direction_matrix,
               }
    return matrices

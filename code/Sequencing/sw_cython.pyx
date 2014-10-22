import numpy as np
cimport numpy as np
cimport cython

DTYPEINT = np.int
ctypedef np.int_t DTYPEINT_t

@cython.boundscheck(False)
def local_alignment(char* query,
                    char* target,
                    int match,
                    int mismatch,
                    int indel,
                   ):
    cdef unsigned int row, col, next_col, next_row, max_row, max_col
    cdef int match_or_mismatch, diagonal, from_left, from_above, new_score, target_index, query_index, max_score
    cdef int num_indels
    shape = (len(query) + 1, len(target) + 1)
    cdef np.ndarray[DTYPEINT_t, ndim=2] score_matrix = np.zeros(shape, DTYPEINT)
    cdef np.ndarray[DTYPEINT_t, ndim=2] row_direction_matrix = np.zeros(shape, DTYPEINT)
    cdef np.ndarray[DTYPEINT_t, ndim=2] col_direction_matrix = np.zeros(shape, DTYPEINT)

    max_score = 0
    max_row = 0
    max_col = 0
    for row in range(1, len(query) + 1):
        for col in range(1, len(target) + 1):
            if query[row - 1] == target[col - 1]:
                match_or_mismatch = match
            elif query[row - 1] == 'N':
                match_or_mismatch = match
            elif target[col - 1] == 'N':
                match_or_mismatch = match
            else:
                match_or_mismatch = mismatch
            diagonal = score_matrix[row - 1, col - 1] + match_or_mismatch
            from_left = score_matrix[row, col - 1] + indel
            from_above = score_matrix[row - 1, col] + indel
            new_score = max(0, diagonal, from_left, from_above)
            score_matrix[row, col] = new_score
            if new_score > max_score:
                max_score = new_score
                max_row = row
                max_col = col
            if new_score == 0:
                pass
            elif new_score == diagonal:
                col_direction_matrix[row, col] = -1
                row_direction_matrix[row, col] = -1
            elif new_score == from_left:
                col_direction_matrix[row, col] = -1
            elif new_score == from_above:
                row_direction_matrix[row, col] = -1

    row, col = max_row, max_col
    if max_score <= 0:
        return 0, []
    alignment = []
    insertions = set()
    deletions = set()
    while True:
        next_col = col + col_direction_matrix[row, col]
        next_row = row + row_direction_matrix[row, col]
        if next_col == col:
            target_index = -1
            deletions.add(row - 1)
        else:
            target_index = col - 1
        if next_row == row:
            query_index = -1
            insertions.add(row - 1 + 0.5)
        else:
            query_index = row - 1
        alignment.append((target_index, query_index))
        row = next_row
        col = next_col
        if score_matrix[next_row, next_col] <= 0:
            break

    alignment = alignment[::-1]
    return max_score, alignment, insertions, deletions

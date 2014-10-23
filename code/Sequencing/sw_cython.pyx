import numpy as np
cimport numpy as np
cimport cython

DTYPEINT = np.int
ctypedef np.int_t DTYPEINT_t

SOFT_CLIPPED = -2
GAP = -1

@cython.boundscheck(False)
def local_alignment(char* query,
                    char* target,
                    int match_bonus,
                    int mismatch_penalty,
                    int indel_penalty,
                   ):
    cdef unsigned int row, col, next_col, next_row, max_row, max_col
    cdef int match_or_mismatch, diagonal, from_left, from_above, new_score, target_index, query_index, max_score
    cdef int num_indels
    shape = (len(query) + 1, len(target) + 1)
    cdef np.ndarray[DTYPEINT_t, ndim=2] score_matrix = np.zeros(shape, DTYPEINT)
    cdef np.ndarray[DTYPEINT_t, ndim=2] row_direction_matrix = np.zeros(shape, DTYPEINT)
    cdef np.ndarray[DTYPEINT_t, ndim=2] col_direction_matrix = np.zeros(shape, DTYPEINT)

    cdef np.ndarray[DTYPEINT_t, ndim=1] query_mappings = np.ones(len(query), DTYPEINT) * SOFT_CLIPPED
    cdef np.ndarray[DTYPEINT_t, ndim=1] target_mappings = np.ones(len(target), DTYPEINT) * SOFT_CLIPPED

    max_score = 0
    max_row = 0
    max_col = 0
    for row in range(1, len(query) + 1):
        for col in range(1, len(target) + 1):
            if query[row - 1] == target[col - 1]:
                match_or_mismatch = match_bonus
            else:
                match_or_mismatch = mismatch_penalty
            diagonal = score_matrix[row - 1, col - 1] + match_or_mismatch
            from_left = score_matrix[row, col - 1] + indel_penalty
            from_above = score_matrix[row - 1, col] + indel_penalty
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
    
    path = []
    insertions = set()
    deletions = set()
    mismatches = set()
    while True:
        next_col = col + col_direction_matrix[row, col]
        next_row = row + row_direction_matrix[row, col]
        if next_col == col:
            target_index = GAP
            insertions.add(row - 1)
        else:
            target_index = col - 1
        if next_row == row:
            query_index = GAP
            deletions.add(col - 1)
        else:
            query_index = row - 1
        
        if target_index != GAP:
            target_mappings[target_index] = query_index
        if query_index != GAP:
            query_mappings[query_index] = target_index
        if target_index != GAP and query_index != GAP and query[query_index] != target[target_index]:
            mismatches.add((query_index, target_index))

        path.append((query_index, target_index))
        
        row = next_row
        col = next_col
        if score_matrix[next_row, next_col] <= 0:
            break

    path = path[::-1]

    alignment = {'score': max_score,
                 'path': path,
                 'query_mappings': query_mappings,
                 'target_mappings': target_mappings,
                 'insertions': insertions,
                 'deletions': deletions,
                 'mismatches': mismatches,
                }

    return alignment

def semi_global_alignment(char* query,
                          char* target,
                          int match_bonus,
                          int mismatch_penalty,
                          int indel_penalty,
                         ):
    ''' Alignments must begin at 0, 0 and must use all of query. '''
    cdef int row, col, next_col, next_row, max_row, max_col
    cdef int match_or_mismatch, diagonal, from_left, from_above, new_score, target_index, query_index, max_score
    cdef int num_indels
    shape = (len(query) + 1, len(target) + 1)
    cdef np.ndarray[DTYPEINT_t, ndim=2] score_matrix = np.zeros(shape, DTYPEINT)
    cdef np.ndarray[DTYPEINT_t, ndim=2] row_direction_matrix = np.zeros(shape, DTYPEINT)
    cdef np.ndarray[DTYPEINT_t, ndim=2] col_direction_matrix = np.zeros(shape, DTYPEINT)
    
    cdef np.ndarray[DTYPEINT_t, ndim=1] query_mappings = np.ones(len(query), DTYPEINT) * SOFT_CLIPPED
    cdef np.ndarray[DTYPEINT_t, ndim=1] target_mappings = np.ones(len(target), DTYPEINT) * SOFT_CLIPPED

    max_score = 0

    for row in range(1, len(query) + 1):
        score_matrix[row, 0] = score_matrix[row - 1, 0] + indel_penalty
        row_direction_matrix[row, 0] = -1
    for col in range(1, len(target) + 1):
        score_matrix[0, col] = score_matrix[0, col - 1] + indel_penalty
        col_direction_matrix[0, col] = -1

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
            score_matrix[row, col] = new_score
            
            if new_score == diagonal:
                col_direction_matrix[row, col] = -1
                row_direction_matrix[row, col] = -1
            elif new_score == from_left:
                col_direction_matrix[row, col] = -1
            elif new_score == from_above:
                row_direction_matrix[row, col] = -1

    row = len(query)
    col = score_matrix[row].argmax()
    max_score = score_matrix[row, col]

    path = []
    insertions = set()
    deletions = set()
    mismatches = set()
    while True:
        next_col = col + col_direction_matrix[row, col]
        next_row = row + row_direction_matrix[row, col]
        if next_col == col:
            target_index = -1
            insertions.add(row - 1)
        else:
            target_index = col - 1
        if next_row == row:
            query_index = -1
            deletions.add(col - 1)
        else:
            query_index = row - 1
        
        if target_index != GAP:
            target_mappings[target_index] = query_index
        if query_index != GAP:
            query_mappings[query_index] = target_index
        if target_index != GAP and query_index != GAP and query[query_index] != target[target_index]:
            mismatches.add((query_index, target_index))

        path.append((query_index, target_index))
        row = next_row
        col = next_col
        if row == 0 and col == 0:
            break

    path = path[::-1]

    alignment = {'score': max_score,
                 'path': path,
                 'query_mappings': query_mappings,
                 'target_mappings': target_mappings,
                 'insertions': insertions,
                 'deletions': deletions,
                 'mismatches': mismatches,
                }

    return alignment

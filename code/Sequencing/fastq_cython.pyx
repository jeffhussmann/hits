import numpy as np
cimport cython

cdef int SANGER_OFFSET_typed = 33
SANGER_OFFSET = SANGER_OFFSET_typed

@cython.boundscheck(False)
def process_read(char* seq, char* qual, long[:, ::1] q_array, long[:, ::1] c_array):
    cdef unsigned int i, q, b, seq_length
    cdef float average_q = 0

    seq_length = len(seq)
    for i in range(seq_length):
        # Automatic type conversion means ord() is unneccesary
        q = qual[i] - SANGER_OFFSET_typed
        q_array[i, q] += 1
        average_q += q
        
        b = seq[i]
        c_array[i, b] += 1

    if seq_length != 0:
        average_q /= seq_length

    return average_q

def RSQCI_length(char *quals, int length):
    cdef int i
    for i in range(length):
        if quals[length - 1 - i] - SANGER_OFFSET_typed != 2:
            return i
    return length

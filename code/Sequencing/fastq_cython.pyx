import numpy as np
cimport numpy as np
cimport cython

DTYPEINT = np.int
ctypedef np.int_t DTYPEINT_t

SANGER_OFFSET = 33
SOLEXA_OFFSET = 64
SOLEXA_TO_SANGER_SHIFT = SOLEXA_OFFSET - SANGER_OFFSET
MAX_QUAL = 93
MAX_EXPECTED_QUAL = 41

base_order = 'TCAGN.'

def quality_and_complexity(reads, int max_read_length):
    ''' Unnecessary code duplication. '''
    cdef int SANGER_OFFSET_typed = SANGER_OFFSET
    cdef int i, q, b, read_length
    cdef char* seq
    cdef char* qual
    cdef char* bases = base_order
    cdef np.ndarray[DTYPEINT_t, ndim=2] q_array = np.zeros(shape=(max_read_length,
                                                                  MAX_EXPECTED_QUAL + 1,
                                                                 ),
                                                           dtype=DTYPEINT,
                                                          )
    cdef np.ndarray[DTYPEINT_t, ndim=2] c_array = np.zeros(shape=(max_read_length,
                                                                  256,
                                                                 ),
                                                           dtype=DTYPEINT,
                                                          )
    
    for read in reads:
        read_length = len(read.seq)
        seq = read.seq # To avoid 'Obtaining char* from temporary Python value'
        qual = read.qual
        for i in range(read_length):
            # Automatic type conversion means ord() is unneccesary
            q = qual[i] - SANGER_OFFSET_typed
            q_array[i, q] += 1
            
            b = seq[i]
            c_array[i, b] += 1
        
    # To avoid a lookup at every single base, c_array is 2*max_read_length x 256.
    # This pulls out only the columns corresponding to possible base
    # identities. 
    c_array = np.vstack([c_array.T[b] for b in bases]).T
    
    return q_array, c_array

def quality_and_complexity_paired(read_pairs, int max_read_length):
    ''' Given read_pairs, extracts the distributions of quality values seen
        and the base composition across all reads at each cycle.
        Paired iter needed to allow efficient filtering based on properties
        of both reads in pair.
    '''    
    cdef int SANGER_OFFSET_typed = SANGER_OFFSET
    cdef int i, q, b, read_length
    cdef char* seq
    cdef char* qual
    cdef char* bases = 'TCAGN.'
    cdef np.ndarray[DTYPEINT_t, ndim=2] q_array = np.zeros(shape=(2 * max_read_length,
                                                                  MAX_EXPECTED_QUAL + 1,
                                                                 ),
                                                           dtype=DTYPEINT,
                                                          )
    cdef np.ndarray[DTYPEINT_t, ndim=2] c_array = np.zeros(shape=(2 * max_read_length,
                                                                  256,
                                                                 ),
                                                           dtype=DTYPEINT,
                                                          )
    
    for (_, R1_seq, R1_qual), (_, R2_seq, R2_qual) in read_pairs:
        read_length = len(R1_seq)
        seq = R1_seq # To avoid 'Obtaining char* from temporary Python value'
        qual = R1_qual
        for i in range(read_length):
            # Automatic type conversion means ord() is unneccesary
            q = qual[i] - SANGER_OFFSET_typed
            q_array[i, q] += 1
            
            b = seq[i]
            c_array[i, b] += 1
        
        read_length = len(R2_seq)
        seq = R2_seq
        qual = R2_qual
        for i in range(read_length):
            q = qual[i] - SANGER_OFFSET_typed
            q_array[i + max_read_length, q] += 1
            
            b = seq[i]
            c_array[i + max_read_length, b] += 1
      
    # To avoid a lookup at every single base, c_array is 2*max_read_length x 256.
    # This pulls out only the columns corresponding to possible base
    # identities. 
    c_array = np.vstack([c_array.T[b] for b in bases]).T
    
    return q_array, c_array

def RSQCI_length(char *quals, int length):
    for i in range(length):
        if quals[length - 1 - i] - SANGER_OFFSET != 2:
            return i
    return length

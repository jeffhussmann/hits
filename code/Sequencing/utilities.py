from itertools import izip, islice, groupby, cycle
import subprocess
import re
import numpy as np
from Bio.Seq import Seq
try:
    import progressbar
except ImportError:
    progressbar = None

identity = lambda x: x

def reverse_complement(sequence_string):
    rc_string = str(Seq(sequence_string).reverse_complement())
    return rc_string

def complement(sequence_string):
    c_string = str(Seq(sequence_string).complement())
    return c_string

base_order = 'ACGTN-.'
base_to_index = {b: i for i, b in enumerate(base_order)}
base_to_complement_index = {b: i for i, b in enumerate(complement(base_order))}

complement_index = {i: base_to_complement_index[b] for i, b in enumerate(base_order)}

def counts_to_array(counts, dim=1):
    ''' Converts a dictionary with integer keys into an array. '''
    if dim == 1:
        if counts:
            biggest = max(counts)
            array = np.array([counts[i] for i in range(biggest + 1)])
        else:
            array = np.array([0])
    elif dim == 2:
        if counts:
            larger_of_two = map(max, counts)
            biggest = max(larger_of_two)
            array = np.zeros(shape=(biggest + 1, biggest + 1), dtype=np.int)
            for key, count in counts.items():
                array[key] = count
        else:
            array = np.array([[0]])
    return array

def mean_from_histogram(histogram):
    if histogram.sum() == 0:
        return 0
    else:
        mean = np.true_divide(np.dot(histogram, np.arange(len(histogram))), histogram.sum())
    return mean

def group_by(iterable, key=None):
    ''' Groups iterable into lists of consecutive elements that are transformed
        into the same value key.
    '''
    groups = groupby(iterable, key)
    group_lists = ((value, list(iterator)) for value, iterator in groups)
    return group_lists

def round_robin(iterables):
    ''' Modified from recipe on itertools doc page credited to George Sakkis.
    '''
    pending = len(iterables)
    nexts = cycle(iter(it).next for it in iterables)
    while pending:
        try:
            for next in nexts:
                yield next()
        except StopIteration:
            pending -= 1
            # The next pending yields from nexts will be the remaining
            # non-empty iterators.
            nexts = cycle(islice(nexts, None, pending))

def line_count(file_name):
    ''' Return the number of lines in file_name by calling wc. '''
    wc_output = subprocess.check_output(['wc', '-l', file_name])
    count, _ = wc_output.strip().split()
    count = int(count)
    return count

def progress_bar(max_val, iterable=None):
    if progressbar == None:
        # If the module wasn't imported
        if iterable:
            return iterable
        else:
            return xrange(max_val)
    else:
        max_str = str(len(str(max_val)))
        format_string = '%(value)' + max_str + 'd / %(max)d'
        widgets = [progressbar.Bar('='),
                   ' ',
                   progressbar.FormatLabel(format_string),
                   ' ',
                   progressbar.ETA(),
                  ]
        bar = progressbar.ProgressBar(widgets=widgets, maxval=max_val)
        if iterable:
            return bar(iterable)
        else:
            return bar(xrange(max_val))

def file_line_groups(name_number_pairs):
    """ Given a list of (file_name, number) pairs, yields tuple
    of groups of (number) lines with newlines stripped from each (file_name). """
    
    numbers = [number for file_name, number in name_number_pairs]
    fhs = [open(file_name) for file_name_number in name_number_pairs]
    
    fh_groups = [izip(*[fh]*number) for fh, number in zip(fhs, numbers)]
    for i, groups in enumerate(izip(*fh_groups)):
        groups = [[line.rstrip() for line in group] for group in groups]
        yield groups

def hamming_distance(s1, s2):
    return sum(1 for c1, c2 in izip(s1, s2) if c1 != c2)

def pairwise(s):
    """ Returns the elements of s in overlapping pairs. """
    return [(s[i - 1], s[i]) for i in range(1, len(s))]

def all_consecutive(s):
    """ Returns True if the elements of s when sorted are consecutive integers. """
    for (x, y) in pairwise(sorted(s)):
        if y - x != 1:
            return False
    return True

def decompose_homopolymer_sequence(fragment):
    """ Decomposes a string into sequences of homopolymer chars and counts."""
    sequence = []
    counts = []
    current_char = fragment[0]
    current_length = 1
    for i in range(1, len(fragment)):
        if fragment[i] == current_char:
            current_length += 1
        else:
            sequence.append(current_char)
            counts.append(current_length)
            current_char = fragment[i]
            current_length = 1
    sequence.append(current_char)
    counts.append(current_length)
    return (sequence, counts)

def reconstitute_homopolymer_sequence(sequence, counts):
    """ Expands sequences of homopolymer chars and counts into a string."""
    reconstituted = []
    for letter, count in zip(sequence, counts):
        reconstituted.append(letter*count)
    return "".join(reconstituted)

def process_annotations(read):
    # flow_indexes is a list of indexes of flows corresponding to each called base
    flow_indexes = np.cumsum(np.array(read.annotations['flow_index']))
    # clip_qual_left is (1-based in SFF but converted to 0-based by BioPython) 
    # the position of the first base after the clipping point.
    # Subtract 1 to get the index of the last base in the tag, 
    # look up in flow_indexes, subtract 1 to get 0-based flow index.
    read.annotations['last_flow_in_tag'] = flow_indexes[read.annotations['clip_qual_left'] - 1] - 1 
    flows = map(lambda x: x / 100., read.annotations['flow_values'])
    flows = np.array(flows)
    # Remove the 1 nuc from the tag from the first hp (TODO: do some tag+adapters end in a non-1 hp?)
    flows[read.annotations['last_flow_in_tag']] -= 1
    read.annotations['flows'] = flows
    read.annotations['last_overall_flow'] = flow_indexes[-1] - 1 
    # The NCBI format specification says clip_qual_left is the 1-based index of the last good base,
    # BioPython gives the 0-based index of the first bad base, which is the same thing,
    # but inconsistent with the spirit of how it handles clip_qual_left.
    read.annotations['last_good_flow'] = flow_indexes[read.annotations['clip_qual_right'] - 1] - 1
    read.annotations['last_good_base'] = read.annotations['clip_qual_right'] - 1 

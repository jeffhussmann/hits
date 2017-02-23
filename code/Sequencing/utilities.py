from __future__ import division
from itertools import izip, islice, groupby, cycle, product
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

base_order = 'ACGTN-'
base_to_index = {b: i for i, b in enumerate(base_order)}
base_to_complement_index = {b: i for i, b in enumerate(complement(base_order))}

complement_index = {i: base_to_complement_index[b] for i, b in enumerate(base_order)}

IUPAC = {'A': {'A'},
         'C': {'C'},
         'G': {'G'},
         'T': {'T'},
         'M': {'A', 'C'},
         'R': {'A', 'G'},
         'W': {'A', 'T'},
         'S': {'C', 'G'},
         'Y': {'C', 'T'},
         'K': {'G', 'T'},
         'V': {'A', 'C', 'G'},
         'H': {'A', 'C', 'T'},
         'D': {'A', 'G', 'T'},
         'B': {'C', 'G', 'T'},
         'N': {'G', 'A', 'T', 'C'},
        }

def counts_to_array(counts, dim=1):
    ''' Converts a dictionary with integer keys into an array. '''
    if dim == 1:
        if len(counts) > 0:
            biggest = max(counts)
            array = np.array([counts[i] for i in range(biggest + 1)])
        else:
            array = np.array([0])
    elif dim == 2:
        if len(counts) > 0:
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

def progress_bar(iterable, max_val=None):
    if progressbar == None:
        # If the module wasn't imported
        return iterable
    else:
        if max_val is None:
            max_val = len(iterable)

        max_str = str(len(str(max_val)))
        format_string = '%(value)' + max_str + 'd / %(max)d'
        widgets = [progressbar.Bar('='),
                   ' ',
                   progressbar.FormatLabel(format_string),
                   ' ',
                   progressbar.ETA(),
                  ]
        bar = progressbar.ProgressBar(widgets=widgets, maxval=max_val)
        return bar(iterable)

def pairwise(s):
    """ Returns the elements of s in overlapping pairs. """
    return [(s[i - 1], s[i]) for i in range(1, len(s))]

def all_consecutive(s):
    """ Returns True if the elements of s when sorted are consecutive integers. """
    for (x, y) in pairwise(sorted(s)):
        if y - x != 1:
            return False
    return True

def empirical_cdf(values):
    ''' From stackoverflow. '''
    sorted_values = np.sort(values)
    cumulative = np.true_divide(np.arange(len(sorted_values)), len(sorted_values))
    return sorted_values, cumulative

def mers(k, include_N=False):
    if include_N:
        chars = 'TCAGN'
    else:
        chars = 'TCAG'
    return (''.join(mer) for mer in product(chars, repeat=k))

def smooth(ys, window):
    smoothed = ys.astype(float)
    for i in range(window, len(ys) - window):
        smoothed[i] = sum(ys[i - window:i + window + 1]) / float(2 * window + 1)
    return smoothed

def reverse_dictionary(d):
    r = {v: k for k, v in d.iteritems()}
    return r

def split_nonempty(string, delim):
    return [f for f in string.split(delim) if f]

from collections import defaultdict
from numbers import Number

from . import sam

def are_disjoint(first, second):
    if first.is_empty or second.is_empty:
        return True
    else:
        return first.start > second.end or second.start > first.end

def are_adjacent(first, second):
    return first.start == second.end + 1 or second.start == first.end + 1

class Interval(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.is_empty = (end < start)

    @classmethod
    def from_feature(self, feature):
        ''' only really necessary to set is_empty '''
        return Interval(feature.start, feature.end)

    @classmethod
    def empty(self):
        return Interval(-1, -2)
    
    @classmethod
    def from_slice(self, sl):
        if sl.start == None:
            start = 0
        else:
            start = sl.start

        if sl.stop == None:
            end = np.inf
        else:
            end = sl.stop - 1 # Note the -1

        return Interval(start, end)

    def __or__(self, other):
        if are_disjoint(self, other):
            left, right = sorted([self, other])
            if are_adjacent(self, other):
                intervals = [Interval(left.start, right.end)]
            else:
                intervals = [left, right]
        else:
            intervals = [Interval(min(self.start, other.start), max(self.end, other.end))]
            
        return DisjointIntervals(intervals)
    
    def __and__(self, other):
        if isinstance(other, DisjointIntervals):
            # Defer to definition in DisjointIntervals
            return other & self
        elif are_disjoint(self, other):
            return []
        else:
            return Interval(max(self.start, other.start), min(self.end, other.end))
        
    def __contains__(self, other):
        if isinstance(other, Interval):
            # is a strict sub-interval of
            return (other.start >= self.start and other.end <= self.end) and (self != other)
        elif isinstance(other, Number):
            return self.start <= other <= self.end
        else:
            raise ValueError(other)
        
    @property    
    def comparison_key(self):
        return self.start, self.end

    @property
    def total_length(self):
        return len(self)
    
    def __lt__(self, other):
        return self.comparison_key < other.comparison_key
    
    def __repr__(self):
        return '[{0:,} - {1:,}]'.format(self.start, self.end)
    
    def __key(self):
        return (self.start, self.end)

    def __eq__(self, other):
        return self.__key() == other.__key()

    def __hash__(self):
        return hash(self.__key())
    
    def __ne__(self, other):
        return not self == other

    def __len__(self):
        if self.is_empty:
            return 0
        else:
            return self.end - self.start + 1

    def __sub__(self, other):
        if isinstance(other, DisjointIntervals):
            survived_all = DisjointIntervals([self])
            for other_interval in other:
                survived_this = self - other_interval
                survived_all = survived_all & survived_this
            return survived_all

        else:
            left = Interval(self.start, min(self.end, other.start - 1))
            right = Interval(max(self.start, other.end + 1), self.end)
            disjoint = DisjointIntervals([left, right])

            if disjoint.total_length == 0:
                return Interval(-1, -2)
            elif len(disjoint.intervals) == 1:
                return disjoint.intervals[0]
            else:
                return disjoint
    
class DisjointIntervals(object):
    def __init__(self, intervals):
        self.intervals = sorted([i for i in intervals if i.end >= i.start])
        self.is_empty = (len(self.intervals) == 0)
        
    def __len__(self):
        return len(self.intervals)
    
    @property
    def start(self):
        if len(self.intervals) == 0:
            return None
        else:
            return min(interval.start for interval in self.intervals)
    
    @property
    def end(self):
        if len(self.intervals) == 0:
            return None
        else:
            return max(interval.end for interval in self.intervals)
    
    def __repr__(self):
        return '{{{}}}'.format(', '.join(map(str, self.intervals)))
    
    def __getitem__(self, sl):
        return self.intervals[sl]
    
    def __or__(self, other):
        if isinstance(other, DisjointIntervals):
            everything = self
            for other_interval in other:
                everything = everything | other_interval
            return everything
        else:
            disjoints = []
            
            for interval in self.intervals:
                union = interval | other
                if len(union) > 1:
                    disjoints.append(interval)
                else:
                    other = union[0]
                    
            disjoints.append(other)
            
            return DisjointIntervals(sorted(disjoints))
    
    def __and__(self, other):
        if isinstance(other, DisjointIntervals):
            survived_some = DisjointIntervals([])
            for other_interval in other:
                survived_this = self & other_interval
                survived_some = survived_some | survived_this
            return survived_some
        else:
            intersections = []
            for interval in self.intervals:
                intersection = interval & other
                if intersection:
                    intersections.append(intersection)
                
            return DisjointIntervals(intersections)
    
    def __eq__(self, other):
        return self.intervals == other.intervals

    def __iter__(self):
        return iter(self.intervals)

    def __hash__(self):
        return hash(self.intervals)
    
    def __ne__(self, other):
        return not self == other
    
    def __contains__(self, other):
        if isinstance(other, Number):
            other = Interval(other, other)

        return (self | other) == self

    @property
    def total_length(self):
        return sum(len(i) for i in self.intervals)
    
    def __sub__(self, other):
        if isinstance(other, DisjointIntervals):
            raise NotImplementedError
        else:
            pieces = []
            for interval in self.intervals:
                piece = interval - other
                if isinstance(piece, Interval):
                    pieces.append(piece)
                elif isinstance(piece, DisjointIntervals):
                    pieces.extend(piece)

            return DisjointIntervals(pieces)

def get_covered(alignment):
    if alignment is None or alignment.is_unmapped:
        return Interval(-1, -2)
    else:
        return Interval(*sam.query_interval(alignment))

def make_disjoint(intervals):
    disjoint = DisjointIntervals([])
    for interval in intervals:
        disjoint = disjoint | interval
    return disjoint

def get_disjoint_covered(alignments):
    intervals = [get_covered(al) for al in alignments if al is not None]
    covered = make_disjoint(intervals)
    return covered

def remove_nested(alignments):
    unnecessary = set()
    covered_list = [get_covered(al) for al in alignments]
    for i, left in enumerate(covered_list):
        for j, right in enumerate(covered_list):
            if i == j:
                continue
            if left in right:
                unnecessary.add(i)
    necessary = [al for i, al in enumerate(alignments) if i not in unnecessary]
    return necessary

def make_parsimonious(alignments):
    initial_covered = get_disjoint_covered(alignments)
    
    no_nested = remove_nested(alignments)
    interval_to_als = defaultdict(list)
    for al in no_nested:
        interval_to_als[get_covered(al)].append(al)
        
    unique_intervals = sorted(interval_to_als, key=len, reverse=True)
    remaining = unique_intervals
        
    contributes = []
    for possibly_exclude in unique_intervals:
        exclude_one = [intvl for intvl in remaining if intvl != possibly_exclude]
        now_covered = make_disjoint(exclude_one)
        if initial_covered != now_covered:
            contributes.append(possibly_exclude)
        else:
            remaining = exclude_one
            
    parsimonious = []
    for interval in contributes:
        parsimonious.extend(interval_to_als[interval])
    
    if get_disjoint_covered(parsimonious) != initial_covered:
        raise ValueError
        
    return parsimonious

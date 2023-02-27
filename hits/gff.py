import pprint
import urllib

from hits import utilities

def parse_attribute_string(attribute_string):
    if attribute_string == '.':
        parsed = {}
    else:
        fields = attribute_string.split(';')
        pairs = [field.split('=') for field in fields]
        parsed = {name: urllib.parse.unquote(value).strip('"') for name, value in pairs}

    return parsed

def make_attribute_string(attribute):
    entries = []
    for key, value in sorted(attribute.items()):
        key = urllib.parse.quote(str(key), safe='')
        value = urllib.parse.quote(str(value), safe='')
        entry = f'{key}={value}'
        entries.append(entry)

    attribute_string = ';'.join(entries)
    return attribute_string

class Feature(object):
    def __init__(self, line=None):
        if line == None:
            # Allow __init__ to be called with no arguments to allow the
            # @classmethod constructor below.
            return

        fields = line.strip().split('\t')
        
        self.seqname = fields[0]
        self.source = fields[1]
        self.feature = fields[2]

        # Note conversion to 0-based indexing.
        self.start = int(fields[3]) - 1
        self.end = int(fields[4]) - 1

        self.score = fields[5]
        self.strand = fields[6]
        self.frame = fields[7]
        if self.frame != '.':
            self.frame = int(self.frame)
        self.attribute_string = fields[8]
        
        self.parse_attribute_string()

        self.parent = None
        self.children = set()
    
    @classmethod
    def from_fields(cls,
                    seqname='.',
                    source='.',
                    feature='.',
                    start='.',
                    end='.',
                    score='.',
                    strand='.',
                    frame='.',
                    attribute_string='.',
                    ID=None,
                   ):
        obj = cls()
        obj.seqname = seqname
        obj.source = source
        obj.feature = feature
        obj.start = start
        obj.end = end
        obj.score = score
        obj.strand = strand
        obj.frame = frame
        obj.attribute_string = attribute_string
        obj.parse_attribute_string()
        if ID is not None:
            obj.attribute['ID'] = ID
        return obj

    @classmethod
    def sequence_edge(cls, seqname, position):
        obj = cls()
        obj.seqname = seqname
        obj.feature = 'edge'
        obj.start = position
        obj.end = position
        obj.strand = '.'
        obj.source = '.'
        obj.score = '.'
        obj.frame = '.'
        obj.attribute = {'ID': f'edge_{seqname}_{position}'}
        obj.populate_attribute_string()
        return obj

    def parse_attribute_string(self):
        self.attribute = parse_attribute_string(self.attribute_string)

    def populate_attribute_string(self):
        self.attribute_string = make_attribute_string(self.attribute)

    def populate_connections(self, id_to_object):
        parent_id = self.attribute.get('Parent')
        if parent_id:
            self.parent = id_to_object[parent_id]
            self.parent.children.add(self)

    @property
    def ID(self):
        return self.attribute.get('ID')

    @property
    def descendants(self):
        descendants = set(self.children)

        for child in self.children:
            descendants.update(child.descendants)

        return descendants

    def print_family(self, level=0):
        indent = '\t' * min(level, 3)

        print(f'{indent}{self}')

        if level == 0:
            pprint.pprint(self.attribute)

        for child in sorted(self.children):
            child.print_family(level=level + 1)
    
    def __str__(self):
        self.populate_attribute_string()
        fields = (self.seqname,
                  self.source,
                  self.feature,
                  str(self.start + 1), # Note conversion back to 1-based indexing.
                  str(self.end + 1),
                  self.score,
                  self.strand,
                  str(self.frame),
                  self.attribute_string,
                 )
        line = '\t'.join(fields)
        return line
    
    @property
    def pasteable(self):
        return f'{self.seqname}:{self.start}-{self.end}'
    
    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)

    @property
    def comparison_key(self):
        key = (self.seqname,
               self.start,
               self.end,
               self.feature,
               self.strand,
              )
        return key

    def __lt__(self, other):
        return self.comparison_key < other.comparison_key

    def is_contained_in(self, other):
        return self.seqname == other.seqname and \
               self.start >= other.start and \
               self.end <= other.end

    def __len__(self):
        return self.end - self.start + 1

    def sequence(self, reference_sequences):
        ''' Given dictionary of reference_sequences which includes an entry for 
        self.seqname, returns stranded sequence for this feature.
        '''
        seq = reference_sequences[self.seqname][self.start:self.end + 1]
        if self.strand == '-':
            seq = utilities.reverse_complement(seq)

        return seq

def populate_all_connections(features):
    for f in features:
        f.children = set()
        f.parent = None

    id_to_object = {f.attribute['ID']: f for f in features if 'ID' in f.attribute}

    for f in features:
        f.populate_connections(id_to_object)

def get_all_features(gff_fn, populate_connections=True):
    def relevant_lines(gff_fn):
        # Ignore any line starting with '#' and all lines after any line
        # starting with '##FASTA'
        with open(gff_fn) as gff_fh:
            for line in gff_fh:
                if line.startswith('##FASTA'):
                    break
                elif line.startswith('#'):
                    continue
                else:
                    yield line

    all_features = [Feature(line) for line in relevant_lines(gff_fn)]

    if populate_connections:
        populate_all_connections(all_features)

    return all_features

def get_top_level_features(features):
    populate_all_connections(features)
    top_level_features = [f for f in features if f.parent == None]
    return top_level_features

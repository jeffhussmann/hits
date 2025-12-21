from Bio.SeqFeature import SeqFeature, FeatureLocation

import hits.utilities

class Enzyme:
    def __init__(self, name, recognition_site, cut_offsets):
        self.name = name
        self.recognition_site = recognition_site
        self.cut_offsets = cut_offsets

        self.overhang = self.cut_offsets['+'] - self.cut_offsets['-']
        self.overhang_length = abs(self.overhang)

    def forward_matches(self, seq):
        return hits.utilities.find_all_substring_starts(seq.upper(), self.recognition_site)

    def reverse_matches(self, seq):
        return hits.utilities.find_all_substring_starts(seq.upper(), hits.utilities.reverse_complement(self.recognition_site))

    def total_matches(self, seq):
        return len(self.forward_matches(seq)) + len(self.reverse_matches(seq))

    def cut_afters(self, seq):
        cut_afters = {
            '+': set(),
            '-': set(),
        }

        for forward_match in self.forward_matches(seq):
            cut_afters['+'].add(forward_match + self.cut_offsets['+'])
            cut_afters['-'].add(forward_match + self.cut_offsets['-'])
            
        for reverse_match in self.reverse_matches(seq):
            cut_afters['-'].add(reverse_match + len(self.recognition_site) - 1 - self.cut_offsets['+'] - 1)
            cut_afters['+'].add(reverse_match + len(self.recognition_site) - 1 - self.cut_offsets['-'] - 1)

        return cut_afters

    def overhang_features(self, seq):
        ''' Return SeqFeatures of overhangs produced by digestion,
        annotated on the strand retained by side not containing the
        recognition site (i.e. the side retained in a IIS cloning
        strategy).
        '''
        all_overhangs = []
        features = []

        if self.overhang_length > 0:

            for forward_match in self.forward_matches(seq):
                cuts = (
                    forward_match + self.cut_offsets['+'],
                    forward_match + self.cut_offsets['-'],
                )

                if self.overhang < 0:
                    strand = 1
                else:
                    strand = -1

                overhang = (*cuts, strand)

                all_overhangs.append(overhang)
                
            for reverse_match in self.reverse_matches(seq):
                cuts = (
                    reverse_match + len(self.recognition_site) - 1 - self.cut_offsets['+'] - 1,
                    reverse_match + len(self.recognition_site) - 1 - self.cut_offsets['-'] - 1,
                )

                if self.overhang > 0:
                    strand = 1
                else:
                    strand = -1

                overhang = (*cuts, strand)

                all_overhangs.append(overhang)

            for plus_cut_after, minus_cut_after, strand in all_overhangs:
                left_cut_after, right_cut_after = sorted([plus_cut_after, minus_cut_after])
                start = left_cut_after + 1 # by definition of cut_after
                end = right_cut_after + 1 # because FeatureLocation is end-exclusive

                overhang_seq = seq[start:end]

                if strand == -1:
                    overhang_seq = hits.utilities.reverse_complement(overhang_seq)

                name = f'overhang_{overhang_seq}'
                feature = SeqFeature(location=FeatureLocation(start,
                                                              end,
                                                              strand=strand,
                                                             ),
                                     id=name,
                                     type='misc_feature',
                                     qualifiers={
                                         'label': name,
                                         'overhang_seq': overhang_seq,
                                         'overhang_seq_rc': hits.utilities.reverse_complement(overhang_seq),
                                     },
                                    )
                features.append(feature)

        return features

    def digest(self, seq):
        ''' If seq contains exactly one forward match and exactly one reverse match,
            return 3 digested fragments each containing all full overhangs.
            If seq contains exactly one total forward or reverse match,
            middle fragment will be empty.
        '''
        forward_matches = self.forward_matches(seq)
        reverse_matches = self.reverse_matches(seq)

        # Last condition is for palindromic sites
        if len(forward_matches) == 1 and len(reverse_matches) == 1 and forward_matches[0] != reverse_matches[0]:
            forward_match = forward_matches[0]
            reverse_match = reverse_matches[0]

            if forward_match > reverse_match:
                right_start = forward_match + min(self.cut_offsets.values()) + 1
                right_seq = seq[right_start:]

                left_end = reverse_match + len(self.recognition_site) - 1 - min(self.cut_offsets.values())
                left_seq = seq[:left_end]

                middle_end = right_start + abs(self.overhang_length)
                middle_start = left_end - abs(self.overhang_length)
                middle_seq = seq[middle_start:middle_end]

            else:
                left_end = forward_match + max(self.cut_offsets.values()) + 1
                left_seq = seq[:left_end]

                right_start = reverse_match + len(self.recognition_site) - 1 - max(self.cut_offsets.values())
                right_seq = seq[right_start:]

                middle_end = right_start + abs(self.overhang_length)
                middle_start = left_end - abs(self.overhang_length)
                middle_seq = seq[middle_start:middle_end]

        elif (((len(forward_matches), len(reverse_matches)) in {(1, 0), (0, 1)}) or
              (len(forward_matches) == 1 and len(reverse_matches) == 1 and forward_matches[0] == reverse_matches[0])
        ):

            all_cut_afters = set.union(*self.cut_afters(seq).values())
            left_end = max(all_cut_afters) + 1
            left_seq = seq[:left_end]

            right_start = min(all_cut_afters) + 1
            right_seq = seq[right_start:]

            middle_seq = ''

        else:
            raise ValueError('unexpected number of restriction sites')

        return left_seq, middle_seq, right_seq

enzymes = {
    'BsmBI': Enzyme(
        'BsmBI',
        'CGTCTC',
        {
            '+': 6,
            '-': 10,
        },
    ),
    'NruI': Enzyme(
        'NruI',
        'TCGCGA',
        {
            '+': 2,
            '-': 2,
        },
    ),
    'XhoI': Enzyme(
        'XhoI',
        'CTCGAG',
        {
            '+': 0,
            '-': 4,
        },
    ),
    'BbsI': Enzyme(
        'BbsI',
        'GAAGAC',
        {
            '+': 7,
            '-': 11,
        },
    ),
}

def enzyme_from_possible_string(possible_string):
    if isinstance(possible_string, str):
        enzyme = enzymes[possible_string]
    else:
        enzyme = possible_string

    if not isinstance(enzyme, Enzyme):
        raise ValueError(enzyme)

    return enzyme

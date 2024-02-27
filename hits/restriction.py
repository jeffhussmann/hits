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

    def digest(self, seq):
        ''' If seq contains exactly one forward match and exactly one reverse match,
            return 3 digested fragments each containing all full overhangs.
            If seq contains exactly one total forward or reverse match,
            middle fragment will be empty.
        '''
        forward_matches = self.forward_matches(seq)
        reverse_matches = self.reverse_matches(seq)

        if len(forward_matches) == 1 and len(reverse_matches) == 1:
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

        elif (len(forward_matches), len(reverse_matches)) in {(1, 0), (0, 1)}:
            if len(forward_matches) == 0:
                match = forward_matches[0]
            else:
                match = reverse_matches[0]

            left_end = match + max(self.cut_offsets.values()) + 1
            left_seq = seq[:left_end]

            right_start = match + len(self.recognition_site) - 1 - max(self.cut_offsets.values())
            right_seq = seq[right_start:]

            middle_seq = ''

        return left_seq, middle_seq, right_seq

    def assemble(self, vector, oligo, digest_oligo=True):
        left_vector, _, right_vector = self.digest(vector)

        if digest_oligo:
            _, middle_oligo, _ = self.digest(oligo)
        else:
            middle_oligo = oligo

        if self.overhang_length > 0:
            # Golden gate
            if left_vector[-self.overhang_length:] != middle_oligo[:self.overhang_length]:
                raise ValueError('incompatible overhang on left')

            if middle_oligo[-self.overhang_length:] != right_vector[:self.overhang_length]:
                raise ValueError('incompatible overhang on right')

            assembled = left_vector[:-self.overhang_length] + middle_oligo + right_vector[self.overhang_length:]

        else:
            # Gibson
            def find_HA_length(left, right):
                matching_lengths = [HA_length for HA_length in range(10, min(len(left), len(right))) if left[-HA_length:] == right[:HA_length]]
                if len(matching_lengths) == 0:
                    raise ValueError
                else:
                    return max(matching_lengths)

            left_HA = find_HA_length(left_vector, oligo)
            right_HA = find_HA_length(oligo, right_vector)

            assembled = left_vector[:-left_HA] + oligo + right_vector[right_HA:]

        return assembled

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
}

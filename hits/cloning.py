import hits.sw
import hits.restriction
import hits.utilities

def gibson_assemble(PCR_template, PCR_primers, backbone, restriction_enzyme):
    restriction_enzyme = hits.restriction.enzyme_from_possible_string(restriction_enzyme)

    PCR_product = hits.sw.amplify_sequence_with_primers(PCR_template, PCR_primers, both_orientations=True)

    cut_afters = restriction_enzyme.cut_afters(backbone)

    # for leftmost, go up to cut on top strand.
    last_three_prime_on_left = min(cut_afters['+'])
    left_backbone = backbone[:last_three_prime_on_left + 1]

    # find last place in left that 15mer prefix of insert matches
    if PCR_product[:15] not in left_backbone:
        PCR_product = hits.utilities.reverse_complement(PCR_product)
        if PCR_product[:15] not in left_backbone:
            raise ValueError('homology arm not detected')

    left_HA_start = left_backbone.index(PCR_product[:15])
    if not PCR_product.startswith(left_backbone[left_HA_start:]):
        raise ValueError('homology arm doesn\'t match')
    
    # for rightmost, go from cut on bottom strand
    first_three_prime_on_right = max(cut_afters['-']) + 1
    right_backbone = backbone[first_three_prime_on_right:]

    # find first place in right that 15mer suffix of insert matches
    if PCR_product[-15:] not in right_backbone:
        raise ValueError('homology arm not detected')

    right_HA_end = right_backbone.index(PCR_product[-15:]) + len(PCR_product[-15:])
    if not PCR_product.endswith(right_backbone[:right_HA_end]):
        raise ValueError('homology arm doesn\'t match')

    return left_backbone[:left_HA_start] + PCR_product + right_backbone[right_HA_end:]

def goldengate_assemble(backbone, insert, restriction_enzyme='BsmBI', digest_insert=True):
    left_vector, _, right_vector, left_vector_overhang_length, right_vector_overhang_length = possibly_double_digest(backbone, restriction_enzyme)

    if digest_insert:
        _, middle_insert, _, left_insert_overhang_length, right_insert_overhang_length = possibly_double_digest(insert, restriction_enzyme)
    else:
        middle_insert, left_insert_overhang_length, right_insert_overhang_length = insert

    if left_vector_overhang_length == 0 or right_vector_overhang_length == 0:
        raise ValueError('no overhang')

    def attempt_assembly(flip_insert):
        if flip_insert:
            oriented_insert = hits.utilities.reverse_complement(middle_insert)
            oriented_left_insert_overhang_length, oriented_right_insert_overhang_length = right_insert_overhang_length, left_insert_overhang_length
        else:
            oriented_insert = middle_insert
            oriented_left_insert_overhang_length, oriented_right_insert_overhang_length = left_insert_overhang_length, right_insert_overhang_length

        left_compatible = (left_vector[-left_vector_overhang_length:] == oriented_insert[:oriented_left_insert_overhang_length])
        right_compatible = (oriented_insert[-oriented_right_insert_overhang_length:] == right_vector[:right_vector_overhang_length])

        if left_compatible and right_compatible:
            assembled = left_vector[:-left_vector_overhang_length] + oriented_insert + right_vector[right_vector_overhang_length:]
        else:
            assembled = None

        return assembled

    assemblies = []

    for flip_insert in [True, False]:
        assembled = attempt_assembly(flip_insert)
        if assembled is not None:
            assemblies.append(assembled)

    if len(assemblies) == 0:
        raise ValueError('no compatible overhang orientation')

    elif len(assemblies) == 2:
        raise ValueError('ambiguous overhang orientation')

    else:
        assembled = assemblies[0]

    return assembled

def double_digest(sequence, enzymes):
    first_enzyme, second_enzyme = [hits.restriction.enzyme_from_possible_string(s) for s in enzymes]

    left, middle, right = first_enzyme.digest(sequence)
    if middle != '':
        raise ValueError

    try:
        second_left, second_middle, second_right = second_enzyme.digest(left)
        right_enzyme, left_enzyme = first_enzyme, second_enzyme
        left, middle, right = [second_left, second_right, right]
        
    except ValueError:
        second_left, second_middle, second_right = second_enzyme.digest(right)
        left_enzyme, right_enzyme = first_enzyme, second_enzyme
        left, middle, right = [left, second_left, second_right]

    if second_middle != '':
        raise ValueError
    
    return left, middle, right, left_enzyme.overhang_length, right_enzyme.overhang_length

def possibly_double_digest(sequence, enzymes):
    if isinstance(enzymes, (tuple, list)):
        if len(enzymes) != 2:
            raise ValueError

        left, middle, right, left_overhang_length, right_overhang_length = double_digest(sequence, enzymes)

    else:
        enzyme = hits.restriction.enzyme_from_possible_string(enzymes)

        left, middle, right = enzyme.digest(sequence)
        left_overhang_length, right_overhang_length = enzyme.overhang_length, enzyme.overhang_length

    return left, middle, right, left_overhang_length, right_overhang_length
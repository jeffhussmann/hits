import hits.sw
import hits.restriction
import hits.utilities

def gibson_assemble(PCR_template, PCR_primers, backbone, restriction_enzyme):
    if isinstance(restriction_enzyme, str):
        restriction_enzyme = hits.restriction.enzymes[restriction_enzyme]

    PCR_product = hits.sw.amplify_sequence_with_primers(PCR_template, PCR_primers)

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
    if isinstance(restriction_enzyme, str):
        restriction_enzyme = hits.restriction.enzymes[restriction_enzyme]

    left_vector, _, right_vector = restriction_enzyme.digest(backbone)

    if digest_insert:
        _, middle_insert, _ = restriction_enzyme.digest(insert)
    else:
        middle_insert = insert

    if restriction_enzyme.overhang_length == 0:
        raise ValueError('no overhang')

    if left_vector[-restriction_enzyme.overhang_length:] != middle_insert[:restriction_enzyme.overhang_length]:
        raise ValueError('incompatible overhang on left')

    if middle_insert[-restriction_enzyme.overhang_length:] != right_vector[:restriction_enzyme.overhang_length]:
        raise ValueError('incompatible overhang on right')

    assembled = left_vector[:-restriction_enzyme.overhang_length] + middle_insert + right_vector[restriction_enzyme.overhang_length:]

    return assembled
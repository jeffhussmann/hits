import random

import pysam
import hits.sam as sam

MATCH = sam.BAM_CMATCH
DEL = sam.BAM_CDEL
INS = sam.BAM_CINS
SOFT_CLIP = sam.BAM_CSOFT_CLIP

# Cropping

def test_simple():
    al = pysam.AlignedSegment()
    al.reference_start = 100
    al.cigar = [(MATCH, 100)]
    al.seq = 'A' * 100

    cropped = pysam.AlignedSegment()
    cropped.reference_start = 150
    cropped.cigar = [(SOFT_CLIP, 50), (MATCH, 10), (SOFT_CLIP, 40)]
    cropped.seq = 'A' * 100

    output = sam.crop_al_to_ref_int(al, 150, 159)
    assert (output == cropped)

    output = sam.crop_al_to_ref_int(al, 150, 158)
    assert (output != cropped)

def test_ends_in_deletion():
    al = pysam.AlignedSegment()
    al.reference_start = 100
    al.cigar = [(MATCH, 50), (DEL, 10), (MATCH, 50)]
    al.seq = 'A' * 100
    
    cropped = pysam.AlignedSegment()
    cropped.reference_start = 100
    cropped.cigar = [(MATCH, 50), (SOFT_CLIP, 50)]
    cropped.seq = 'A' * 100
    
    output = sam.crop_al_to_ref_int(al, 100, 155)
    assert (output == cropped)
    
    output = sam.crop_al_to_ref_int(al, 100, 160)
    assert (output != cropped)

def test_ends_just_after_deletion():
    al = pysam.AlignedSegment()
    al.reference_start = 100
    al.cigar = [(MATCH, 50), (DEL, 10), (MATCH, 50)]
    al.seq = 'A' * 100
    
    cropped = pysam.AlignedSegment()
    cropped.reference_start = 100
    cropped.cigar = [(MATCH, 50), (DEL, 10), (MATCH, 1), (SOFT_CLIP, 49)]
    cropped.seq = 'A' * 100
    
    output = sam.crop_al_to_ref_int(al, 100, 160)
    assert (output == cropped)
    
    output = sam.crop_al_to_ref_int(al, 100, 155)
    assert (output != cropped)

def test_starts_in_deletion():
    al = pysam.AlignedSegment()
    al.reference_start = 100
    al.cigar = [(MATCH, 50), (DEL, 10), (MATCH, 50)]
    al.seq = 'A' * 100
    
    cropped = pysam.AlignedSegment()
    cropped.reference_start = 160
    cropped.cigar = [(SOFT_CLIP, 50), (MATCH, 50)]
    cropped.seq = 'A' * 100
    
    output = sam.crop_al_to_ref_int(al, 155, 210)
    assert (output == cropped)
    
    output = sam.crop_al_to_ref_int(al, 149, 210)
    assert (output != cropped)

def test_starts_just_before_deletion():
    al = pysam.AlignedSegment()
    al.reference_start = 100
    al.cigar = [(MATCH, 50), (DEL, 10), (MATCH, 50)]
    al.seq = 'A' * 100
    
    cropped = pysam.AlignedSegment()
    cropped.reference_start = 149
    cropped.cigar = [(SOFT_CLIP, 49), (MATCH, 1), (DEL, 10), (MATCH, 50)]
    cropped.seq = 'A' * 100
    
    output = sam.crop_al_to_ref_int(al, 149, 210)
    assert (output == cropped)
    
    output = sam.crop_al_to_ref_int(al, 150, 210)
    assert (output != cropped)

def test_ends_just_before_insertion():
    al = pysam.AlignedSegment()
    al.reference_start = 100
    al.cigar = [(MATCH, 50), (INS, 10), (MATCH, 50)]
    al.seq = 'A' * 110
    
    cropped = pysam.AlignedSegment()
    cropped.reference_start = 100
    cropped.cigar = [(MATCH, 50), (SOFT_CLIP, 60)]
    cropped.seq = 'A' * 110
    
    output = sam.crop_al_to_ref_int(al, 100, 149)
    assert (output == cropped)
    
    output = sam.crop_al_to_ref_int(al, 100, 150)
    assert (output != cropped)

def test_ends_just_after_insertion():
    al = pysam.AlignedSegment()
    al.reference_start = 100
    al.cigar = [(MATCH, 50), (INS, 10), (MATCH, 50)]
    al.seq = 'A' * 110
    
    cropped = pysam.AlignedSegment()
    cropped.reference_start = 100
    cropped.cigar = [(MATCH, 50), (INS, 10), (MATCH, 1), (SOFT_CLIP, 49)]
    cropped.seq = 'A' * 110
    
    output = sam.crop_al_to_ref_int(al, 100, 150)
    assert (output == cropped)
    
    output = sam.crop_al_to_ref_int(al, 100, 149)
    assert (output != cropped)

def test_starts_just_before_insertion():
    al = pysam.AlignedSegment()
    al.reference_start = 100
    al.cigar = [(MATCH, 50), (INS, 10), (MATCH, 50)]
    al.seq = 'A' * 110
    
    cropped = pysam.AlignedSegment()
    cropped.reference_start = 149
    cropped.cigar = [(SOFT_CLIP, 49), (MATCH, 1), (INS, 10), (MATCH, 50)]
    cropped.seq = 'A' * 110
    
    output = sam.crop_al_to_ref_int(al, 149, 210)
    assert (output == cropped)
    
    output = sam.crop_al_to_ref_int(al, 150, 210)
    assert (output != cropped)

def test_starts_just_after_insertion():
    al = pysam.AlignedSegment()
    al.reference_start = 100
    al.cigar = [(MATCH, 50), (INS, 10), (MATCH, 50)]
    al.seq = 'A' * 110
    
    cropped = pysam.AlignedSegment()
    cropped.reference_start = 150
    cropped.cigar = [(SOFT_CLIP, 60), (MATCH, 50)]
    cropped.seq = 'A' * 110
    
    output = sam.crop_al_to_ref_int(al, 150, 210)
    assert (output == cropped)
    
    output = sam.crop_al_to_ref_int(al, 149, 210)
    assert (output != cropped)

# Merging

def test_merge_deletion():
    random.seed(0)

    bases = list('TCAG') * 100
    random.shuffle(bases)
    ref_seq = ''.join(bases)

    ref_seqs = {'random': ref_seq}
    header = pysam.AlignmentHeader.from_references(['random'], [len(ref_seq)])

    read_seq = ref_seq[100:150] + ref_seq[200:250]

    left_al = pysam.AlignedSegment(header)
    left_al.reference_name = 'random'
    left_al.reference_start = 100
    left_al.cigar = [(MATCH, 50), (SOFT_CLIP, 50)]
    left_al.seq = read_seq

    right_al = pysam.AlignedSegment(header)
    right_al.reference_name = 'random'
    right_al.reference_start = 200
    right_al.cigar = [(SOFT_CLIP, 50), (MATCH, 50)]
    right_al.seq = read_seq

    merged = sam.merge_adjacent_alignments(left_al, right_al, ref_seqs)

    expected_cigar = [(MATCH, 50), (DEL, 50), (MATCH, 50)]

    assert (merged.cigar == expected_cigar)

    # Alignments that dont point towards each other shouldn't be merged.
    left_al.reference_start = 200
    right_al.reference_start = 100

    merged = sam.merge_adjacent_alignments(left_al, right_al, ref_seqs)

    assert (merged is None)
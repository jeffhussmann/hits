from Sequencing import utilities, fastq, genomes, mapping_tools, sw, sam
from Sequencing import fasta
from Sequencing.Serialize import counts
import pysam
import string
from itertools import izip, izip_longest, chain, islice
from Circles import periodicity
import numpy as np
import os.path

def mapping_to_alignment(mapping, sam_file, base_lookup):
    ''' Convert a mapping represented by a pysam.AlignedRead into an alignment. '''
    path = []
    mismatches = set()
    deletions = set()

    for read_i, ref_i in sam.aligned_pairs_exclude_soft_clipping(mapping):
        if read_i != None:
            if ref_i == None:
                ref_i = sw.GAP
            else:
                read_base = mapping.seq[read_i]
                ref_base = base_lookup(mapping.tid, ref_i).upper()
                if read_base != ref_base:
                    mismatches.add((read_i, ref_i))

            path.append((read_i, ref_i))
        else:
            deletions.add(ref_i)

    alignment = {'path': path,
                 'XO': mapping.opt('XO'),
                 'XM': mapping.opt('XM'),
                 'is_reverse': mapping.is_reverse,
                 'mismatches': mismatches,
                 'deletions': deletions,
                 'rname': sam_file.get_reference_name(mapping.tid),
                 'query': mapping.seq,
                }
    return alignment

def produce_bowtie2_alignments(reads,
                               index_prefix,
                               genome_dir,
                               score_min,
                              ):

    bowtie2_options = {'local': True,
                       'report_up_to': 10,
                       'seed_mismatches': 1,
                       'seed_interval_function': 'C,1,0',
                       'seed_length': 10,
                      }

    sam_file, mappings = mapping_tools.map_bowtie2(index_prefix,
                                                   reads=reads,
                                                   custom_binary=True,
                                                   score_min=score_min,
                                                   yield_mappings=True,
                                                   **bowtie2_options)

    base_lookup = genomes.build_base_lookup(genome_dir, sam_file)

    mapping_groups = utilities.group_by(mappings, lambda m: m.qname)
    
    for qname, group in mapping_groups:
        group = sorted(group, key=lambda m: (m.tid, m.pos))
        alignments = [mapping_to_alignment(mapping, sam_file, base_lookup)
                      for mapping in group if not mapping.is_unmapped]
        yield qname, alignments

def get_local_alignments(read, targets):
    seq = read.seq
    seq_rc = utilities.reverse_complement(read.seq)
    all_alignments = []
    for target in targets:
        min_score = min(20, 2 * len(target.seq))
        for query, is_reverse in [(seq, False), (seq_rc, True)]:
            alignments = sw.generate_alignments(query,
                                                target.seq,
                                                'local',
                                                min_score=min_score,
                                                max_alignments=3,
                                               )
            for alignment in alignments:
                if alignment['score'] >= 0.7 * 2 * len(alignment['path']):
                    alignment['query'] = query
                    alignment['rname'] = target.name
                    alignment['is_reverse'] = is_reverse
                    all_alignments.append(alignment)

    return all_alignments

def get_edge_alignments(read, targets):
    seq = read.seq
    seq_rc = utilities.reverse_complement(read.seq)
    all_alignments = []
    min_score = 10
    for target in targets:
        for query, is_reverse in [(seq, False), (seq_rc, True)]:
            alignments = sw.generate_alignments(query,
                                                target.seq,
                                                'unpaired_adapter',
                                                min_score=min_score,
                                               )
            for alignment in alignments:
                if alignment['score'] >= 2 * len(alignment['path']):
                    alignment['query'] = query
                    alignment['rname'] = target.name
                    alignment['is_reverse'] = is_reverse
                    all_alignments.append(alignment)

    return all_alignments

def up_to_first_space(string):
    beginning = string.split(' ')[0]
    return beginning

def produce_sw_alignments(reads, genome_dirs, extra_targets, max_to_report=5):
    targets = set()

    for genome_dir in genome_dirs:
        fasta_fns = genomes.get_all_fasta_file_names(genome_dir)
        for fasta_fn in fasta_fns:
            targets.update(fasta.reads(fasta_fn))

    targets.update(extra_targets)

    for read in reads:
        alignments = get_local_alignments(read, targets) + get_edge_alignments(read, targets)
        # bowtie2 only retains up to the first space in a qname, so do the same
        # here to allow qnames to be compared
        alignments = sorted(alignments, key=lambda a: a['score'], reverse=True)
        alignments = alignments[:max_to_report]

        sanitized_name = up_to_first_space(read.name)
        yield sanitized_name, alignments

def produce_representations(alignment_groups_list):
    for alignment_group_list in izip_longest(*alignment_groups_list):
        representations = []
        qnames = set()
        for (qname, alignment_group) in alignment_group_list:
            qnames.add(qname)
            representations.extend(map(represent_alignment, alignment_group))

        if len(qnames) > 1:
            print alignment_group_list
            raise ValueError('Attempted to flatten alignment groups with different qnames:', qnames)

        qname = qnames.pop()

        yield qname, representations

def find_best_offset(R1_seq, R2_rc_seq, periodic=False):
    first_between_fractions = periodicity.compute_between_fractions(R1_seq, R2_rc_seq) 
    first_offset = np.argmax(first_between_fractions)
    first_offset_fraction = first_between_fractions[first_offset]

    second_between_fractions = periodicity.compute_between_fractions(R2_rc_seq, R1_seq) 
    second_offset = np.argmax(second_between_fractions)
    second_offset_fraction = second_between_fractions[second_offset]

    if periodic:
        within_fractions = (periodicity.compute_within_fractions(R1_seq) + periodicity.compute_within_fractions(R2_rc_seq)) / 2
        period, period_fraction = periodicity.best_period(within_fractions)
    else:
        period_fraction = 0

    # Heuristic - if both offset fractions are very high, use the one with the
    # most overlap, which is the one with the smaller offset.
    if first_offset_fraction > 0.9 and second_offset_fraction > 0.9:
        if first_offset < second_offset:
            which = 'first'
        else:
            which = 'second'
    elif first_offset_fraction >= second_offset_fraction:
        which = 'first'
    else:
        which = 'second'

    if which == 'first':
        offset = first_offset
        if period_fraction > 0.7:
            offset = offset % period
    else:
        offset = -second_offset
        if period_fraction > 0.7:
            offset = -((-offset) % period)

    best_offset_fraction = max(first_offset_fraction, second_offset_fraction)
    return offset, best_offset_fraction

def represent_alignment(alignment):
    ''' Returns text representation of the portion of the read mapped by a local
    mapping and the genomic region mapped to. Also returns a key for sorting
    such representations by how far from the left edge they begin.
    '''
    read_positions = {}
    read_to_ref = {}
    query_length = len(alignment['query'])
    for read_i, ref_i in alignment['path']:
        if read_i != sw.GAP:
            if alignment['is_reverse']:
                absolute_read_position = query_length - 1 - read_i
                match = '<'
            else:
                absolute_read_position = read_i
                match = '>'
            if ref_i != sw.GAP:
                if (read_i, ref_i) in alignment['mismatches']:
                    match = 'x'

                if ref_i - 1 in alignment['deletions']:
                    match = '\\'
                elif ref_i + 1 in alignment['deletions']:
                    match = '/'

                read_positions[absolute_read_position] = match
                read_to_ref[absolute_read_position] = ref_i
            else:
                read_positions[absolute_read_position] = '-'

    leftmost_read = min(read_positions)
    rightmost_read = max(read_positions)
    read_positions[leftmost_read] = '|'
    read_positions[rightmost_read] = '|'
    
    before_left = ' '*leftmost_read
    after_right = ' '*(query_length - rightmost_read - 1)

    width = rightmost_read - leftmost_read + 1

    ref_name_string = string.center(alignment['rname'], width)
    ref_name_line = before_left + ref_name_string + after_right
    
    # Create a line marking the reference positions of the boundaries of the
    # mapping. If there is room, summarize the mismatches and indel of the
    # mapping in between.
    if alignment['XM'] > 0 or alignment['XO'] > 0:
        mismatch_string = ' XM={0}, XO={1} '.format(alignment['XM'],
                                                    alignment['XO'],
                                                   )
    else:
        mismatch_string = ''

    mismatch_chars = list(string.center(mismatch_string, width, '-'))

    right_ref_edge_string = '{0:,}'.format(read_to_ref[rightmost_read])
    left_ref_edge_string = '{0:,}'.format(read_to_ref[leftmost_read])

    left_edge_slice = slice(1, len(left_ref_edge_string) + 1)
    right_edge_slice = slice(-(len(right_ref_edge_string) + 1), -1)

    overlaps_left = any(c != '-' for c in mismatch_chars[left_edge_slice])
    overlaps_right = any(c != '-' for c in mismatch_chars[right_edge_slice])

    if overlaps_left or overlaps_right:
        # If the verbose version overlaps, try to make a smaller string.
        mismatch_string = ' {0} {1} '.format(alignment['XM'],
                                             alignment['XO'],
                                            )
        mismatch_chars = list(string.center(mismatch_string, width, '-'))

        overlaps_left = any(c != '-' for c in mismatch_chars[left_edge_slice])
        overlaps_right = any(c != '-' for c in mismatch_chars[right_edge_slice])

        if overlaps_left or overlaps_right:
            mismatch_string = ''
    
    mismatch_chars = list(string.center(mismatch_string, width))
    extent_chars = mismatch_chars

    left_start, left_stop, _ = left_edge_slice.indices(len(extent_chars))
    right_start, right_stop, _ = right_edge_slice.indices(len(extent_chars))
    if left_stop < right_start:
        extent_chars[left_edge_slice] = list(left_ref_edge_string)
        extent_chars[right_edge_slice] = list(right_ref_edge_string)
    extent_line = before_left + ''.join(extent_chars) + after_right
    
    read_positions_chars = [' ']*query_length
    for p in read_positions:
        read_positions_chars[p] = read_positions[p]
    read_positions_line = ''.join(read_positions_chars)

    lines = (ref_name_line, extent_line, read_positions_line)
    
    return frozenset(read_positions), lines

def combine_representations(first_read_positions,
                            first_lines,
                            second_read_positions,
                            second_lines,
                           ):
    ''' Combiner to be used during the collapsing process. '''
    combined_read_positions = first_read_positions | second_read_positions
    arrays = [[' ' for p in range(len(first_lines[0]))] for l in range(len(first_lines))]
    for p in first_read_positions:
        for l, line in enumerate(first_lines):
            arrays[l][p] = line[p]
    for p in second_read_positions:
        for l, line in enumerate(second_lines):
            arrays[l][p] = line[p]
    combined_lines = [''.join(array) for array in arrays]
    return combined_read_positions, combined_lines

def leftmost_position(representation):
    read_positions, lines = representation
    return min(read_positions)

def positions_covered(representation):
    read_positions, lines = representation
    return len(read_positions)

def collapse_representations(representations):
    ''' Push nonoverlapping representation down into the same line. '''
    if not representations:
        return representations

    # Sort by leftmost position prior to collapsing
    representations = sorted(set(representations), key=leftmost_position)

    collapsed = [representations[0]]
    for read_positions, lines in representations[1:]:
        was_collapsed = False
        for t, (target_positions, target_lines) in enumerate(collapsed):
            if read_positions & target_positions == set():
                collapsed[t] = combine_representations(read_positions, lines, target_positions, target_lines)
                was_collapsed = True
                break

        if not was_collapsed:
            collapsed.append((read_positions, lines))
    
    # Sort by total positions covered after collapsing
    collapsed = sorted(collapsed, key=positions_covered)

    return collapsed

def lowercase_below_qual_threshold(seq, qual, threshold):
    ''' Returns seq with characters made lowercase at any position for which
    qual is below threshold.
    '''
    seq = list(seq)
    qual = fastq.decode_sanger(qual)
    for p, (s, q) in enumerate(zip(seq, qual)):
        if q <= threshold:
            seq[p] = s.lower()
    return ''.join(seq)

def visualize_unpaired_alignments(get_reads,
                                  sw_genome_dirs,
                                  extra_targets,
                                  bowtie2_targets,
                                  output_fn,
                                 ):

    if isinstance(get_reads, str):
        R1_fn = get_reads
        def get_reads():
            return islice(fastq.reads(R1_fn), 1000)

    R1_alignment_groups_list = []

    for genome_dir, index_prefix, score_min in bowtie2_targets:
        R1_alignment_groups = produce_bowtie2_alignments(get_reads(),
                                                         index_prefix,
                                                         genome_dir,
                                                         score_min,
                                                        )
        R1_alignment_groups_list.append(R1_alignment_groups)
        
    R1_sw_alignment_groups = produce_sw_alignments(get_reads(),
                                                   sw_genome_dirs,
                                                   extra_targets,
                                                  )
    R1_alignment_groups_list.append(R1_sw_alignment_groups)
    
    R1_representation_groups = produce_representations(R1_alignment_groups_list)

    everything = [get_reads(),
                  R1_representation_groups,
                 ]

    with open(output_fn, 'w') as output_fh:
        for R1, (R1_qname, R1_representations) in izip_longest(*everything):
            if up_to_first_space(R1.name) != R1_qname:
                raise ValueError('Iters out of sync', R1.name, R1_qname)

            R1_seq = lowercase_below_qual_threshold(R1.seq, R1.qual, 20)
            
            R1_representations = collapse_representations(R1_representations)

            output_fh.write(R1.name + '\n\n')

            for _, lines in R1_representations:
                for line in lines:
                    output_fh.write(line + '\n')

            output_fh.write(R1_seq + '\n')
            output_fh.write('\n\n')

def visualize_paired_end_mappings(get_read_pairs,
                                  sw_genome_dirs,
                                  extra_targets,
                                  bowtie2_targets,
                                  output_fn,
                                 ):

    R1_alignment_groups_list = []
    R2_alignment_groups_list = []

    if isinstance(get_read_pairs, tuple):
        R1_fn, R2_fn = get_read_pairs
        def get_R1_reads():
            return islice(fastq.reads(R1_fn), 100)
        
        def get_R2_rc_reads():
            return islice(fastq.reverse_complement_reads(R2_fn), 100)
    else:
        def get_R1_reads():
            read_pairs = islice(get_read_pairs(), 100)
            return (R1 for R1, R2 in read_pairs)
        
        def get_R2_rc_reads():
            read_pairs = islice(get_read_pairs(), 100)
            return (fastq.Read(R2.name, utilities.reverse_complement(R2.seq), R2.qual[::-1]) for R1, R2 in read_pairs)

    for genome_dir, index_prefix, score_min in bowtie2_targets:
        R1_alignment_groups = produce_bowtie2_alignments(get_R1_reads(),
                                                         index_prefix,
                                                         genome_dir,
                                                         score_min,
                                                        )
        R1_alignment_groups_list.append(R1_alignment_groups)
        
        # Design decisions made in the parsing make it easier if R2 reads are
        # reverse complemented before mapping.
        R2_alignment_groups = produce_bowtie2_alignments(get_R2_rc_reads(),
                                                         index_prefix,
                                                         genome_dir,
                                                         score_min,
                                                        )
        R2_alignment_groups_list.append(R2_alignment_groups)

    R1_sw_alignment_groups = produce_sw_alignments(get_R1_reads(),
                                                   sw_genome_dirs,
                                                   extra_targets,
                                                  )
    R1_alignment_groups_list.append(R1_sw_alignment_groups)
    
    # Design decisions made in the parsing make it easier if R2 reads are
    # reverse complemented before mapping.
    R2_sw_alignment_groups = produce_sw_alignments(get_R2_rc_reads(),
                                                   sw_genome_dirs,
                                                   extra_targets,
                                                  )
    R2_alignment_groups_list.append(R2_sw_alignment_groups)


    R1_representation_groups = produce_representations(R1_alignment_groups_list)
    R2_representation_groups = produce_representations(R2_alignment_groups_list)

    R1_reads = get_R1_reads()
    R2_rc_reads = get_R2_rc_reads()

    everything = [R1_reads,
                  R2_rc_reads,
                  R1_representation_groups,
                  R2_representation_groups,
                 ]
    
    with open(output_fn, 'w') as output_fh:
        for R1, R2_rc, (R1_qname, R1_representations), (R2_qname, R2_representations) in izip_longest(*everything):
            if up_to_first_space(R1.name) != R1_qname:
                raise ValueError('Iters out of sync', R1.name, R1_qname)
            if up_to_first_space(R2_rc.name) != R2_qname:
                raise ValueError('Iters out of sync', R2_rc.name, R2_qname)

            offset, best_offset_fraction = find_best_offset(R1.seq, R2_rc.seq)
            if best_offset_fraction > 0.7:
                gap_line = ''
                if offset >= 0:
                    R1_shift = ''
                    R2_shift = ' '*offset
                else:
                    R1_shift = ' '*(-offset)
                    R2_shift = ''
            else:
                gap_line = '\n'
                R1_shift = ''
                R2_shift = ''
            
            R1_seq = lowercase_below_qual_threshold(R1.seq, R1.qual, 20)
            R2_rc_seq = lowercase_below_qual_threshold(R2_rc.seq, R2_rc.qual, 20)
            
            R1_representations = collapse_representations(R1_representations)
            R2_representations = collapse_representations(R2_representations)

            output_fh.write(R1.name + '\n')
            output_fh.write(R2_rc.name + '\n')

            for _, lines in R1_representations:
                for line in lines:
                    output_fh.write(R1_shift + line + '\n')

            output_fh.write(R1_shift + R1_seq + '\n')
            output_fh.write(gap_line)
            output_fh.write(R2_shift + R2_rc_seq + '\n')

            for _, lines in R2_representations[::-1]:
                for line in lines[::-1]:
                    output_fh.write(R2_shift + line + '\n')

            output_fh.write('\n\n')

if __name__ == '__main__1':
    bowtie2_targets = [#('/home/jah/genomes/RF_oligos', '/home/jah/bowtie2/RF_oligos', 'C,20,0'),
                       #('/home/jah/genomes/saccharomyces_cerevisiae', '/home/jah/bowtie2/saccharomyces_cerevisiae', 'C,20,0'),
                      ]
    sw_genome_dirs = ['/home/jah/genomes/truseq',
                      '/home/jah/genomes/circle_hairpins',
                      #'/home/jah/genomes/RF_oligos',
                     ]
    extra_targets = [fasta.Read('NFBC12', 'GGCTAC'),
                     fasta.Read('R', 'TGATCTCAGATCGAAAGAAGCATGGTTGTTGTTTCTGTTAGTGTAAGCAAGCGGTTTGAAAAAGAGCGCCATGAATGACTTC'),
                    ]
    
    R1_fn = '/home/jah/projects/mutations/experiments/miseq_UT_2014_11_07/data/small_R1.fastq'
    R2_fn = '/home/jah/projects/mutations/experiments/miseq_UT_2014_11_07/data/small_R2.fastq'
    results_dir = '/home/jah/projects/mutations/experiments/miseq_UT_2014_11_07/results'
    output_fn = '/home/jah/projects/mutations/experiments/miseq_UT_2014_11_07/results/visualized_mappings.txt'
  
    with open(output_fn, 'w') as output_fh:
        visualize_local_mappings(R1_fn,
                                 R2_fn,
                                 sw_genome_dirs,
                                 extra_targets,
                                 bowtie2_targets,
                                 results_dir,
                                )

if __name__ == '__main__1':
    bowtie2_targets = [('/home/jah/genomes/pCEP4_plus_e_coli', '/home/jah/bowtie2/pCEP4_plus_e_coli', 'C,20,0'),
                      ]
    sw_genome_dirs = ['/home/jah/genomes/truseq',
                      '/home/jah/genomes/circle_hairpins',
                     ]
    extra_targets = [fasta.Read('NFBC13', 'AGTCAA'),
                    ]
    
    R1_fn = '/home/jah/projects/mutations/experiments/miseq_UT_2014_09_02/CH/data/R1_contains.fastq'
    R2_fn = '/home/jah/projects/mutations/experiments/miseq_UT_2014_09_02/CH/data/R2_contains.fastq'
    output_fn = '/home/jah/projects/mutations/experiments/miseq_UT_2014_09_02/CH/results/new_visualized_mappings.txt'
  
    with open(output_fn, 'w') as output_fh:
        visualize_local_mappings(R1_fn,
                                 R2_fn,
                                 sw_genome_dirs,
                                 extra_targets,
                                 bowtie2_targets,
                                 output_fn,
                                )

if __name__ == '__main__1':
    bowtie2_targets = [('/home/jah/genomes/saccharomyces_cerevisiae', '/home/jah/bowtie2/saccharomyces_cerevisiae', 'C,20,0'),
                      ]
    sw_genome_dirs = ['/home/jah/genomes/truseq',
                      #'/home/jah/projects/ribosomes/data/guydosh_markers/',
                     ]
    extra_targets = [fasta.Read('smRNA_linker', 'CTGTAGGCACCATCAAT'),
                    ]
    
    R1_fn = '/home/jah/projects/ribosomes/experiments/guydosh_cell/wild-type_CHX/data/SRR1042853.fastq'
    
    output_fn = '/home/jah/projects/ribosomes/experiments/guydosh_cell/wild-type_CHX/results/wild-type_CHX_structures.txt'

    def get_reads():
        return islice(fastq.reads(R1_fn), 1000)


    visualize_unpaired_alignments(get_reads,
                                  sw_genome_dirs,
                                  extra_targets,
                                  bowtie2_targets,
                                  output_fn,
                                 )

if __name__ == '__main__1':
    bowtie2_targets = [('/home/jah/genomes/pCEP4_plus_contaminants/', '/home/jah/bowtie2/pCEP4_plus_contaminants', 'C,20,0'),
                       ('/home/jah/genomes/hg19/', '/home/jah/bowtie2/hg19', 'C,50,0'),
                      ]
    sw_genome_dirs = ['/home/jah/genomes/truseq',
                     ]
    extra_targets = [fasta.Read('UTBC52', 'AACATA'),
                    ]
    
    R1_unmapped_fn = '/home/jah/scratch/projects/lnt/experiments/hiseq_UT_2015_02_03/S2/results.12.00/S2_R1_unmapped.fastq'
    R2_unmapped_fn = '/home/jah/scratch/projects/lnt/experiments/hiseq_UT_2015_02_03/S2/results.12.00/S2_R2_unmapped.fastq'
    #R1_unmapped_fn = '/home/jah/scratch/projects/lnt/experiments/hiseq_UT_2015_02_03/S2/results.12.00/test_R1.fastq'
    #R2_unmapped_fn = '/home/jah/scratch/projects/lnt/experiments/hiseq_UT_2015_02_03/S2/results.12.00/test_R2.fastq'
    #R1_unmapped_fn = '/home/jah/projects/lnt/experiments/hiseq_UT_2015_02_03/S2/data/small_R1.fastq'
    #R2_unmapped_fn = '/home/jah/projects/lnt/experiments/hiseq_UT_2015_02_03/S2/data/small_R2.fastq'
    output_fn = '/home/jah/scratch/projects/lnt/experiments/hiseq_UT_2015_02_03/S2/results.12.00/'

    visualize_paired_end_mappings(R1_unmapped_fn,
                                  R2_unmapped_fn,
                                  sw_genome_dirs,
                                  extra_targets,
                                  bowtie2_targets,
                                  output_fn,
                                 )

if __name__ == '__main__':
    bowtie2_targets = [('/home/jah/genomes/saccharomyces_cerevisiae', '/home/jah/bowtie2/saccharomyces_cerevisiae', 'C,20,0'),
                      ]
    sw_genome_dirs = ['/home/jah/genomes/truseq',
                      #'/home/jah/projects/ribosomes/data/guydosh_markers/',
                     ]
    extra_targets = [fasta.Read('smRNA_linker', 'CTGTAGGCACCATCAAT'),
                    ]
    
    R1_fn = '/home/jah/projects/ribosomes/experiments/pelechano_cell/by_chx_1_r/data/SRR1919065.fastq'
    R2_fn = '/home/jah/projects/ribosomes/experiments/pelechano_cell/by_chx_1_p/data/SRR1646644_2.fastq'
    
    output_fn = '/home/jah/projects/ribosomes/experiments/pelechano_cell/by_chx_1_p/results/test.txt'

    visualize_unpaired_alignments(R1_fn,
                                  sw_genome_dirs,
                                  extra_targets,
                                  bowtie2_targets,
                                  output_fn,
                                 )

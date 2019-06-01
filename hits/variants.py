import hashlib
import itertools
import logging

from collections import defaultdict, Counter
from functools import partial

import scipy.optimize
import scipy.stats
import numpy as np
import pysam

from . import sam, fastq, utilities, genomes, paired_end
from .utilities import base_order, base_to_index, base_to_complement_index

change_order = [('A', 'G'),
                ('T', 'C'),
                ('G', 'A'),
                ('C', 'T'),
                ('G', 'C'),
                ('C', 'G'),
                ('A', 'C'),
                ('T', 'G'),
                ('A', 'T'),
                ('T', 'A'),
                ('G', 'T'),
                ('C', 'A'),
                ('A', 'A'),
                ('T', 'T'),
                ('G', 'G'),
                ('C', 'C'),
                ('N', 'N'),
                ('-', '-'),
                ('N', 'A'),
                ('N', 'C'),
                ('N', 'T'),
                ('N', 'G'),
                ('N', '-'),
                ('A', 'N'),
                ('C', 'N'),
                ('T', 'N'),
                ('G', 'N'),
                ('-', 'N'),
                ('-', 'A'),
                ('-', 'C'),
                ('-', 'T'),
                ('-', 'G'),
                ('A', '-'),
                ('C', '-'),
                ('T', '-'),
                ('G', '-'),
               ]

coords_order = [(base_to_index[ref_char], base_to_index[read_char])
                for ref_char, read_char in change_order]

def clopper_pearson(x, n, alpha=0.05):
    if n == 0:
        return 0., 0.
    elif x == 0:
        cdf = lambda t: scipy.stats.binom.cdf(x, n, t) - alpha
        u = scipy.optimize.brentq(cdf, 0, 1)
        l = 0.
        mle = 0.
    elif x == n:
        sf = lambda t: scipy.stats.binom.sf(x - 1, n, t) - alpha
        u = 1.
        l = scipy.optimize.brentq(sf, 0, 1)
        mle = 1.
    else:
        cdf = lambda t: scipy.stats.binom.cdf(x, n, t) - (alpha / 2)
        sf = lambda t: scipy.stats.binom.sf(x - 1, n, t) - (alpha / 2)
        u = scipy.optimize.brentq(cdf, 0, 1)
        l = scipy.optimize.brentq(sf, 0, 1)
        mle = float(x) / n
    return mle - l, u - mle

def get_ineligible_positions(file_name, just_single_positions=False):
    ''' Loads a set of positions at which variants should be ignored (i.e.,
        known SNPs) from file_name. Each line of file_name should consist of
        either
        ref_seq_name
        to indiciate an entire rname that is ineligible, or
        ref_seq_name,position
        to indicate a single ineligible position, or
        ref_seq_name,start-end
        to indicate a range of ineligible positions (i.e. rDNA in yeast).
    '''
    ineligible_positions = set()
    ineligible_rnames = set()
    for line in open(file_name):
        fields = line.strip().split(',')
        if len(fields) == 1:
            ineligible_rnames.add(fields[0])
        elif len(fields) == 2:
            ref_seq_name, position = line.strip().split(',')
            range_bounds = position.split('-')
            if len(range_bounds) == 2:
                if not just_single_positions:
                    start, end = map(int, range_bounds)
                    for p in range(start, end):
                        ineligible_positions.add((ref_seq_name, p))
            else:
                ineligible_positions.add((ref_seq_name, int(position)))
    return ineligible_positions, ineligible_rnames

def hash_ineligible_positions(file_name):
    ''' Produces SHA1 hash of file_name contents. Assembled from various
        stackoverflow answers.
    '''
    sha1 = hashlib.sha1()
    with open(file_name) as fh:
        for chunk in iter(lambda: fh.read(sha1.block_size), ''):
            sha1.update(chunk)
    return sha1.hexdigest()

def count_mismatches(sam_lines,
                     max_read_length,
                     ineligible_positions,
                     ineligible_rnames=[],
                     max_qual=fastq.MAX_QUAL,
                    ):
    ''' max_read_length - a hack to make this have the same argument signature
                          as count_pe_mismatches
        ineligible_positions - set of reference positions to ignore
    '''
    circles_type_shape = (max_read_length / 2 + 1,
                          max_qual + 1,
                          len(base_order),
                          len(base_order),
                         )
    type_counts = np.zeros(shape=circles_type_shape, dtype=int)

    mismatch_ref_positions = []
    
    num_reads_skipped = 0
    total_no_indel_bases = 0
    total_no_SNP_bases = 0

    def process_alignment(sam_line):
        no_SNP_bases = 0
        # Flipping strands is kind of arbitrary for circular products.
        if sam_line['strand'] == '+':
            index_lookup = base_to_index
        else:
            index_lookup = base_to_complement_index

        alignment = sam.produce_alignment(sam_line)
        ref_seq_name = sam_line['RNAME']
        read_length = len(sam_line['SEQ'])
       
        for ref_char, read_char, qual, ref_pos, read_pos in alignment:
            # Skip known SNP positions
            if ref_seq_name in ineligible_rnames or (ref_seq_name, ref_pos) in ineligible_positions:
                continue
            no_SNP_bases += 1

            from_left = read_pos
            from_right = read_length - 1 - read_pos
            from_edge = min(from_left, from_right)
            
            ref_index = index_lookup[ref_char]
            read_index = index_lookup[read_char]
            coords = (from_edge, qual, ref_index, read_index)
            type_counts[coords] += 1
            
            if ref_char != read_char:
                #if qual >= 60:
                #    print sam_line
                #    print ref_char, read_char, qual, read_pos

                close_to_edge = from_edge < 20
                if not close_to_edge and qual > 90:
                    mismatch_ref_position = (ref_seq_name,
                                             ref_pos,
                                             ref_char,
                                             read_char,
                                             1,
                                            )
                    mismatch_ref_positions.append(mismatch_ref_position)

        return no_SNP_bases

    parsed_lines = (sam.parse_line(line) for line in sam_lines)
    
    for parsed_line in parsed_lines:
        skip = False
        if parsed_line['mapped'] == False:
            skip = True
        
        if parsed_line['MAPQ'] == 0:
            skip = True

        if sam.contains_indel(parsed_line) or sam.contains_soft_clipping(parsed_line):
            skip = True
        
        if skip:
            num_reads_skipped += 1
            continue

        total_no_indel_bases += len(parsed_line['SEQ'])
        
        no_SNP_bases = process_alignment(parsed_line)
        total_no_SNP_bases += no_SNP_bases

    mismatch_info = {'type_counts': type_counts,
                     'ref_positions': mismatch_ref_positions,
                     'num_reads_skipped': num_reads_skipped,
                     'no_indel_bases_count': total_no_indel_bases,
                     'no_SNP_bases_count': total_no_SNP_bases,
                    }
    return mismatch_info

def count_pe_mismatches(sam_file,
                        mappings,
                        max_read_length,
                        ineligible_positions,
                        ineligible_rnames=[],
                        max_qual=fastq.MAX_QUAL,
                        threshold_qual=fastq.MAX_QUAL,
                       ):
    pe_type_shape = (2 * max_read_length,
                     max_qual + 1,
                     len(base_order),
                     len(base_order),
                    )
    type_counts = np.zeros(shape=pe_type_shape, dtype=int)
    
    mismatch_ref_positions = []

    num_reads_skipped = 0
    total_no_indel_bases = 0
    total_no_SNP_bases = 0

    mismatches_per_read = Counter()
    mismatch_qnames = Counter()
    indel_qnames = Counter()
    
    def process_alignment(mapping):
        mismatches_this_read = 0
        no_SNP_bases = 0
        
        if mapping.is_reverse:
            strand = '-'
        else:
            strand = '+'

        # R1s are privileged, R2s need to be flipped.
        if mapping.is_read1:
            if strand == '+':
                index_lookup = base_to_index
            else:
                index_lookup = base_to_complement_index
        elif mapping.is_read2:
            if strand == '+':
                index_lookup = base_to_complement_index
            else:
                index_lookup = base_to_index

        alignment = sam.produce_alignment(mapping)
        rname = sam_file.getrname(mapping.tid)
        
        if rname not in ineligible_rnames:
            for ref_char, read_char, qual, ref_pos, read_pos in alignment:
                # Skip known SNP positions
                if (rname, ref_pos) in ineligible_positions:
                    continue
                no_SNP_bases += 1
                
                # Mappings may have overlapping regions truncate. We still want
                # from_left and from_right to reflect distances from the
                # untruncated edges.
                if strand == '+':
                    from_left = read_pos
                    from_right = max_read_length - 1 - read_pos
                    position = from_left
                else:
                    from_left = read_pos - (mapping.query_alignment_end - max_read_length)
                    from_right = mapping.query_alignment_end - 1 - read_pos
                    position = from_right
                
                if mapping.is_read2:
                    position = position + max_read_length
                
                ref_index = index_lookup[ref_char]
                read_index = index_lookup[read_char]
                coords = (position, qual, ref_index, read_index)
                type_counts[coords] += 1
                
                if ref_char != read_char:
                    from_edge = min(from_left, from_right)
                    close_to_edge = from_edge < 10
                    if not close_to_edge and qual >= threshold_qual:
                        mismatch_ref_position = (rname,
                                                 ref_pos,
                                                 ref_char,
                                                 read_char,
                                                 1,
                                                )
                        mismatch_ref_positions.append(mismatch_ref_position)
                        mismatch_qnames[mapping.qname] += 1
                        mismatches_this_read += 1

        return no_SNP_bases, mismatches_this_read

    mapping_pairs = itertools.izip(*[mappings]*2)

    for R1_mapping, R2_mapping in mapping_pairs:
        skip = False

        if R1_mapping.is_unmapped or R2_mapping.is_unmapped:
            skip = True
        
        if R1_mapping.tid != R2_mapping.tid:
            skip = True
        
        # Ignore, for now, reads with indels or clipping
        if sam.contains_indel_pysam(R1_mapping) or sam.contains_indel_pysam(R2_mapping):
            if sam.contains_indel_pysam(R1_mapping):
                if sam.indel_distance_from_edge(R1_mapping.cigar) > 10:
                    indel_qnames[R1_mapping.qname] += 1
            else:
                if sam.indel_distance_from_edge(R2_mapping.cigar) > 10:
                    indel_qnames[R2_mapping.qname] += 1

            skip = True
        if sam.contains_soft_clipping_pysam(R1_mapping) or sam.contains_soft_clipping_pysam(R2_mapping):
            skip = True
        
        if skip:
            num_reads_skipped += 1
            continue
       
        if R1_mapping.seq != None:
            total_no_indel_bases += len(R1_mapping.seq)
        if R2_mapping.seq != None:
            total_no_indel_bases += len(R2_mapping.seq)
        
        R1_no_SNP_bases, R1_mismatches_this_read = process_alignment(R1_mapping)
        total_no_SNP_bases += R1_no_SNP_bases
        R2_no_SNP_bases, R2_mismatches_this_read = process_alignment(R2_mapping)
        total_no_SNP_bases += R2_no_SNP_bases

        mismatches_this_pair = R1_mismatches_this_read + R2_mismatches_this_read
        mismatches_per_read[mismatches_this_pair] += 1
        #if mismatches_this_pair > 1:
        #    print R1_parsed['QNAME'], mismatches_this_pair

    mismatches_per_read = utilities.counts_to_array(mismatches_per_read)
        
    mismatch_info = {'type_counts': type_counts,
                     'ref_positions': mismatch_ref_positions,
                     'num_reads_skipped': num_reads_skipped,
                     'no_indel_bases_count': total_no_indel_bases,
                     'no_SNP_bases_count': total_no_SNP_bases,
                     'mismatches_per_read': mismatches_per_read,
                     'mismatch_qnames': mismatch_qnames,
                     'indel_qnames': indel_qnames,
                    }
    return mismatch_info

def compute_average_qualities(position_type_counts):
    position_qs = position_type_counts.sum(axis=(2, 3))
    position_average_qs = map(utilities.mean_from_histogram, position_qs)
    return position_average_qs

def compute_fractions_miscalled_as(position_type_counts, min_q=0):
    ''' Given an array of type counts whose indices are
    (position in read, quality score, reference identity, read identity),
    computes the fraction of all base calls at each position that were
    misidentified as each identity.
    '''
    num_positions, _, _, _ = position_type_counts.shape
    all_rates = np.zeros((num_positions, 6))
    for p in range(num_positions):
        position_ms = position_type_counts[p][min_q:].sum(axis=0)
        calls_of_each_type = position_ms.sum(axis=0)
        good_calls = position_ms.diagonal()
        bad_calls = calls_of_each_type - good_calls
        total_calls = calls_of_each_type.sum()
        # rates is the fraction of all calls at this position that were of each
        # type when they shouldn't have been
        rates = np.true_divide(bad_calls, max(1, total_calls))
        all_rates[p] = rates[:6]
    return all_rates

def compute_overall_rate(type_counts, min_q=None, max_q=None):
    ''' Return information about overall mismatch rates in type_counts. '''
    if min_q == None:
        min_q = len(type_counts) - 1
    if max_q == None:
        max_q = len(type_counts) - 1

    # Compute expected overall mismatch rate given quality distribution
    q_counts = type_counts[min_q:max_q + 1].sum(axis=1).sum(axis=1)
    es = np.array([10**(-q / 10.) for q in range(min_q, max_q + 1)])
    if q_counts.sum() == 0:
        expected_rate = 0
    else:
        expected_rate = np.dot(q_counts, es) / q_counts.sum()

    # Compute actual overall mismatch rate and error bars
    counts = type_counts[min_q:max_q + 1].sum(axis=0)
    total_bases = int(counts.sum())
    total_matches = counts.diagonal().sum()
    total_mismatches = total_bases - total_matches
    if total_bases != 0:
        overall_rate = float(total_mismatches) / total_bases
    else:
        overall_rate = 0. 
    overall_err = clopper_pearson(total_mismatches, total_bases)
    overall_err = np.array(overall_err).T

    rate_info = {
        'total_bases': total_bases,
        'total_mismatches': total_mismatches,
        'overall_rate': overall_rate,
        'overall_err': overall_err,
        'expected_rate': expected_rate,
    }
    return rate_info

def compute_rates(type_counts, min_q, max_q):
    ''' Computes all categories of mismatch rates from type_counts. '''
    counts = type_counts[min_q:max_q + 1].sum(axis=0)
    # counts is a 6x6 matrix of (ref_char, read_char) counts in the 
    # eligible quality range
    
    mismatch_counts = []
    rates = []
    errs = []
    for ref_coord, read_coord in coords_order[:12]:
        ref_total = sum(counts[ref_coord])
        
        if ref_total != 0:
            rate = float(counts[ref_coord, read_coord]) / ref_total
        else:
            rate = 0.
        rates.append(rate)
        
        mismatch_count = counts[ref_coord, read_coord]
        mismatch_counts.append(mismatch_count)
        
        err = clopper_pearson(counts[ref_coord, read_coord], ref_total)
        errs.append(err)

    errs = np.array(errs).T  # plt.bar expects this to be 2xN

    rate_info = {'counts': mismatch_counts,
                 'rates': rates,
                 'errs': errs,
                }
    rate_info.update(compute_overall_rate(type_counts, min_q, max_q))
    
    return rate_info

def trim_paired_end_edges(mismatch_counts, trim_length):
    num_positions, _, _, _ = mismatch_counts.shape
    read_length = num_positions / 2
    min_distance_from_edge = trim_length
    far_from_edge_slice = slice(min_distance_from_edge,
                                read_length - min_distance_from_edge,
                               )
    R1_types = mismatch_counts[:read_length]
    R2_types = mismatch_counts[read_length:]
    R1_eligible_types = R1_types[far_from_edge_slice].sum(axis=0)
    R2_eligible_types = R2_types[far_from_edge_slice].sum(axis=0)
    eligible_types = R1_eligible_types + R2_eligible_types
    return eligible_types

def amplicon_mismatches(reads, expected_seq, primers, max_qual):
    per_read_counts = Counter()
    type_shape = (len(expected_seq), max_qual + 1, 6, 6)
    type_counts = np.zeros(shape=type_shape, dtype=int)

    expected_indexes = [base_to_index[b] for b in expected_seq]

    for read in reads:
        per_read_count = 0
        if len(read.seq) == len(expected_seq):
            read_indexes = [base_to_index[b] for b in read.seq]
            read_quals = fastq.decode_sanger(read.qual)
            columns = zip(read_indexes, read_quals, expected_indexes)
            for p, (read_index, q, expected_index) in enumerate(columns):
                if read_index != expected_index and q >= max_qual - 10 and p > len(primers['forward']) and p < len(expected_seq) - len(primers['reverse']):
                    per_read_count += 1
            per_read_counts[per_read_count] += 1
            if per_read_count > 4:
                continue
            for p, (read_index, q, expected_index) in enumerate(columns):
                type_counts[p, q, expected_index, read_index] += 1

    return type_counts, per_read_counts

def reference_mismatches(bam_fn,
                         genome_dir,
                         rname,
                         start=0,
                         end=None,
                         max_qual=fastq.MAX_EXPECTED_QUAL,
                         min_distance_from_edge=10,
                         need_to_split=True,
                        ):
    genome_index = genomes.get_genome_index(genome_dir)
    bam_file = pysam.Samfile(bam_fn)
    region_fetcher = genomes.build_region_fetcher(genome_dir,
                                                  load_references=True,
                                                  sam_file=bam_file,
                                                 )

    if end == None:
        end = genome_index[rname].length - 1 

    length = end - start + 1
    
    shape = (2, # forward and reverse strands
             length,
             max_qual + 1,
             len(utilities.base_order),
             len(utilities.base_order),
            )
    position_type_counts = np.zeros(shape, int)

    total_mappings = 0
    
    sequence = region_fetcher(bam_file.gettid(rname), start, end + 1).upper()

    mappings = bam_file.fetch(reference=rname, start=start, end=end)
    if need_to_split:
        split_mappings = (paired_end.split_combined_mapping(mapping, remove_soft_clipped=True)
                          for mapping in mappings
                         )
        mappings = itertools.chain.from_iterable(split_mappings)

    for mapping in mappings:
        total_mappings += 1
        if total_mappings % 10000 == 0:
            logging.info('{0:,} mappings processed'.format(total_mappings))
        #if sam.contains_soft_clipping_pysam(mapping):
        #    continue
        if mapping.seq == None:
            # Was part of an entirely-overlapping pair
            continue

        decoded_qual = fastq.decode_sanger(mapping.qual)
        
        aligned_pairs = mapping.aligned_pairs
        for pair_i, (read_i, ref_i) in enumerate(aligned_pairs):
            if read_i == None:
                # Deletion
                read_base = '-'
                # For purposes of determining how far from the edge of the read
                # the deletion is, use the read positions on either side of it.
                before_i = pair_i
                while True:
                    before_i -= 1
                    if before_i < 0:
                        break

                    previous_read_i, _ = aligned_pairs[before_i]
                    if previous_read_i != None:
                        break

                if before_i < 0:
                    logging.info(str(mapping))
                    break

                after_i = pair_i
                while True:
                    after_i += 1
                    if after_i >= len(aligned_pairs):
                        break

                    next_read_i, _ = aligned_pairs[after_i]
                    if next_read_i != None:
                        break

                if after_i >= len(aligned_pairs):
                    logging.info(str(mapping))
                    break

                from_left = previous_read_i
                from_right = mapping.qlen - 1 - next_read_i

                read_q = max_qual # Not sure if this is the right design
            else:
                from_left = read_i
                from_right = mapping.qlen - 1 - read_i
                
                read_base = mapping.seq[read_i]
                read_q = decoded_qual[read_i]
            
            from_edge = min(from_left, from_right)
            if from_edge < min_distance_from_edge:
                continue

            if ref_i == None:
                # Insertion
                ref_base = '-'

                # Assign the insertion to the ref position of the previous 
                before_i = pair_i
                while True:
                    before_i -= 1
                    if before_i < 0:
                        break

                    _, previous_ref_i = aligned_pairs[before_i]
                    if previous_ref_i != None:
                        break

                if before_i < 0:
                    logging.info(str(mapping))
                    break

                ref_i = previous_ref_i
            else:
                ref_base = region_fetcher(mapping.tid, ref_i, ref_i + 1).upper()

            if mapping.is_reverse:
                strand = 1
            else:
                strand = 0

            if not (start <= ref_i < end):
                continue

            coords = (strand,
                      ref_i - start,
                      read_q,
                      utilities.base_to_index[ref_base],
                      utilities.base_to_index[read_base],
                     )
            position_type_counts[coords] += 1

    logging.info('{0:,} total mappings processed'.format(total_mappings))

    return position_type_counts, sequence

import pysam
import shutil
from itertools import izip
from functools import partial
from Sequencing import external_sort, sam, utilities
from . import (array_1d,
               array_2d, 
               array_3d, 
               coverage, 
               ref_positions, 
               expression, 
               log, 
               mismatches, 
               fast_mismatches, 
               counts, 
               codon_counts,
               RPKMs,
              )
from . import read_positions_hdf5 as read_positions

_concatenate_formats = {'fastq',
                        'fasta',
                        'concatenate',
                       }

_format_to_library = {'array_1d': array_1d,
                      'array_2d': array_2d,
                      'array_3d': array_3d,
                      'coverage':   coverage,
                      'read_positions': read_positions,
                      'expression': expression,
                      'log':        log,
                      'mismatches': mismatches,
                      'fast_mismatches': fast_mismatches,
                      'ref_positions':  ref_positions,
                      'counts': counts,
                      'codon_counts': codon_counts,
                      'RPKMs': RPKMs,
                     }

def _concatenate(input_file_names, output_file_name):
    with open(output_file_name, 'w') as output_file:
        for input_file_name in input_file_names:
            shutil.copyfileobj(open(input_file_name), output_file)

def _bind_map_reduce(file_format):
    library = _format_to_library[file_format]
    
    def bound_map_reduce(input_file_names, output_file_name):
        processed_inputs = (library.read_file(fn) for fn in input_file_names)
        merged = reduce(library.combine_data, processed_inputs)
        library.write_file(merged, output_file_name)

    return bound_map_reduce

def _merge_sorted_bam_files(input_file_names, merged_file_name):
    pysam.merge('-f', merged_file_name, *input_file_names)
    pysam.index(merged_file_name)

def _merge_sam_files(input_file_names, merged_file_name, are_sorted=False):
    ''' Merges a list of sam files.
        Requires all input files to have the same @SQ lines.
    '''
    sq_lines = None
    for file_name in input_file_names:
        these_sq_lines = sam.get_sq_lines(file_name)
        if sq_lines == None:
            sq_lines = these_sq_lines
        else:
            if sq_lines != these_sq_lines:
                raise ValueError('@SQ lines do not agree')

    with open(merged_file_name, 'w') as merged_file:
        for sq_line in sq_lines:
            merged_file.write(sq_line)

        input_files = [sam.open_to_reads(fn) for fn in input_file_names]
        if are_sorted:
            for line in external_sort.merge(input_files):
                merged_file.write(line)
        else:
            for input_file in input_files:
                shutil.copyfileobj(input_file, merged_file)

def _merge_interleaved_sam_files(input_file_names, merged_file_name):
    input_files = [pysam.Samfile(fn) for fn in sorted(input_file_names)]
    merged_file = pysam.Samfile(merged_file_name, 'wh', template=input_files[0])
    read_groups = [utilities.group_by(input_file, lambda r: r.qname)
                   for input_file in input_files]
    for qname, group in utilities.round_robin(read_groups):
        for line in group:
            merged_file.write(line)

def write_file(data, file_name, file_format):
    write_function = _format_to_library[file_format].write_file
    write_function(data, file_name)

def append(data, file_name, file_format):
    append_function = _format_to_library[file_format].append
    append_function(data, file_name)

def read_file(file_name, file_format):
    read_function = _format_to_library[file_format].read_file
    data = read_function(file_name)
    return data

def merge_files(input_file_names, output_file_name, file_format):
    if file_format in _concatenate_formats:
        merger = _concatenate
    elif file_format in _format_to_library:
        merger = _bind_map_reduce(file_format)
    elif file_format == 'bam':
        merger = _merge_sorted_bam_files
    elif file_format == 'sam_unsorted':
        merger = partial(_merge_sam_files, are_sorted=False)
    elif file_format == 'sam_sorted':
        merger = partial(_merge_sam_files, are_sorted=True)
    elif file_format == 'sam_interleaved':
        merger = _merge_interleaved_sam_files

    merger(input_file_names, output_file_name)

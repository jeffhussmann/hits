import os
import subprocess
from glob import glob
import Bio.SeqIO

def build_bowtie2_index(index_prefix, sequence_file_names):
    bowtie2_build_command = ['bowtie2-build',
                             ','.join(sequence_file_names),
                             index_prefix,
                            ]
    subprocess.check_call(bowtie2_build_command)

def map_bwa(file_name, genome, sai_file_name, sam_file_name, error_file_name, threads=1):
    ''' Map reads in file_name to genome with bwa. '''
    index_location = '{0}/{1}'.format(bwa_indexes_location, genome)
    with open(sai_file_name, 'w') as sai_file, open(error_file_name, 'w') as error_file:
        bwa_aln_command = ['bwa', 'aln',
                           '-t', str(threads),
                           index_location,
                           file_name,
                          ]
        subprocess.check_call(bwa_aln_command, 
                              stdout=sai_file,
                              stderr=error_file,
                             )
    
    with open(sam_file_name, 'w') as sam_file, open(error_file_name, 'a') as error_file:
        bwa_samse_command = ['bwa', 'samse',
                             '-n', '100',
                             index_location,
                             sai_file_name,
                             file_name,
                            ]
        subprocess.check_call(bwa_samse_command,
                              stdout=sam_file,
                              stderr=error_file,
                             )
    
def map_paired_bwa(R1_file_name, R2_file_name,
                   genome,
                   R1_sai_file_name, R2_sai_file_name,
                   sam_file_name,
                   error_file_name,
                   threads=1,
                  ):
    ''' Map paired end reads in R1_file_name and R2_file_name to genome with bwa. '''
    index_location = '{0}/{1}'.format(bwa_indexes_location, genome)
    with open(R1_sai_file_name, 'w') as R1_sai_file, open(error_file_name, 'w') as error_file:
        bwa_aln_R1_command = ['bwa', 'aln', '-t', str(threads), index_location, R1_file_name]
        subprocess.check_call(bwa_aln_R1_command,
                              stdout=R1_sai_file,
                              stderr=error_file,
                             )
    with open(R2_sai_file_name, 'w') as R2_sai_file, open(error_file_name, 'w') as error_file:
        bwa_aln_R2_command = ['bwa', 'aln', '-t', str(threads), index_location, R2_file_name]
        subprocess.check_call(bwa_aln_R2_command,
                              stdout=R2_sai_file,
                              stderr=error_file,
                             )

    bwa_sampe_command = ['bwa', 'sampe',
                         '-n', '100',
                         '-r', '@RG\tID:0\tSM:6',
                         '-a', '1000',
                         index_location,
                         R1_sai_file_name, R2_sai_file_name,
                         R1_file_name, R2_file_name,
                        ]
    with open(sam_file_name, 'w') as sam_file, open(error_file_name, 'a') as error_file:
        subprocess.check_call(bwa_sampe_command,
                              stdout=sam_file,
                              stderr=error_file,
                             )

def map_bowtie2(R1_file_name,
                index_prefix,
                sam_file_name,
                R2_file_name=None,
                error_file_name='/dev/null',
                **kwargs):
    ''' Map reads to index_prefix with bowtie2. '''

    kwarg_to_bowtie2_argument = [
        ('aligned_reads_file_name',   ['--al', kwargs.get('aligned_reads_file_name')]),
        ('unaligned_reads_file_name', ['--un', kwargs.get('unaligned_reads_file_name')]),
        ('aligned_pairs_file_name',   ['--al-conc', kwargs.get('aligned_pairs_file_name')]),
        ('unaligned_pairs_file_name', ['--un-conc', kwargs.get('unaligned_pairs_file_name')]),
        ('suppress_unaligned_SAM',    ['--no-unal']),
        ('memory_mapped_IO',          ['--mm']),
        ('local',                     ['--local']),
        ('ignore_quals',              ['--ignore-quals']),
        ('threads',                   ['--threads', str(kwargs.get('threads')), '--reorder']),
        ('seed_length',               ['-L', str(kwargs.get('seed_length'))]),
        ('seed_failures',             ['-D', str(kwargs.get('seed_failures'))]),
        ('reseed',                    ['-R', str(kwargs.get('reseed'))]),
        ('seed_mismatches',           ['-N', str(kwargs.get('seed_mismatches'))]),
        ('seed_interval_function',    ['-i', kwargs.get('seed_interval_function')]),
        ('gbar',                      ['--gbar', str(kwargs.get('gbar'))]),
        ('report_all',                ['-a']),
        ('report_up_to',              ['-k', str(kwargs.get('report_up_to'))]),
        ('fasta_input',               ['-f']),
        ('report_timing',             ['-t']),
        ('omit_unmapped',             ['--no-unal']),
        ('min_insert_size',           ['-I', str(kwargs.get('min_insert_size'))]),
        ('max_insert_size',           ['-X', str(kwargs.get('max_insert_size'))]),
        ('forward_forward',           ['--ff']),
        ('score_min',                 ['--score-min', kwargs.get('score_min')]),
    ]

    bowtie2_command = ['bowtie2']
    
    for kwarg, bowtie2_argument in kwarg_to_bowtie2_argument:
        if kwarg in kwargs:
            value = kwargs.pop(kwarg)
            # If set to false, don't add the argument.
            if value:
                bowtie2_command.extend(bowtie2_argument)

    bowtie2_command.extend(['-x', index_prefix])

    if R2_file_name:
        bowtie2_command.extend(['-1', R1_file_name])
        bowtie2_command.extend(['-2', R2_file_name])
    else:
        bowtie2_command.extend(['-U', R1_file_name])

    bowtie2_command.extend(['-S', sam_file_name])
    
    # kwargs are getting popped, so if anything is left, then it was something
    # that wasn't being looked for
    if len(kwargs) > 0:
        raise ValueError('Unknown keyword argument', kwargs)

    with open(error_file_name, 'w') as error_file:
        subprocess.check_call(bowtie2_command, stderr=error_file)

def map_bowtie2_paired(R1_file_name,
                       R2_file_name,
                       index_prefix,
                       sam_file_name,
                       **kwargs):
    kwargs['R2_file_name'] = R2_file_name
    map_bowtie2(R1_file_name, index_prefix, sam_file_name, **kwargs)

def map_tophat(reads_file_names,
               bowtie2_index,
               gtf_file_name,
               transcriptome_index,
               tophat_dir,
               num_threads=1,
               no_sort=False,
              ):

    joined_reads_names = ','.join(reads_file_names)

    options = [
        '--GTF', gtf_file_name,
        '--no-novel-juncs',
        '--num-threads', str(num_threads),
        '--output-dir', tophat_dir,
        '--transcriptome-index', transcriptome_index,
        '--report-secondary-alignments',
        '--read-realign-edit-dist', '0',
    ]
    if no_sort:
        options.append('--no-sort-bam')

    tophat_command = ['tophat2'] + options + [bowtie2_index, joined_reads_names]
    # tophat maintains its own logs of everything that is written to the
    # console, so discard output.
    subprocess.check_output(tophat_command, stderr=subprocess.STDOUT)

def map_tophat_paired(R1_fn,
                      R2_fn,
                      bowtie2_index,
                      gtf_file_name,
                      transcriptome_index,
                      tophat_dir,
                      num_threads=1,
                     ):
    tophat_command = ['tophat2',
                      '--GTF', gtf_file_name,
                      '--no-novel-juncs',
                      '--num-threads', str(num_threads),
                      '--output-dir', tophat_dir,
                      '--transcriptome-index', transcriptome_index,
                      bowtie2_index,
                      R1_fn,
                      R2_fn,
                     ]
    # tophat maintains its own logs of everything that is written to the
    # console, so discard output.
    subprocess.check_output(tophat_command, stderr=subprocess.STDOUT)

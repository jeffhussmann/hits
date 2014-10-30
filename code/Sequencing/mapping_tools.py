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

def map_bowtie2(index_prefix, 
                R1_file_name,
                R2_file_name,
                output_file_name,
                error_file_name='/dev/null',
                custom_binary=False,
                bam_output=False,
                unpaired_Reads=None,
                paired_Reads=None,
                **options):
    ''' Map reads to index_prefix with bowtie2. '''

    kwarg_to_bowtie2_argument = [
        ('aligned_reads_file_name',   ['--al', options.get('aligned_reads_file_name')]),
        ('unaligned_reads_file_name', ['--un', options.get('unaligned_reads_file_name')]),
        ('aligned_pairs_file_name',   ['--al-conc', options.get('aligned_pairs_file_name')]),
        ('unaligned_pairs_file_name', ['--un-conc', options.get('unaligned_pairs_file_name')]),
        ('suppress_unaligned_SAM',    ['--no-unal']),
        ('omit_secondary_seq',        ['--omit-sec-seq']),
        ('memory_mapped_IO',          ['--mm']),
        ('local',                     ['--local']),
        ('ignore_quals',              ['--ignore-quals']),
        ('threads',                   ['--threads', str(options.get('threads')), '--reorder']),
        ('seed_length',               ['-L', str(options.get('seed_length'))]),
        ('seed_failures',             ['-D', str(options.get('seed_failures'))]),
        ('reseed',                    ['-R', str(options.get('reseed'))]),
        ('seed_mismatches',           ['-N', str(options.get('seed_mismatches'))]),
        ('seed_interval_function',    ['-i', options.get('seed_interval_function')]),
        ('gbar',                      ['--gbar', str(options.get('gbar'))]),
        ('report_all',                ['-a']),
        ('report_up_to',              ['-k', str(options.get('report_up_to'))]),
        ('fasta_input',               ['-f']),
        ('report_timing',             ['-t']),
        ('omit_unmapped',             ['--no-unal']),
        ('min_insert_size',           ['-I', str(options.get('min_insert_size'))]),
        ('max_insert_size',           ['-X', str(options.get('max_insert_size'))]),
        ('forward_forward',           ['--ff']),
        ('score_min',                 ['--score-min', options.get('score_min')]),
        ('maximum_ambiguous',         ['--n-ceil', options.get('maximum_ambiguous')]),
        ('ambiguous_penalty',         ['--np', str(options.get('ambiguous_penalty'))]),
        ('allow_dovetail',            ['--dovetail']),
    ]

    if custom_binary:
        bowtie2_command = ['/home/jah/src/bowtie2-dev/bowtie2']
    else:
        bowtie2_command = ['bowtie2']

    for kwarg, bowtie2_argument in kwarg_to_bowtie2_argument:
        if kwarg in options:
            value = options.pop(kwarg)
            # If set to false, don't add the argument.
            if value:
                bowtie2_command.extend(bowtie2_argument)

    # options are getting popped, so if anything is left, then it was
    # something that wasn't being looked for
    if len(options) > 0:
        raise ValueError('Unknown keyword argument', options)

    bowtie2_command.extend(['-x', index_prefix])

    if R2_file_name != None:
        bowtie2_command.extend(['-1', R1_file_name])
        bowtie2_command.extend(['-2', R2_file_name])
    else:
        bowtie2_command.extend(['-U', R1_file_name])

    if unpaired_Reads and paired_Reads:
        raise RuntimeError('Can\'t give unpaired_Reads and paired_Reads')

    if paired_Reads and not custom_binary:
        raise RuntimeError('Can\'t used named pipes for paired Reads without custom binary because of buffer size mismatch')

    using_fifos = (unpaired_Reads != None or paired_Reads != None)

    if using_fifos:
        if os.path.exists(R1_file_name):
            raise RuntimeError('R1_file_name already exists')
        os.mkfifo(R1_file_name)

        if R2_file_name:
            if os.path.exists(R2_file_name):
                raise RuntimeError('R1_file_name already exists')
            os.mkfifo(R2_file_name)

    try:
        if not bam_output:
            bowtie2_command.extend(['-S', output_file_name])
            
            with open(error_file_name, 'w') as error_file:
                bowtie2_process = subprocess.Popen(bowtie2_command,
                                                   stderr=error_file,
                                                  )
        else:
            with open(error_file_name, 'w') as error_file:
                bowtie2_process = subprocess.Popen(bowtie2_command,
                                                   stdout=subprocess.PIPE,
                                                   stderr=error_file,
                                                  )
                bam_command = ['samtools',
                               'view',
                               '-bu',
                               '-o',
                               output_file_name,
                               '-',
                              ]
                samtools_process = subprocess.Popen(bam_command,
                                                    stdin=bowtie2_process.stdout,
                                                   )
                bowtie2_process.stdout.close()

        if unpaired_Reads:
            with open(R1_file_name, 'w') as R1_fh:
                for R1 in unpaired_Reads:
                    R1_fh.write(str(R1))
        elif paired_Reads:
            with open(R1_file_name, 'w') as R1_fh, open(R2_file_name, 'w') as R2_fh:
                for R1, R2 in paired_Reads:
                    R1_fh.write(str(R1))
                    R2_fh.write(str(R2))

        if not bam_output:
            bowtie2_process.wait()
            if bowtie2_process.returncode != 0:
                raise subprocess.CalledProcessError(bowtie2_command,
                                                    bowtie2_process.returncode,
                                                   )
        else:
            samtools_process.communicate()
            if samtools_process.returncode != 0:
                raise subprocess.CalledProcessError(bowtie2_command,
                                                    samtools.returncode,
                                                   )
    finally:
        if using_fifos:
            os.remove(R1_file_name)
            if R2_file_name:
                os.remove(R2_file_name)

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

import os
import tempfile
import subprocess
import threading
import Queue
from Sequencing import fastq
from Sequencing import sam

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

class ThreadWriter(threading.Thread):
    def __init__(self, fifo_name, reads):
        threading.Thread.__init__(self)
        self.daemon = True
        self.fifo_name = fifo_name
        self.reads = reads
        self.start()

    def run(self):
        with open(self.fifo_name, 'w') as fifo_fh:
            for read in self.reads:
                fifo_fh.write(str(read))

class PairedThreadWriter(threading.Thread):
    def __init__(self, R1_fifo_name, R2_fifo_name, read_pairs):
        threading.Thread.__init__(self)
        self.daemon = True
        self.R1_fifo_name = R1_fifo_name
        self.R2_fifo_name = R2_fifo_name
        self.read_pairs = read_pairs
        self.start()

    def run(self):
        with open(self.R1_fifo_name, 'w') as R1_fifo_fh, open(self.R2_fifo_name, 'w') as R2_fifo_fh:
            for R1, R2 in self.read_pairs:
                R1_fifo_fh.write(str(R1))
                R2_fifo_fh.write(str(R2))

class ThreadReader(threading.Thread):
    def __init__(self, fifo_name):
        threading.Thread.__init__(self)
        self.daemon = True
        self.fifo_name = fifo_name
        self.queue = Queue.Queue(maxsize=1000)
        self.start()
        
    def run(self):
        with open(self.fifo_name, 'r') as fifo_fh:
            for read in fastq.reads(fifo_fh):
                self.queue.put(read)

        self.queue.put(None)

    def output(self):
        for read in iter(self.queue.get, None):
            yield read

class PairedThreadReader(threading.Thread):
    def __init__(self, R1_fifo_name, R2_fifo_name):
        threading.Thread.__init__(self)
        self.daemon = True
        self.R1_fifo_name = R1_fifo_name
        self.R2_fifo_name = R2_fifo_name
        self.queue = Queue.Queue(maxsize=1000)
        self.start()
        
    def run(self):
        with open(self.R1_fifo_name) as R1_fifo_fh, open(self.R2_fifo_name) as R2_fifo_fh:
            for read_pair in fastq.read_pairs(R1_fifo_fh, R2_fifo_fh):
                self.queue.put(read_pair)

        self.queue.put(None)

    def output(self):
        for read_pair in iter(self.queue.get, None):
            yield read_pair

def prep_fifos(is_paired):
    temp_dir = tempfile.mkdtemp()
    template_fn = '{0}/R%_fifo'.format(temp_dir)
    R1_fn = '{0}/R1_fifo'.format(temp_dir)
    os.mkfifo(R1_fn)

    if not is_paired:
        R2_fn = None
    else:
        R2_fn = '{0}/R2_fifo'.format(temp_dir)
        os.mkfifo(R2_fn)

    return temp_dir, template_fn, R1_fn, R2_fn

def cleanup_fifos(temp_dir, R1_fn, R2_fn):
    os.remove(R1_fn)
    if R2_fn != None:
        os.remove(R2_fn)
    os.rmdir(temp_dir)

def prep_bowtie2_command(index_prefix,
                         R1_fn,
                         R2_fn,
                         custom_binary,
                         **options):
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

    if R2_fn != None:
        bowtie2_command.extend(['-1', R1_fn])
        bowtie2_command.extend(['-2', R2_fn])
    else:
        bowtie2_command.extend(['-U', R1_fn])

    return bowtie2_command

def launch_bowtie2_pipeline(bowtie2_command,
                            output_file_name,
                            error_file_name,
                            bam_output,
                           ):
    error_file = open(error_file_name, 'w')
    if not bam_output:
        bowtie2_command.extend(['-S', output_file_name])
        
        bowtie2_process = subprocess.Popen(bowtie2_command,
                                           stderr=error_file,
                                          )
        last_process = bowtie2_process
    else:
        bowtie2_process = subprocess.Popen(bowtie2_command,
                                           stdout=subprocess.PIPE,
                                           stderr=error_file,
                                          )
        view_command = ['samtools', 'view', '-ubh', '-']
        sort_command = ['samtools', 'sort',
                        '-T', output_file_name,
                        '-o', output_file_name,
                        '-',
                       ]
        view_process = subprocess.Popen(view_command,
                                        stdin=bowtie2_process.stdout,
                                        stdout=subprocess.PIPE,
                                       )
        sort_process = subprocess.Popen(sort_command,
                                        stdin=view_process.stdout,
                                       )
        bowtie2_process.stdout.close()
        view_process.stdout.close()
        last_process = sort_process

    return last_process

def map_bowtie2(index_prefix,
                R1_fn,
                R2_fn,
                output_file_name,
                error_file_name='/dev/null',
                custom_binary=False,
                bam_output=False,
                unpaired_Reads=None,
                paired_Reads=None,
                **options):
    if unpaired_Reads and paired_Reads:
        raise RuntimeError('Can\'t give unpaired_Reads and paired_Reads')

    if paired_Reads and not custom_binary:
        raise RuntimeError('Can\'t used named pipes for paired Reads without custom binary because of buffer size mismatch')

    using_input_fifos = (unpaired_Reads != None or paired_Reads != None)
    is_paired = (R2_fn != None or paired_Reads != None)

    try:
        if using_input_fifos:
            temp_input_dir, _, R1_fn, R2_fn = prep_fifos(is_paired)
        
        bowtie2_command = prep_bowtie2_command(index_prefix,
                                               R1_fn,
                                               R2_fn,
                                               custom_binary,
                                               **options)

        last_process = launch_bowtie2_pipeline(bowtie2_command,
                                               output_file_name,
                                               error_file_name,
                                               bam_output,
                                              )

        if using_input_fifos:
            if not is_paired:
                writer = ThreadWriter(R1_fn, unpaired_Reads)
            else:
                writer = PairedThreadWriter(R1_fn, R2_fn, paired_Reads)

        last_process.wait()
        if last_process.returncode != 0:
            raise subprocess.CalledProcessError(bowtie2_command,
                                                last_process.returncode,
                                               )
        if bam_output:
            sam.index_bam(output_file_name)
    finally:
        if using_input_fifos:
            cleanup_fifos(temp_input_dir, R1_fn, R2_fn)

def map_bowtie2_yield_unaligned(index_prefix,
                                R1_fn,
                                R2_fn,
                                output_file_name,
                                error_file_name='/dev/null',
                                custom_binary=False,
                                bam_output=False,
                                unpaired_Reads=None,
                                paired_Reads=None,
                                **options):
    if unpaired_Reads and paired_Reads:
        raise RuntimeError('Can\'t give unpaired_Reads and paired_Reads')

    if paired_Reads and not custom_binary:
        raise RuntimeError('Can\'t used named pipes for paired Reads without custom binary because of buffer size mismatch')

    using_input_fifos = (unpaired_Reads != None or paired_Reads != None)
    is_paired = (R2_fn != None or paired_Reads != None)

    try:
        if using_input_fifos:
            temp_input_dir, _, R1_fn, R2_fn = prep_fifos(is_paired)
        
        temp_output_dir, unal_template_fn, unal_R1_fn, unal_R2_fn = prep_fifos(is_paired)
        if not is_paired:
            options['unaligned_reads_file_name'] = unal_R1_fn
        else:
            options['unaligned_pairs_file_name'] = unal_template_fn

        bowtie2_command = prep_bowtie2_command(index_prefix,
                                               R1_fn,
                                               R2_fn,
                                               custom_binary,
                                               **options)

        last_process = launch_bowtie2_pipeline(bowtie2_command,
                                               output_file_name,
                                               error_file_name,
                                               bam_output,
                                              )

        if using_input_fifos:
            if not is_paired:
                writer = ThreadWriter(R1_fn, unpaired_Reads)
            else:
                writer = PairedThreadWriter(R1_fn, R2_fn, paired_Reads)

        if not is_paired:
            reader = ThreadReader(unal_R1_fn)
            for read in reader.output():
                yield read
        else:
            reader = PairedThreadReader(unal_R1_fn,
                                        unal_R2_fn,
                                       )
            for read_pair in reader.output():
                yield read_pair
        
        last_process.wait()
        if last_process.returncode != 0:
            raise subprocess.CalledProcessError(bowtie2_command,
                                                last_process.returncode,
                                               )
        if bam_output:
            sam.index_bam(output_file_name)
    finally:
        if using_input_fifos:
            cleanup_fifos(temp_input_dir, R1_fn, R2_fn)

        cleanup_fifos(temp_output_dir, unal_R1_fn, unal_R2_fn)

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

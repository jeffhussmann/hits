import argparse
import subprocess
import time
import os
import shutil
import re
import logging
from . import split_file
from . import launcher
from Sequencing import Serialize

def extend_stages(whole_stages, specific_stages):
    for whole_stage, specific_stage in zip(whole_stages, specific_stages):
        whole_stage.extend(specific_stage)

class MapReduceExperiment(object):
    specific_results_files = [
        ('summary', Serialize.log, '{name}_summary.txt'),
        ('log', None, 'log.txt'),
    ]
    specific_figure_files = []
    specific_outputs = []
    specific_work = []
    specific_cleanup = []

    def __init__(self, **kwargs):
        self.name = kwargs['name']
        self.work_prefix = kwargs['work_prefix'].rstrip('/')
        self.scratch_prefix = kwargs['scratch_prefix'].rstrip('/')
        self.relative_results_dir = kwargs['relative_results_dir'].strip('/')
        self.num_pieces = kwargs['num_pieces']
        self.which_piece = kwargs['which_piece']
        
        suffix = split_file.generate_suffix(self.num_pieces, self.which_piece)

        self.scratch_results_dir = '{0}/{1}{2}'.format(self.scratch_prefix,
                                                       self.relative_results_dir,
                                                       suffix,
                                                      )
        self.work_results_dir = '{0}/{1}'.format(self.work_prefix,
                                                 self.relative_results_dir,
                                                )
        
        if self.which_piece == -1:
            # Sentinel value that indicates this is the merged experiment.
            # Only create the directory in this instance to avoid race
            # conditions.
            if not os.path.isdir(self.work_results_dir):
                os.makedirs(self.work_results_dir)
        else:
            # Only need to create the scratch directory if this ISN't the merged
            # experiment.
            if not os.path.isdir(self.scratch_results_dir):
                os.makedirs(self.scratch_results_dir)

        # Build file name templates and lists of outputs, work, and cleanup from
        # the specific lists of these things specified by each member in the
        # inheritance chain.
        class_hierarchy = [c for c in self.__class__.__mro__[::-1] if c != object]
        self.results_files = []
        self.figure_files = []
        self.outputs = [[] for _ in range(self.num_stages)]
        self.work = [[] for _ in range(self.num_stages)]
        self.cleanup = [[] for _ in range(self.num_stages)]
        for c in class_hierarchy:
            self.results_files.extend(c.specific_results_files)
            self.figure_files.extend(c.specific_figure_files)
            extend_stages(self.outputs, c.specific_outputs)
            extend_stages(self.work, c.specific_work)
            extend_stages(self.cleanup, c.specific_cleanup)
    
        for stage in range(self.num_stages):
            for kind in ['summary_stage', 'timing']:
                key = '{0}_{1}'.format(kind, stage)
                tail_template = '{{name}}_{0}_{1}.txt'.format(kind, stage)
                self.results_files.append((key, Serialize.log, tail_template))
                self.outputs[stage].append(key)

        self.cleanup[-1].append('consolidate_summaries')

        self.make_file_names()
        self.summary = []

        if self.which_piece == -1:
            piece_string = ' '
        else:
            piece_string = ' {0:{length}} / {1} '.format(self.which_piece, self.num_pieces, length=len(str(self.num_pieces)))

        format_string = '%(asctime)s{0}%(message)s'.format(piece_string)

        logging.basicConfig(filename=self.file_names['log'],
                            level=logging.INFO,
                            format=format_string, 
                           )

    @classmethod
    def from_description_file_name(cls, description_file_name, num_pieces=1, which_piece=-1):
        description = parse_description(description_file_name)
        return cls(num_pieces=num_pieces, which_piece=which_piece, **description)

    def make_file_names(self):
        self.file_names = {}
        self.merged_file_names = {}
        self.file_types = {}

        for key, serialize_type, tail_template in self.results_files:
            file_tail = tail_template.format(name=self.name)
            self.file_names[key] = '{0}/{1}'.format(self.scratch_results_dir, file_tail)
            self.merged_file_names[key] = '{0}/{1}'.format(self.work_results_dir, file_tail)
            self.file_types[key] = serialize_type

        self.figure_file_names = {}
        for key, tail_template in self.figure_files:
            file_tail = tail_template.format(name=self.name)
            self.figure_file_names[key] = '{0}/{1}'.format(self.work_results_dir, file_tail)

        if self.which_piece == -1:
            self.file_names = self.merged_file_names

    def consolidate_summaries(self):
        stage_file_names = [self.file_names['summary_stage_{0}'.format(stage)]
                            for stage in range(self.num_stages)]
        Serialize.log.consolidate_stages(stage_file_names, self.file_names['summary'])

    def write_file(self, key, data):
        file_format = self.file_types[key]
        file_format.write_file(data, self.file_names[key])

    def read_file(self, key, merged=False, **kwargs):
        if merged == False:
            file_name = self.file_names[key]
        else:
            file_name = self.merged_file_names[key]
        
        file_format = self.file_types[key]
        data = file_format.read_file(file_name, **kwargs)

        return data

    def do_work(self, stage):
        logging.info('Beginning work for stage {0}'.format(stage))

        times = []
        for function_name in self.work[stage]:
            logging.info('Starting function {0}'.format(function_name))
            start_time = time.time()
            self.__getattribute__(function_name)()
            end_time = time.time()
            times.append((function_name, end_time - start_time))

        self.write_file('timing_{0}'.format(stage), times)
        self.write_file('summary_stage_{0}'.format(stage), self.summary)

    def do_cleanup(self, stage):
        times = []
        for function_name in self.cleanup[stage]:
            start_time = time.time()
            self.__getattribute__(function_name)()
            end_time = time.time()
            times.append((function_name, end_time - start_time))
        
        Serialize.log.append(times, self.merged_file_names['timing_{0}'.format(stage)])

def controller(ExperimentClass, script_path):
    args = parse_arguments()
    if args.subparser_name == 'launch':
        launch(args, script_path, ExperimentClass.num_stages)
    elif args.subparser_name == 'process':
        process(args, ExperimentClass)
    elif args.subparser_name == 'finish':
        finish(args, ExperimentClass)

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--job_dir',
                        required=True,
                       )

    subparsers = parser.add_subparsers(dest='subparser_name')

    parser_process = subparsers.add_parser('process')
    parser_process.add_argument('--num_pieces', type=int,
                                help='how many pieces',
                                default=1,
                               )
    parser_process.add_argument('--which_piece', type=int,
                                help='which piece this is',
                                default=0,
                               )
    parser_process.add_argument('--stage', type=int,
                                help='stage',
                                default=0,
                               )
    
    parser_launch = subparsers.add_parser('launch')
    parser_launch.add_argument('--num_pieces', type=int,
                               default=1,
                              )

    parser_finish = subparsers.add_parser('finish')
    parser_finish.add_argument('--num_pieces', type=int,
                                help='how many pieces',
                                default=1,
                               )
    parser_finish.add_argument('--stage', type=int,
                                help='stage',
                                default=0,
                               )

    args = parser.parse_args()
    return args

def parse_description(description_fn):
    description = dict(line.strip().split() for line in open(description_fn)
                       if not line.startswith('#'))
    return description

def launch(args, script_path, num_stages):
    description_file_name = '{0}/description.txt'.format(args.job_dir)
    description = parse_description(description_file_name)

    def make_process_command(args, which_piece, stage):
        command = ['python', script_path,
                   '--job_dir', args.job_dir,
                   'process',
                   '--num_pieces', str(args.num_pieces),
                   '--which_piece', str(which_piece),
                   '--stage', str(stage),
                  ]
        command_string = ' '.join(command) + '\n'
        return command_string
    
    def make_finish_command(args, stage):
        command = ['python', script_path,
                   '--job_dir', args.job_dir,
                   'finish',
                   '--num_pieces', str(args.num_pieces),
                   '--stage', str(stage),
                  ]
        command_string = ' '.join(command) + '\n'
        return command_string

    process_file_names = []
    finish_file_names = []
    job_names = []

    for stage in range(num_stages):
        job_name = '{0}_{1}_{2}'.format(description['name'],
                                        args.num_pieces,
                                        stage,
                                       )
        job_names.append(job_name)

        process_file_name = '{0}/process_{1}_stage_{2}'.format(args.job_dir,
                                                               args.num_pieces,
                                                               stage,
                                                              )
        process_file_names.append(process_file_name)

        finish_file_name = '{0}/finish_{1}_stage_{2}'.format(args.job_dir,
                                                             args.num_pieces,
                                                             stage,
                                                            )
        finish_file_names.append(finish_file_name)
                              
        with open(process_file_name, 'w') as process_file:
            for which_piece in range(args.num_pieces):
                line = make_process_command(args, which_piece, stage)
                process_file.write(line)

        with open(finish_file_name, 'w') as finish_file:
            line = make_finish_command(args, stage)
            finish_file.write(line)

    # Launch with qsub on lonestar, sbatch on stampede, parallel otherwise.
    hostname = os.environ.get('HOSTNAME', '')
    if 'tacc' in hostname:
        if 'ls4' in hostname:
            submitter = 'qsub'
            def get_job_id(output):
                return re.search(r'Your job (\d+)', output).group(1)
        elif 'stampede' in hostname:
            submitter = 'sbatch'
            def get_job_id(output):
                return re.search(r'Submitted batch job (\d+)', output).group(1)

        starting_path = os.getcwd()
        os.chdir(args.job_dir)
        launcher_file_name = launcher.create(
            job_names[0],
            process_file_names[0],
            time='00:20:00',
            optional_finish='bash {0}'.format(finish_file_names[0]),
        )
        output = subprocess.check_output([submitter, launcher_file_name])
        this_job_id = get_job_id(output)
        print 'Launched stage {0} with jid {1}'.format(stage, this_job_id)

        for stage in range(1, num_stages):
            previous_job_id = this_job_id
            launcher_file_name = launcher.create(
                job_names[stage],
                process_file_names[stage],
                time='00:20:00',
                optional_finish='bash {0}'.format(finish_file_names[stage]),
                hold_jid=previous_job_id,
            )
            output = subprocess.check_output([submitter, launcher_file_name])
            this_job_id = get_job_id(output)
            print 'Launched stage {0} with jid {1}, holding on {2}'.format(stage,
                                                                           this_job_id,
                                                                           previous_job_id,
                                                                          )
        os.chdir(starting_path)
    else:
        for stage in range(num_stages):
            print 'Launched stage {0} with parallel'.format(stage)
            subprocess.check_call('parallel < {0}'.format(process_file_names[stage]), shell=True)
            subprocess.check_call('bash {0}'.format(finish_file_names[stage]), shell=True)

def process(args, ExperimentClass):
    description_file_name = '{0}/description.txt'.format(args.job_dir)
    description = parse_description(description_file_name)

    experiment = ExperimentClass(num_pieces=args.num_pieces,
                                 which_piece=args.which_piece,
                                 **description)
    experiment.do_work(args.stage)

def finish(args, ExperimentClass):
    description_file_name = '{0}/description.txt'.format(args.job_dir)
    description = parse_description(description_file_name)
    
    merged = ExperimentClass(num_pieces=args.num_pieces,
                             which_piece=-1,
                             **description)
    pieces = [ExperimentClass(num_pieces=args.num_pieces,
                              which_piece=which_piece,
                              **description)
              for which_piece in range(args.num_pieces)]

    merge_times = []
    for key in merged.outputs[args.stage]:
        piece_file_names = [piece.file_names[key] for piece in pieces]
        merged_file_name = merged.merged_file_names[key]

        start_time = time.time()
        Serialize.merge_files(piece_file_names,
                              merged_file_name,
                              merged.file_types[key],
                             )
        end_time = time.time()
        merge_times.append(('Merging {}'.format(key), end_time - start_time))

    Serialize.log.append(merge_times, merged.merged_file_names['timing_{0}'.format(args.stage)])
    merged.do_cleanup(args.stage)

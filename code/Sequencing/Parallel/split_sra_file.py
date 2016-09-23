import subprocess32 as subprocess
import os
import contextlib
import Sequencing.Parallel
import download_GSE
import sys

def piece(srr_fn, num_pieces, which_piece, paired=False):
    root, _ = os.path.splitext(srr_fn)
    xml_fn = '{0}.xml'.format(root)
    if not os.path.exists(xml_fn):
        raise ValueError('XML doesn\'t exist for {0}'.format(srr_fn))

    info = download_GSE.parse_run_xml(xml_fn)
    total_spots = info['total_spots']

    bounds = Sequencing.Parallel.get_bounds(total_spots, num_pieces)
    first = bounds[which_piece] + 1
    last = bounds[which_piece + 1]

    with dump_spots(srr_fn, first, last, paired) as lines:
        for line in lines:
            yield line

@contextlib.contextmanager
def dump_spots(srr_fn, first, last, paired):
    if paired:
        name_format = '@$ac.$si.$ri'
    else:
        name_format = '@$ac.$si'

    command = ['fastq-dump',
               '--dumpbase',
               '--minSpotId', str(first),
               '-maxSpotId', str(last),
               '--defline-seq', name_format,
               '--stdout',
               srr_fn,
              ]

    if paired:
        command.insert(1, '--split-spot')

    process = subprocess.Popen(command,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                              )
    lines = iter(process.stdout)
    try:
        yield lines
    finally:
        process.terminate()
        process.stdout.close()
        for line in process.stderr:
            if not (line.startswith('Read') or line.startswith('Written')):
                sys.stderr.write(line)

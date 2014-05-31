import numpy as np
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

include_dirs = [np.get_include()]
ext_modules = [#Extension('adapters_cython', ['Circles/adapters_cython.pyx'], include_dirs=include_dirs),
               #Extension('amplicons_cython', ['Circles/amplicons_cython.pyx'], include_dirs=include_dirs),
               #Extension('barcodes_cython', ['Circles/barcodes_cython.pyx'], include_dirs=include_dirs),
               #Extension('periodicity_cython', ['Circles/periodicity_cython.pyx'], include_dirs=include_dirs),
               #Extension('bayes_base_call_cython', ['Circles/bayes_base_call_cython.pyx'], include_dirs=include_dirs),
               Extension('fastq_cython', ['Sequencing/fastq_cython.pyx'], include_dirs=include_dirs),
               #Extension('seed_extend_cython', ['Circles/seed_extend_cython.pyx'], include_dirs=include_dirs),
              ]

setup(
    name='Sequencing',
    version='0.1',
    author='Jeff Hussmann',
    author_email='jeff.hussmann@gmail.com',
    ext_package='Sequencing',
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
)

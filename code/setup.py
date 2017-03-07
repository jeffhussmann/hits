import numpy as np
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

include_dirs = [np.get_include()]

ext_modules = [Extension('adapters_cython', ['Sequencing/adapters_cython.pyx'], include_dirs=include_dirs),
               Extension('fastq_cython', ['Sequencing/fastq_cython.pyx'], include_dirs=include_dirs),
               Extension('sw_cython', ['Sequencing/sw_cython.pyx'], include_dirs=include_dirs),
              ]

setup(
    name='Sequencing',
    version='0.1',
    author='Jeff Hussmann',
    author_email='jeff.hussmann@gmail.com',
    ext_package='Sequencing',
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
    packages=['Sequencing',
              'Sequencing/Visualize',
              'Sequencing/Visualize/interactive',
             ],
    package_data={'Sequencing/Visualize/interactive': ['*.coffee',
                                                       'example_df.txt',
                                                       '*.json',
                                                      ],
                 },
)

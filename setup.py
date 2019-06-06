from setuptools import setup, Extension

setup(
    name='hits',
    version='0.0.4',

    author='Jeff Hussmann',
    author_email='jeff.hussmann@gmail.com',
    url='https://github.com/jeffhussmann/hits',
    description='utilities for processing high-throughput sequencing experiments',

    setup_requires=['cython'],

    ext_package='hits',
    ext_modules=[
        Extension('adapters_cython', ['hits/adapters_cython.pyx']),
        Extension('fastq_cython', ['hits/fastq_cython.pyx']),
        Extension('sw_cython', ['hits/sw_cython.pyx']),
    ],

    packages=[
        'hits',
        'hits/visualize',
        'hits/visualize/interactive',
    ],
    package_data={
        'hits/visualize/interactive': [
            '*.coffee',
            'example_df.txt',
            '*.json',
        ],
    },

    python_requires='>=3.6',

    install_requires=[
        'biopython',
        'matplotlib',
        'numpy',
        'ipython',
        'ipywidgets',
        'pandas',
        'pillow',
        'pysam',
        'PyYAML',
        'scipy',
    ],
)

from setuptools import setup, Extension

setup(
    name='hits',
    version='0.3.0',

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
            '*.js',
            'example_df.txt',
            '*.json',
        ],
    },

    python_requires='>=3.6',

    install_requires=[
        'biopython>=1.72',
        'bokeh>=1.0.4',
        'ipython>=7.8.0',
        'ipywidgets>=7.4.2',
        'matplotlib>=3.0.2',
        'numpy>=1.15.4',
        'pandas>=0.23.4',
        'pillow>=5.3.0',
        'pysam>=0.15.1',
        'pyyaml>=3.13',
        'scipy>=1.2.1',
        'statsmodels==0.12.1',
    ],
)

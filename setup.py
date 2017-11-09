from setuptools import setup, find_packages
from codecs import open
from os import path
import os
import re
import io

def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pyseer',
    version=find_version("pyseer"),
    description='Sequence Elements Enrichment Analysis (SEER), python implementation',
    long_description=long_description,
    url='https://github.com/mgalardini/pyseer',
    author='Marco Galardini and John Lees',
    author_email='marco@ebi.ac.uk',
    license='GPL3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='gwas bacteria k-mer',
    # TODO: use entry_points instead
    scripts=['pyseer', 'square_mash'],
    #
    install_requires=['numpy',
                      'scipy',
                      'pandas',
                      'statsmodels'],
)

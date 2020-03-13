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
    version=find_version("pyseer/__init__.py"),
    description='Sequence Elements Enrichment Analysis (SEER), python implementation',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/mgalardini/pyseer',
    author='Marco Galardini and John Lees',
    author_email='marco@ebi.ac.uk',
    license='Apache Software License',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='gwas bacteria k-mer',
    packages=['pyseer',
              'pyseer.fastlmm',
              'pyseer.kmer_mapping'],
    entry_points={
        "console_scripts": [
            'pyseer = pyseer.__main__:main',
            'square_mash = pyseer.mash:main',
            'similarity_pyseer = pyseer.similarity:main',
            'scree_plot_pyseer = pyseer.scree_plot:main',
            'phandango_mapper = pyseer.kmer_mapping.phandango_plot:main',
            'annotate_hits_pyseer = pyseer.kmer_mapping.annotate_hits:main',
            'enet_predict_pyseer = pyseer.enet_predict:main'
            ]
    },
    install_requires=['numpy',
                      'scipy',
                      'pandas',
                      'statsmodels>=0.10.0',
                      'scikit-learn',
                      'pysam',
                      'DendroPy',
                      'matplotlib',
                      'pybedtools',
                      'tqdm',
                      'glmnet_python@https://github.com/johnlees/glmnet_python/archive/v1.0.2.zip'],
    dependency_links = ['https://github.com/johnlees/glmnet_python/tarball/v1.0.2#egg=glmnet_python-v1.0.2'],
    test_suite="tests",
)

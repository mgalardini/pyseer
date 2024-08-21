# CHANGELOG
v1.3.12 (August 2024)
- Improvement: `annotate_hits.py` has a new argument to specify more GFF types to use
- Improvement: `count_patterns.py` has a new argument to only print the p-value threshold (thanks to Bamu Damaris)
- Bugfix: refuse to run if certain arguments are not provided
- Docs: fix readthedocs' configuration
- CI: speeding things up using mamba
- CI: workaround to avoid GH actions crashing due to glmnet_py
- Tests: avoid failures due to newer versions of sklearn

v1.3.11 (May 2023)
- Improvement: Rtab files can use a empty string to signify missing values
- Bugfix: avoid crashes wn running wgGWAS with covariates
- Bugfix (summarise_annotations.py): more robust parsing
- Bugfix: avoid crashes when distance matrices are not squared
- Bugfix: properly format output when using `--no-distances`

v1.3.10 (July 2022)
- Docs: more clarity on effect sizes for LMM
- Bugfix (similarity): make sure the script actually uses all variants
- Bugfix (QQ-plot): avoid axis inversion with newer versions of statsmodels (thanks to Julian Libiseller-Egger)
- Bugfix: use kwargs when invoking functions from scikit-learn
- Bugfix (annotate hits): avoid crashes when full path to reference genomes contains non-ASCII chars
- Bugfix: WG models should not be ran with the `--output-patterns` function
- Bugfix: avoid a crash when saving a WG model with covariates
- Bugfix: add the "lineage" header when running a whole genome model

v1.3.9 (June 2021)
- Bugfix: avoid a crash when providing lineages in whole genome mode

v1.3.8 (May 2021)
- Improvement: fall back to Firth regression when encountering a matrix inversion error (thanks to Julian Libiseller-Egger)
- Bugfix: check for zero passing variants in read_all (enet, thanks to Julian Libiseller-Egger)
- Bugfix: use len(all_strains) instead of len(sample_order) to determine shape of sparse matrix in load_all_vars (enet, thanks to Julian Libiseller-Egger) 
- Bugfix: properly report all filtered variants
- Bugfix: --lmm requires either --similarity or a LMM cache

v1.3.7 (May 2021)
- Enhancement: check that provided phenotypes are of numeric type
- Bugfix: properly report all filtered variants
- Bugfix: don't crash if regression fails for missing data
- Bugfix for whole genome regression and lineages
- Bugfix/docs: properly report that the covariates file should have a header
- scripts/summarise_annotations.py can work with unadjusted p-values (new option, thanks to Lindsay Clark)
- scripts/phylogeny_distance.py: no need to reroot twice with updated dendropy
- CI: transition to GitHub actions

v1.3.6 (May 2020)
- Bugfix for missing variants in VCF files (now properly handled)
- Bugfixes for k-mer mapping (lack of annotation, bwa fastmap with many hits
- Explicitly look for duplicated samples in input files
- Removed package from pypi (glmnet_py dependency could not be installed that way)

v1.3.5 (March 2020)
- Bugfix for covariates samples mismatching for LMM
- Bugfix for numeric sample names
- Bugfix for crashes happening during enet variants filtering
- Fix pip and setup.py installs (glmnet_py)
- Remove deprecation warnings for fastlmm

v1.3.4 (December 2019)
- Update installation from pip and setup.py
- Allow multiple intervals for a single burden tests

v1.3.3 (September 2019)
- Fix for multiple reported alignments from bwa

v1.3.2 (September 2019)
- Features, bugfixes and documentation for prediction
- Bugfix for multiple hits/overlapping regions in annotate hits

v1.3.1 (July 2019)
- Handle API changes between different statsmodels versions
- Handle corner cases where no distance matrix is provided but --no-distances should be used

v1.3.0 (June 2019)
- Handle missing variants as missing (with new `--max-missing` option to control maximum allowed)
- Unit testing for enet functionality
- Added option to perform midpoint rooting on tree used by `phylogeny_distance.py`
- Updated docs with a new section about the use of unitigs
- Various bugfixes
   - Checks on phenotypes file
   - Improved BWA command and results handling
   - New optimizer for `fit_null`
   - Avoid deprecation warnings from pandas
   - Report to user if multiple chromosome/contigs are found while mapping
   - Use the correct samples order when using lineages and covariates

v1.2.0 (November 2018)
- Added elastic net functionality
- Avoid exiting with error in python 3.7+
- Allow spaces in Rtab's COG names
- Make unit tests compatible with python 2

v1.1.2 (October 2018)
- Bugfix for sample names interpreted as integers
- Bugfix for user-defined lineages
- Small fixes in docs and README

v1.1.1 (June 2018)
- Further bugfixes for the annotation of kmers
- Platform specific pandas bugfix with LMM
- Test kmers annotation
- Improved README

v1.1.0 (May 2018)
- New option (--no-distances) to run fixed effects associations without a distance matrix
- New option (--print-filtered) allowing to print variants not passing filters
- Check for consistency between lmm and lineage samples
- Add lineage to phandango mapping output
- annotate_hits_pyseer: properly parse GFF3 entries
- annotate_hits_pyseer: map all possible kemrs for reference genomes
- Update docs
- Tests are less dependent on machine precision

v1.0.2 (Mar 2018)
- Better dependencies definition

v1.0.1 (Mar 2018)
- Check for missing values in VCF entries

v1.0.0 (Feb 2018)
- First release

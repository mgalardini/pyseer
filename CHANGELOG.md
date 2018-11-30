# CHANGELOG
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

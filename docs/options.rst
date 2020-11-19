Option reference
================

Usage::

   usage: pyseer [-h] --phenotypes PHENOTYPES
                 [--phenotype-column PHENOTYPE_COLUMN]
                 (--kmers KMERS | --vcf VCF | --pres PRES) [--burden BURDEN]
                 [--distances DISTANCES | --load-m LOAD_M]
                 [--similarity SIMILARITY | --load-lmm LOAD_LMM]
                 [--save-m SAVE_M] [--save-lmm SAVE_LMM]
                 [--mds {classic,metric,non-metric}]
                 [--max-dimensions MAX_DIMENSIONS] [--no-distances]
                 [--continuous] [--lmm] [--wg {enet,rf,blup}] [--lineage]
                 [--lineage-clusters LINEAGE_CLUSTERS]
                 [--lineage-file LINEAGE_FILE] [--sequence-reweighting]
                 [--save-vars SAVE_VARS] [--load-vars LOAD_VARS]
                 [--save-model SAVE_MODEL] [--alpha ALPHA] [--n-folds N_FOLDS]
                 [--min-af MIN_AF] [--max-af MAX_AF] [--max-missing MAX_MISSING]
                 [--filter-pvalue FILTER_PVALUE] [--lrt-pvalue LRT_PVALUE]
                 [--cor-filter COR_FILTER] [--covariates COVARIATES]
                 [--use-covariates [USE_COVARIATES [USE_COVARIATES ...]]]
                 [--print-samples] [--print-filtered]
                 [--output-patterns OUTPUT_PATTERNS] [--uncompressed] [--cpu CPU]
                 [--block_size BLOCK_SIZE] [--version]

   SEER (doi: 10.1038/ncomms12797), reimplemented in python

   optional arguments:
     -h, --help            show this help message and exit

   Phenotype:
     --phenotypes PHENOTYPES
                           Phenotypes file (whitespace separated)
     --phenotype-column PHENOTYPE_COLUMN
                           Phenotype file column to use [Default: last column]

   Variants:
     --kmers KMERS         Kmers file
     --vcf VCF             VCF file. Will filter any non 'PASS' sites
     --pres PRES           Presence/absence .Rtab matrix as produced by roary and
                           piggy
     --burden BURDEN       VCF regions to group variants by for burden testing
                           (requires --vcf). Requires vcf to be indexed

   Distances:
     --distances DISTANCES
                           Strains distance square matrix (fixed or lineage
                           effects)
     --load-m LOAD_M       Load an existing matrix decomposition
     --similarity SIMILARITY
                           Strains similarity square matrix (for --lmm)
     --load-lmm LOAD_LMM   Load an existing lmm cache
     --save-m SAVE_M       Prefix for saving matrix decomposition
     --save-lmm SAVE_LMM   Prefix for saving LMM cache
     --mds {classic,metric,non-metric}
                           Type of multidimensional scaling [Default: classic]
     --max-dimensions MAX_DIMENSIONS
                           Maximum number of dimensions to consider after MDS
                           [Default: 10]
     --no-distances        Allow run without a distance matrix

   Association options:
     --continuous          Force continuous phenotype [Default: binary auto-
                           detect]
     --lmm                 Use random instead of fixed effects to correct for
                           population structure. Requires a similarity matrix
     --wg {enet,rf,blup}   Use a whole genome model for association and
                           prediction. Population structure correction is
                           implicit.
     --lineage             Report lineage effects
     --lineage-clusters LINEAGE_CLUSTERS
                           Custom clusters to use as lineages [Default: MDS
                           components]
     --lineage-file LINEAGE_FILE
                           File to write lineage association to [Default:
                           lineage_effects.txt]

   Whole genome options:
     --sequence-reweighting
                           Use --lineage-clusters to downweight sequences.
     --save-vars SAVE_VARS
                           Prefix for saving variants
     --load-vars LOAD_VARS
                           Prefix for loading variants
     --save-model SAVE_MODEL
                           Prefix for saving model
     --alpha ALPHA         Set the mixing between l1 and l2 penalties [Default:
                           0.0069]
     --n-folds N_FOLDS     Number of folds cross-validation to perform [Default:
                           10]

   Filtering options:
     --min-af MIN_AF       Minimum AF [Default: 0.01]
     --max-af MAX_AF       Maximum AF [Default: 0.99]
     --max-missing MAX_MISSING
                           Maximum missing (vcf/Rtab) [Default: 0.05]
     --filter-pvalue FILTER_PVALUE
                           Prefiltering t-test pvalue threshold [Default: 1]
     --lrt-pvalue LRT_PVALUE
                           Likelihood ratio test pvalue threshold [Default: 1]
     --cor-filter COR_FILTER
                           Correlation filter for elastic net (phenotype/variant
                           correlation quantile at which to start keeping
                           variants) [Default: 0.25]

   Covariates:
     --covariates COVARIATES
                           User-defined covariates file (tab-delimited, with
                           header, first column contains sample names)
     --use-covariates [USE_COVARIATES [USE_COVARIATES ...]]
                           Covariates to use. Format is "2 3q 4" (q for
                           quantitative) [Default: load covariates but don't use
                           them]

   Other:
     --print-samples       Print sample lists [Default: hide samples]
     --print-filtered      Print filtered variants (i.e. fitting errors)
                           [Default: hide them]
     --output-patterns OUTPUT_PATTERNS
                           File to print patterns to, useful for finding pvalue
                           threshold
     --uncompressed        Uncompressed kmers file [Default: gzipped]
     --cpu CPU             Processes [Default: 1]
     --block_size BLOCK_SIZE
                           Number of variants per core [Default: 3000]
     --version             show program's version number and exit

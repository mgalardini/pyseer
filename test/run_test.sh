#!/bin/bash

green="\033[1;32m"
red="\033[1;31m"
reset="\033[0m"

die () {
        echo -e $red"############"$reset
	echo -e $red$1$reset
	echo -e $red"Test failed!"$reset
	echo -e $red"############"$reset
	exit 1
}

# unpack the test data
zcat distances.tsv.gz > distances.tsv

# test all command line options
../pyseer kmers.gz subset.pheno distances.tsv --save-m pop_struct > 1.log 2> 1.err || die "Save population structure"
../pyseer kmers.gz subset.pheno distances.tsv --load-m pop_struct.pkl > 2.log 2> 2.err || die "Load population structure"
../pyseer kmers.gz subset.pheno distances.tsv --filter-pvalue 1E-5 --lrt-pvalue 1E-8 --load-m pop_struct.pkl > 3.log 2> 3.err || die "Basic filters"
../pyseer kmers.gz example.pheno distances.tsv --max-dimensions 3 --min-af 0.4 --max-af 0.6 > 4.log 2> 4.err || die "Binary phenotype w/ AF filtering"
../pyseer kmers.gz subset.cont.pheno distances.tsv --max-dimensions 3 --load-m pop_struct.pkl > 5.log 2> 5.err || die "Continuous phenotype"
../pyseer kmers.gz subset.pheno distances.tsv --max-dimensions 3 --continuous --load-m pop_struct.pkl > 6.log 2> 6.err || die "Force continuous phenotype"
../pyseer kmers.gz subset.pheno distances.tsv --max-dimensions 3 --print-samples --load-m pop_struct.pkl > 7.log 2> 7.err || die "Print samples"
../pyseer kmers.gz subset.pheno distances.tsv --load-m pop_struct.pkl > 8.log 2> 8.err || die "Use whole population structure"
../pyseer kmers.gz subset.pheno distances.tsv --max-dimensions 3 --covariates covariates.txt --use-covariates 2q 3 --load-m pop_struct.pkl > 9.log 2> 9.err || die "Use covariates"
../pyseer kmers.txt subset.pheno distances.tsv --uncompressed --load-m pop_struct.pkl > 10.log 2> 10.err || die "Uncompressed kmers"
../pyseer kmers.gz subset.pheno distances.tsv --cpu 2 --load-m pop_struct.pkl > 11.log 2> 11.err || die "Multiple cores"

# test all command line options (things that should fail or behave weirdly)
../pyseer kmers.txt subset.pheno distances.tsv --load-m pop_struct.pkl > /dev/null 2> /dev/null && die "Uncompressed kmers but no option"
../pyseer kmers.gz subset.pheno distances.tsv --load-m pop_struct.pkl --max-dimensions 1000 > /dev/null 2> /dev/null && die "Too many dimensions requested"
../pyseer kmers.gz subset.pheno distances.tsv --load-m pop_struct.pkl --max-dimensions 0 > /dev/null 2> /dev/null && die "Too few dimensions requested"
../pyseer kmers.gz subset.pheno distances.tsv --load-m pop_struct.pkl --covariates covariates.txt --use-covariates 10 > /dev/null 2> /dev/null && die "Incorrect covariate column ID"
../pyseer kmers.gz subset.pheno distances.tsv --load-m kmers.txt > /dev/null 2> /dev/null && die "Bogus population structure"
../pyseer kmers.gz supersubset.pheno distances.tsv --load-m pop_struct.pkl > /dev/null 2> /dev/null && die "Null model failure"
../pyseer kmers.gz supersubset.cont.pheno distances.tsv --load-m pop_struct.pkl > /dev/null 2> /dev/null || die "Weak results for continuous phenotype"
../pyseer kmers.gz monosubset.pheno distances.tsv --load-m pop_struct.pkl > /dev/null 2> /dev/null && die "Extreme skewed in binary phenotypes"

# Now compare the outputs
for t in $(seq 1 11);
do
  python compare_tests $t.log $t.err baseline/$t.log baseline/$t.err || die "Baseline comparison failed for $t";
done

echo -e $green"All good!"$reset

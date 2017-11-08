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

zcat distances.tsv.gz > distances.tsv

../pyseer kmers.gz subset.pheno distances.tsv --filter-pvalue 1E-5 --lrt-pvalue 1E-8 > 1.log 2> 1.err || die "Basic filters"
../pyseer kmers.gz example.pheno distances.tsv --filter-pvalue 1 --lrt-pvalue 1 --max-dimensions 3 > 2.log 2> 2.err || die "Binary phenotype"
../pyseer kmers.gz subset.cont.pheno distances.tsv --filter-pvalue 1 --lrt-pvalue 1 --max-dimensions 3 > 3.log 2> 3.err || die "Continuous phenotype"
../pyseer kmers.gz subset.pheno distances.tsv --filter-pvalue 1 --lrt-pvalue 1 --max-dimensions 3 --continuous > 4.log 2> 4.err || die "Force continuous phenotype"
../pyseer kmers.gz subset.pheno distances.tsv --filter-pvalue 1 --lrt-pvalue 1 --max-dimensions 3 --print-samples > 5.log 2> 5.err || die "Print samples"
../pyseer kmers.gz subset.pheno distances.tsv --filter-pvalue 1 --lrt-pvalue 1 > 6.log 2> 6.err || die "Use whole population structure"
../pyseer kmers.gz subset.pheno distances.tsv --filter-pvalue 1 --lrt-pvalue 1 --max-dimensions 3 --covariates covariates.txt --use-covariates 2q 3 > 7.log 2> 7.err || die "Use covariates"

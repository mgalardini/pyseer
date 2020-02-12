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
gzip -d -c presence_absence.Rtab.gz > presence_absence.Rtab

# test all command line options
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --distances distances.tsv.gz --save-m pop_struct > 1.log 2> 1.err || die "Save population structure"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --load-m pop_struct.pkl > 2.log 2> 2.err || die "Load population structure"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --filter-pvalue 1E-5 --lrt-pvalue 1E-8 --load-m pop_struct.pkl > 3.log 2> 3.err || die "Basic filters"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes example.pheno --distances distances.tsv.gz --max-dimensions 3 --min-af 0.4 --max-af 0.6 > 4.log 2> 4.err || die "Binary phenotype w/ AF filtering"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --max-dimensions 3 --load-m pop_struct.pkl --phenotype-column continuous > 5.log 2> 5.err || die "Continuous phenotype"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --max-dimensions 3 --continuous --load-m pop_struct.pkl > 6.log 2> 6.err || die "Force continuous phenotype"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --max-dimensions 3 --print-samples --load-m pop_struct.pkl > 7.log 2> 7.err || die "Print samples"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --load-m pop_struct.pkl > 8.log 2> 8.err || die "Use whole population structure"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --max-dimensions 3 --covariates covariates.txt --use-covariates 2q 3 --load-m pop_struct.pkl > 9.log 2> 9.err || die "Use covariates"
python ../pyseer-runner.py --kmers kmers.txt --phenotypes subset.pheno --uncompressed --load-m pop_struct.pkl > 10.log 2> 10.err || die "Uncompressed kmers"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --cpu 2 --load-m pop_struct.pkl > 11.log 2> 11.err || die "Multiple cores"
python ../pyseer-runner.py --vcf variants.vcf.gz --phenotypes subset.pheno --load-m pop_struct.pkl --max-dimensions 3 > 12.log 2> 12.err || die "VCF input"
python ../pyseer-runner.py --vcf variants.vcf.gz --burden burden_regions.txt --phenotypes subset.pheno --load-m pop_struct.pkl --max-dimensions 3 > 13.log 2> 13.err || die "VCF w/ burden test"
python ../pyseer-runner.py --pres presence_absence.Rtab --phenotypes subset.pheno --load-m pop_struct.pkl --max-dimensions 3 > 14.log 2> 14.err || die "Roary/piggy input"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --distances distances.tsv.gz --max-dimensions 3 --mds classic --continuous > 15.log 2> 15.err || die "Classic MDS"
# ugly hack for sklearn random_state
export PYSEERSEED="42"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --distances distances.tsv.gz --max-dimensions 3 --mds metric --continuous > 16.log 2> 16.err || die "Metric MDS"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --distances distances.tsv.gz --max-dimensions 3 --mds non-metric --continuous > 17.log 2> 17.err || die "Non-metric MDS"
#
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --distances distances.tsv.gz --max-dimensions 3 --lineage --lineage-file lineage.txt > 18.log 2> 18.err || die "Lineage effects"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --distances distances.tsv.gz --max-dimensions 3 --lineage --lineage-clusters lineage_clusters.txt > 19.log 2> 19.err || die "Lineage effects w/ user-provided clusters"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --similarity similarity.tsv.gz --lmm --save-lmm lmm.cache > 20.log 2> 20.err || die "LMM saving cache"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --lmm --load-lmm lmm.cache.npz > 21.log 2> 21.err || die "LMM loading cache"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --lmm --load-lmm lmm.cache.npz --lineage --load-m pop_struct.pkl > 22.log 2> 22.err || die "LMM with lineage"
python ../pyseer-runner.py --vcf variants.vcf.gz --phenotypes subset.pheno --lmm --load-lmm lmm.cache.npz > 23.log 2> 23.err || die "LMM with VCF input"
python ../pyseer-runner.py --pres presence_absence.Rtab --phenotypes subset.pheno --lmm --load-lmm lmm.cache.npz > 24.log 2> 24.err || die "LMM with roary/piggy input"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --lmm --load-lmm lmm.cache.npz --covariates covariates.txt --use-covariates 2q 3 > 25.log 2> 25.err || die "LMM with covariates"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --lmm --load-lmm lmm.cache.npz --cpu 2 > 26.log 2> 26.err || die "LMM with multiple CPUs"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --lmm --load-lmm lmm.cache.npz --output-patterns patterns.txt > 27.log 2> 27.err || die "Output patterns"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --no-distances > 28.log 2> 28.err || die "No distances"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --no-distances --use-covariates 3 --covariates covariates.txt > 29.log 2> 29.err || die "No distances and covariates"
python ../pyseer-runner.py --kmers kmers_int.gz --phenotypes subset_int.pheno --distances distances_int.tsv.gz > 30.log 2> 30.err || die "Sample names are all integers"
python ../pyseer-runner.py --vcf variants.vcf.gz --phenotypes subset.pheno --save-vars enet_vcf --wg enet --alpha 1 --cor-filter 0.25 > 31.log 2> 31.err || die "Enet with VCF input"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --wg enet --alpha 1 --cor-filter 0.25 > 32.log 2> 32.err || die "Load Enet with kmers input"
python ../pyseer-runner.py --pres presence_absence.Rtab --phenotypes subset.pheno --wg enet --alpha 1 --cor-filter 0.25 > 33.log 2> 33.err || die "Load Enet with roary/piggy input"
python ../pyseer-runner.py --vcf variants.vcf.gz --phenotypes subset.pheno --load-vars enet_vcf --wg enet --save-model enet_vcf_model --alpha 1 --cor-filter 0.25 > 34.log 2> 34.err || die "Load Enet and save model"
python ../pyseer-runner.py --vcf variants.vcf.gz --burden burden_regions_multiple.txt --phenotypes subset.pheno --load-m pop_struct.pkl --max-dimensions 3 > 35.log 2> 35.err || die "Multiple regions for burden testing"

# test other pyseer commands
python ../scree_plot_pyseer-runner.py distances.tsv.gz --max-dimensions 20 > /dev/null 2> /dev/null || die "Scree plot"
python ../similarity-runner.py samples.txt --kmers kmers.gz > /dev/null 2> /dev/null || die "Similarity w/ kmers"
python ../similarity-runner.py samples.txt --vcf variants.vcf.gz > /dev/null 2> /dev/null || die "Similarity w/ vcf"
python ../similarity-runner.py samples.txt --pres presence_absence.Rtab > /dev/null 2> /dev/null || die "Similarity w/ roary/piggy"
cat mash.tsv | python ../square_mash-runner.py > /dev/null 2> /dev/null || die "Mash squarer"
python ../annotate_hits_pyseer-runner.py significant_kmers.txt references.txt /dev/null > /dev/null 2> /dev/null || die 'Annotate hits'
python ../phandango_mapper-runner.py significant_kmers.txt Spn23F.fa phandango.test.out > /dev/null 2> /dev/null || die 'Phandango mapper'
python ../enet_predict-runner.py --vcf variants.vcf.gz enet_vcf_model.pkl subset.samples_list > /dev/null 2> /dev/null || die "Enet predict"

# test the scripts folder
python ../enet_predict-runner.py --vcf variants.vcf.gz enet_vcf_model.pkl subset.samples_list > /dev/null 2> /dev/null || die "Enet predict"
python ../scripts/count_patterns.py patterns.txt > /dev/null 2> /dev/null || die "Count patterns"
python ../scripts/phylogeny_distance.py tree.nwk > /dev/null 2> /dev/null || die "Tree distances"
python ../scripts/phylogeny_distance.py tree.nwk --lmm > /dev/null 2> /dev/null || die "Tree distances (C)"
python ../scripts/phylogeny_distance.py tree.nwk --topology > /dev/null 2> /dev/null || die "Tree distances (topology)"
python ../pyseer-runner.py --vcf variants.vcf.gz --phenotypes subset.pheno --lmm --similarity similarity.tsv.gz > vcf.lmm.txt 2> /dev/null || die "LMM input for enet model"
python ../scripts/save_model.py --p-cutoff 0.5 vcf.lmm.txt enet.lmm > /dev/null 2> /dev/null || die "Save enet model from LMM input"
python ../enet_predict-runner.py --vcf variants.vcf.gz enet.lmm.pkl subset.samples_list > /dev/null 2> /dev/null || die "Enet predict with LMM model"

# test all command line options (things that should fail or behave weirdly)
python ../pyseer-runner.py --kmers kmers.txt --phenotypes subset.pheno --load-m pop_struct.pkl > /dev/null 2> /dev/null && die "Uncompressed kmers but no option"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --load-m pop_struct.pkl --max-dimensions 1000 > /dev/null 2> /dev/null && die "Too many dimensions requested"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --load-m pop_struct.pkl --max-dimensions 0 > /dev/null 2> /dev/null && die "Too few dimensions requested"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --load-m pop_struct.pkl --covariates covariates.txt --use-covariates 10 > /dev/null 2> /dev/null && die "Incorrect covariate column ID"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --load-m kmers.txt > /dev/null 2> /dev/null && die "Bogus population structure"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes supersubset.pheno --load-m pop_struct.pkl > /dev/null 2> /dev/null && die "Null model failure"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes supersubset.pheno --load-m pop_struct.pkl --phenotype-column continuous > /dev/null 2> /dev/null || die "Weak results for continuous phenotype"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes monosubset.pheno --load-m pop_struct.pkl > /dev/null 2> /dev/null && die "Extremely skewed binary phenotypes"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --no-distances --lmm > /dev/null 2> /dev/null && die "No distances but LMM"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --no-distances --load-m pop_struct.pkl > /dev/null 2> /dev/null && die "No distances but distances provided"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --similarity similarity.tsv.gz --lmm --covariates covariates_missing.txt --use-covariates 2q --save-lmm lmm.cache > /dev/null 2> /dev/null && die "LMM with non-matching covariates"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --lmm --covariates covariates_missing.txt --use-covariates 2q --load-lmm lmm.cache.npz > /dev/null 2> /dev/null && die "LMM reloaded with non-matching covariates"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --similarity similarity.tsv.gz --lmm --covariates covariates_mismatch.txt --use-covariates 2q --save-lmm lmm.cache > /dev/null 2> /dev/null && die "LMM with non-matching covariates"
python ../pyseer-runner.py --kmers kmers.gz --phenotypes subset.pheno --lmm --covariates covariates_mismatch.txt --use-covariates 2q --load-lmm lmm.cache.npz > /dev/null 2> /dev/null && die "LMM reloaded with non-matching covariates"
# TODO: test enet options failures

# Now compare the outputs
for t in $(seq 1 35);
do
  echo "Comparing results and error messages to baseline "$t;
  python compare_tests $t.log $t.err baseline/$t.log baseline/$t.err || die "Baseline comparison failed for $t";
done

echo -e $green"All good!"$reset

#!/bin/sh -l

#SBATCH --nodes=1              ##number of requested nodes (a node is a group of CPUs)
#SBATCH --tasks-per-node=32    ##number of requested CPUs
#SBATCH --time=1-00:00:00      ##time requested
#SBATCH --job-name meta4.0     ##name of job
#SBATCH -A microbiome          ##name of the queue you are using

###############################################################################
#step 1, cd into the proper directory.

cd /scratch/bell/dong316/amplicon4/qiime_out
pwd

#step 2, Load Modules

echo "loading bioinfo"
module load bioinfo
echo "loading qiime2-2021.11"
module load Qiime/2-2021.11

#step 3, Run the code you want

echo "Start time"
date +"%d %B %Y %H:%M:%S"
echo
# echo "##############################"
# echo "Data Import"
# # Make the sequences into a qiime recognizable format

# time qiime tools import \
  # --type 'SampleData[SequencesWithQuality]' \
  # --input-path ../manifest1 \
  # --output-path single-end-demux.qza \
  # --input-format SingleEndFastqManifestPhred33V2

# time qiime tools import \
# --type 'SampleData[PairedEndSequencesWithQuality]' \
# --input-path ../manifest2 \
# --output-path paired-end-demux.qza \
# --input-format PairedEndFastqManifestPhred33V2

# echo "Data Import End"
# echo "##############################" 
# echo

# echo "##############################"
# echo "demux summarize visualization"
# # 1 for single-end sequence, 2 for paired-end sequence

# mv single-end-demux.qza demux1.qza
# mv paired-end-demux.qza demux2.qza

# time qiime demux summarize \
  # --i-data demux1.qza \
  # --o-visualization demux1.qzv
  
# time qiime demux summarize \
  # --i-data demux2.qza \
  # --o-visualization demux2.qzv
  
# echo "visualization end" 
# echo "##############################"
# echo

# echo "##############################"
# echo "merge paired reads using q2-vsearch join-paires"

# time qiime vsearch join-pairs \
  # --i-demultiplexed-seqs demux2.qza \
  # --p-minovlen 20 \
  # --p-minmergelen 200 \
  # --p-maxmergelen 600 \
  # --p-threads 8 \
  # --p-allowmergestagger \
  # --o-joined-sequences demux-joined2.qza
  
# echo "join pair visualization" 

# time qiime demux summarize \
  # --i-data demux-joined2.qza \
  # --o-visualization demux-joined2.qzv
# echo "visualization paired end"
# echo "##############################"
# echo

# echo "##############################"
# echo "Sequence quality control"

# echo "single"
# time qiime quality-filter q-score \
 # --i-demux demux1.qza \
 # --p-min-quality 20 \
 # --o-filtered-sequences demux-filtered1.qza \
 # --o-filter-stats demux-filter-stats1.qza
 
# echo "paired"
# time qiime quality-filter q-score \
  # --i-demux demux-joined2.qza \
  # --p-min-quality 20 \
  # --o-filtered-sequences demux-filtered2.qza \
  # --o-filter-stats demux-joined-filter-stats2.qza
  
# echo "quality-filter summarize single"  
# time qiime demux summarize \
  # --i-data demux-filtered1.qza \
  # --o-visualization demux-filterd1.qzv
  
# echo "quality-filter summarize paired"  
# time qiime demux summarize \
  # --i-data demux-filtered2.qza \
  # --o-visualization demux-filterd2.qzv

# echo "Sequence quality control end"
# echo "##############################"
# echo

# echo "##############################"
# echo "dereplicating"

# cp demux-filtered1.qza seqs1.qza
# cp demux-filtered2.qza seqs2.qza

# echo "derep single"
# time qiime vsearch dereplicate-sequences \
  # --i-sequences seqs1.qza \
  # --o-dereplicated-table table1.qza \
  # --o-dereplicated-sequences rep-seqs1.qza
  
# echo "derep paired"  
# time qiime vsearch dereplicate-sequences \
  # --i-sequences seqs2.qza \
  # --o-dereplicated-table table2.qza \
  # --o-dereplicated-sequences rep-seqs2.qza
# date
# echo "dereplicating end"
# echo "##############################"
# echo

# echo "##############################"
# echo "visualization after derep"

# time qiime feature-table summarize \
  # --i-table table1.qza \
  # --o-visualization table1.qzv \
  # --m-sample-metadata-file ../sample-metadata.tsv

# time qiime feature-table tabulate-seqs \
  # --i-data rep-seqs1.qza \
  # --o-visualization rep-seqs1.qzv
  
# time qiime feature-table summarize \
  # --i-table table2.qza \
  # --o-visualization table2.qzv \
  # --m-sample-metadata-file ../sample-metadata.tsv

# time qiime feature-table tabulate-seqs \
  # --i-data rep-seqs2.qza \
  # --o-visualization rep-seqs2.qzv
  
# echo
# date +"%d %B %Y %H:%M:%S"
# echo "visualization end after derep"
# echo "##############################"
# echo

# echo "##############################"
# echo "reference import"

# cd /depot/john2185/data/dong316/meta_analysis/test

# wget "https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip"
# unpigz Silva_132_release.zip

# SILVA_DIR=/depot/john2185/data/dong316/meta_analysis/test/Silva132

# qiime tools import \
 # --input-path ${SILVA_DIR}/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna \
 # --output-path ${SILVA_DIR}/97_otus.qza \
 # --type 'SampleData[Sequence]'

# qiime tools import \
 # --type 'FeatureData[Taxonomy]' \
 # --input-format HeaderlessTSVTaxonomyFormat \
 # --input-path ${SILVA_DIR}/SILVA_132_QIIME_release/taxonomy/16S_only/97/taxonomy_7_levels.txt \
 # --output-path ${SILVA_DIR}/taxonomy_97.qza

# qiime tools import \
 # --input-path ${SILVA_DIR}/SILVA_132_QIIME_release/trees/97/97_otus.tre \
 # --output-path ${SILVA_DIR}/tree_97.qza \
 # --type 'Phylogeny[Rooted]'
 
# cd -
# cp ${SILVA_DIR}/*97* .

# echo "reference import"
# echo "##############################"

# echo "##############################"
# echo "closed  reference clustering"
# # cluster at 97% identity against silva132 97% OTUs reference database

# echo "single"
# time qiime vsearch cluster-features-closed-reference \
 # --i-table table1.qza \
# --i-sequences rep-seqs1.qza \
 # --i-reference-sequences ../97_otus.qza \
 # --p-perc-identity 0.97 \
 # --p-threads 32 \
 # --o-clustered-table table1-cr-97.qza \
 # --o-clustered-sequences rep-seqs1-cr-97.qza \
 # --o-unmatched-sequences unmatched1-cr-97.qza

# echo "paird"  
# time qiime vsearch cluster-features-closed-reference \
  # --i-table table2.qza \
  # --i-sequences rep-seqs2.qza \
  # --i-reference-sequences ../97_otus.qza \
  # --p-perc-identity 0.97 \
  # --p-threads 32 \
  # --o-clustered-table table2-cr-97.qza \
  # --o-clustered-sequences rep-seqs2-cr-97.qza \
  # --o-unmatched-sequences unmatched2-cr-97.qza
  
# echo
# echo "clustering end"
# echo "##############################"
# echo

# echo "##############################"  
# echo "visualization after closed reference"

# echo "single"
# time qiime feature-table summarize \
  # --i-table table1-cr-97.qza \
  # --o-visualization table1-cr-97.qzv \
  # --m-sample-metadata-file ../sample-metadata.tsv

# time qiime feature-table tabulate-seqs \
  # --i-data rep-seqs1-cr-97.qza \
  # --o-visualization rep-seqs1-cr-97.qzv
  
# echo
 
# echo "paired" 
# time qiime feature-table summarize \
  # --i-table table2-cr-97.qza \
  # --o-visualization table2-cr-97.qzv \
  # --m-sample-metadata-file ../sample-metadata.tsv
  
# time qiime feature-table tabulate-seqs \
  # --i-data rep-seqs2-cr-97.qza \
  # --o-visualization rep-seqs2-cr-97.qzv
# echo

# echo "visualization end after closed reference"
# echo "##############################"
# echo

# echo "##############################" 
# echo "merge"
# time qiime feature-table merge \
  # --i-tables table1-cr-97.qza \
  # --i-tables table2-cr-97.qza \
  # --o-merged-table table3.qza
  
# echo "merge rep seqs"
# time qiime feature-table merge-seqs \
  # --i-data rep-seqs1-cr-97.qza \
  # --i-data rep-seqs2-cr-97.qza \
  # --o-merged-data rep-seqs3.qza
# echo

# echo "merge end"
# echo "##############################" 
# echo

# echo "##############################"
# echo "visualization after merge"

# time qiime feature-table summarize \
  # --i-table table3.qza \
  # --o-visualization table3.qzv \
  # --m-sample-metadata-file ../sample-metadata.tsv

# time qiime feature-table tabulate-seqs \
  # --i-data rep-seqs3.qza \
  # --o-visualization rep-seqs3.qzv
# echo

# echo "visualization end after merge"
# echo "##############################"
# echo

# echo "##############################" 
# echo "Run de novo chimera checking"
# time qiime vsearch uchime-denovo \
  # --i-table table3.qza \
  # --i-sequences rep-seqs3.qza \
  # --output-dir uchime-dn-out
# echo

# echo "Visualize summary stats"
# time qiime metadata tabulate \
  # --m-input-file uchime-dn-out/stats.qza \
  # --o-visualization uchime-dn-out/stats.qzv
# echo

# echo
# echo "Exclude chimeras but retain borderline chimeras"
# time qiime feature-table filter-features \
  # --i-table table3.qza \
  # --m-metadata-file uchime-dn-out/chimeras.qza \
  # --p-exclude-ids \
  # --o-filtered-table uchime-dn-out/table-nonchimeric-w-borderline.qza
  
# time qiime feature-table filter-seqs \
  # --i-data rep-seqs3.qza \
  # --m-metadata-file uchime-dn-out/chimeras.qza \
  # --p-exclude-ids \
  # --o-filtered-data uchime-dn-out/rep-seqs-nonchimeric-w-borderline.qza
  
# echo

# echo "Exclude chimeras but retain borderline chimeras end"
# echo "dechimera end"
# echo "##############################" 
# echo
# echo "##############################" 
# echo "visualization final"
# cp uchime-dn-out/table-nonchimeric-w-borderline.qza table.qza
# cp uchime-dn-out/rep-seqs-nonchimeric-w-borderline.qza rep-seqs.qza

# time qiime feature-table summarize \
  # --i-table table.qza \
  # --o-visualization table.qzv \
  # --m-sample-metadata-file ../sample-metadata.tsv

# time qiime feature-table tabulate-seqs \
  # --i-data rep-seqs.qza \
  # --o-visualization rep-seqs.qzv
# echo

# echo "visualization final end"
# echo "##############################" 
# echo

# echo "##############################"
# echo "rarefaction curve"

# time qiime diversity alpha-rarefaction \
  # --i-table table.qza \
  # --i-phylogeny ../tree_97.qza \
  # --p-max-depth 10000 \
  # --m-metadata-file ../sample-metadata.tsv \
  # --o-visualization alpha-rarefaction-10000.qzv

# time qiime diversity alpha-rarefaction \
  # --i-table table.qza \
  # --i-phylogeny ../tree_97.qza \
  # --p-max-depth 20000 \
  # --m-metadata-file ../sample-metadata.tsv \
  # --o-visualization alpha-rarefaction-20000.qzv

# time qiime diversity alpha-rarefaction \
  # --i-table table.qza \
  # --i-phylogeny ../tree_97.qza \
  # --p-max-depth 40000 \
  # --m-metadata-file ../sample-metadata.tsv \
  # --o-visualization alpha-rarefaction-40000.qzv

# echo "rarefaction curve end"
# echo "##############################"
# echo

# echo "##############################"
# echo "filter data"

# # samples with a total frequency less than 2000 will be filtered.
# time qiime feature-table filter-samples \
  # --i-table table.qza \
  # --p-min-frequency 2000 \
  # --o-filtered-table sample-frequency-filtered-table.qza

# # remove features with a total abundance (summed across all samples) of less than 10 as follows.
# time qiime feature-table filter-features \
  # --i-table sample-frequency-filtered-table.qza \
  # --p-min-frequency 10 \
  # --o-filtered-table feature-frequency-filtered-table.qza

# # features that are present in only a single sample will be filtered from a feature table.
# time qiime feature-table filter-features \
  # --i-table feature-frequency-filtered-table.qza \
  # --p-min-samples 2 \
  # --o-filtered-table sample-contingency-filtered-table.qza 

# # samples that contain only a few features will be filtered from a feature table.  
# time qiime feature-table filter-samples \
  # --i-table sample-contingency-filtered-table.qza \
  # --p-min-features 10 \
  # --o-filtered-table feature-contingency-filtered-table.qza

# # table-out.qza can be used to DMM cluster directly 
# mv feature-contingency-filtered-table.qza table-out.qza

# echo "filter data end"
# echo "##############################"
# echo

# echo "##############################"
# echo "calculate diversity metrics"

# to keep 90% samples  
# qiime diversity core-metrics-phylogenetic \
  # --i-phylogeny ../tree_97.qza \
  # --i-table table.qza \
  # --p-sampling-depth 7550 \
  # --m-metadata-file ../sample-metadata.tsv \
  # --p-n-jobs-or-threads 32 \
  # --output-dir core-metrics-results

# cd core-metrics-results

# time qiime diversity alpha \
  # --i-table rarefied_table.qza \
  # --p-metric goods_coverage \
  # --o-alpha-diversity goods_coverage.qza

# # alpha correlation
# time qiime diversity alpha-correlation \
  # --i-alpha-diversity faith_pd_vector.qza \
  # --m-metadata-file ../../sample-metadata.tsv \
  # --o-visualization alpha-corr-faith.qzv
  
# time qiime diversity alpha-correlation \
  # --i-alpha-diversity shannon_vector.qza \
  # --m-metadata-file ../../sample-metadata.tsv \
  # --o-visualization alpha-corr-shannon.qzv  

# time qiime diversity alpha-correlation \
  # --i-alpha-diversity observed_features_vector.qza \
  # --m-metadata-file ../../sample-metadata.tsv \
  # --o-visualization alpha-corr-observed.qzv 
  
# time qiime diversity alpha-correlation \
  # --i-alpha-diversity evenness_vector.qza \
  # --m-metadata-file ../../sample-metadata.tsv \
  # --o-visualization alpha-corr-evenness.qzv

# # metadata distance matrix
# time qiime metadata distance-matrix \
  # --m-metadata-file ../../sample-metadata.tsv \
  # --m-metadata-column age_day \
  # --o-distance-matrix sample-metadata-age.qza

# # correlation between metadata distance matrix and beta distance matrix
# time qiime diversity mantel \
  # --i-dm1 weighted_unifrac_distance_matrix.qza \
  # --i-dm2 sample-metadata-age.qza \
  # --p-method spearman \
  # --p-intersect-ids True \
  # --p-label1 weighted_unifrac \
  # --p-label2 age \
  # --o-visualization weighted_unifrac_age.qzv

# time qiime diversity bioenv \
  # --i-distance-matrix weighted_unifrac_distance_matrix.qza \
  # --m-metadata-file ../../sample-metadata.tsv \
  # --o-visualization weighted-unifrac-bioenv.qzv  
  
# time qiime diversity mantel \
  # --i-dm1 unweighted_unifrac_distance_matrix.qza \
  # --i-dm2 sample-metadata-age.qza \
  # --p-method spearman \
  # --p-intersect-ids True \
  # --p-label1 unweighted_unifrac \
  # --p-label2 age \
  # --o-visualization unweighted_unifrac_age.qzv 

# time qiime diversity bioenv \
  # --i-distance-matrix unweighted_unifrac_distance_matrix.qza \
  # --m-metadata-file ../../sample-metadata.tsv \
  # --o-visualization unweighted-unifrac-bioenv.qzv  
  
# time qiime diversity mantel \
  # --i-dm1 bray_curtis_distance_matrix.qza \
  # --i-dm2 sample-metadata-age.qza \
  # --p-method spearman \
  # --p-intersect-ids True \
  # --p-label1 bray_curtis \
  # --p-label2 age \
  # --o-visualization bray_curtis_age.qzv 
  
# time qiime diversity bioenv \
  # --i-distance-matrix bray_curtis_distance_matrix.qza \
  # --m-metadata-file ../../sample-metadata.tsv \
  # --o-visualization bray_curtis-bioenv.qzv
  
# time qiime diversity mantel \
  # --i-dm1 jaccard_distance_matrix.qza \
  # --i-dm2 sample-metadata-age.qza \
  # --p-method spearman \
  # --p-intersect-ids True \
  # --p-label1 jaccard \
  # --p-label2 age \
  # --o-visualization jaccard_age.qzv 
  
# time qiime diversity bioenv \
  # --i-distance-matrix jaccard_distance_matrix.qza \
  # --m-metadata-file ../../sample-metadata.tsv \
  # --o-visualization jaccard-bioenv.qzv

# echo "calculate diversity metrics end"
# echo "##############################"
# echo

# cd ..
# time qiime taxa barplot \
  # --i-table table-out.qza \
  # --i-taxonomy ../taxonomy_97.qza \
  # --m-metadata-file ../sample-metadata.tsv \
  # --o-visualization taxa-bar-plots.qzv

# echo "collapse"

# time qiime taxa collapse \
  # --i-table table.qza \
  # --i-taxonomy ../taxonomy_97.qza \
  # --p-level 6 \
  # --o-collapsed-table table-level6.qza

# time qiime taxa collapse \
  # --i-table table.qza \
  # --i-taxonomy ../taxonomy_97.qza \
  # --p-level 5 \
  # --o-collapsed-table table-level5.qza

# time qiime taxa collapse \
  # --i-table table.qza \
  # --i-taxonomy ../taxonomy_97.qza \
  # --p-level 4 \
  # --o-collapsed-table table-level4.qza

# time qiime taxa collapse \
  # --i-table table.qza \
  # --i-taxonomy ../taxonomy_97.qza \
  # --p-level 3 \
  # --o-collapsed-table table-level3.qza

# time qiime taxa collapse \
  # --i-table table.qza \
  # --i-taxonomy ../taxonomy_97.qza \
  # --p-level 2 \
  # --o-collapsed-table table-level2.qza
  
time qiime feature-table summarize \
  --i-table table-level6.qza \
  --o-visualization table-level6.qzv \
  --m-sample-metadata-file ../sample-metadata.tsv

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ../tree_97.qza \
  --i-table table-level6.qza \
  --p-sampling-depth 7550 \
  --m-metadata-file ../sample-metadata.tsv \
  --p-n-jobs-or-threads 32 \
  --output-dir core-metrics-results-level6  
 

echo "End time"
date +"%d %B %Y %H:%M:%S"
#to run: sbatch <file.name>

#!/usr/bin/env bash

source activate qiime2


cd  fastq_files
for i in *; do mv $i $(echo $i | awk -F'[_]' '{print $1}')_ASD_$(echo $i | awk -F'[_]' '{print $2 "_" $3  "_" $4}'); done
for i in *; do mv $i ${i}_001.fastq; done
for i in *; do gzip $i; done
cd ../

mkdir forward_reads
cp fastq_files/*R1* forward_reads

#import demultiplexed data
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path forward_reads/ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-single-end.qza
qiime demux summarize --i-data demux-single-end.qza --o-visualization demux.qzv

#run dada2
qiime dada2 denoise-single --i-demultiplexed-seqs demux-single-end.qza --o-table table --o-representative-sequences rep-seqs --p-trunc-len 230 --p-n-threads 100 --verbose --p-chimera-method consensus --output-dir data2-single
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file metadata.txt


#make tree of ref sequences
qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza
qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza


#CHOOSE YOUR SAMPLING DEPTH BY LOOKING AT table.qzv!!!
DEPTH=1300

#run diversity metrics
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza  --i-table table.qza --p-sampling-depth $DEPTH --output-dir core-metrics-results --m-metadata-file metadata.txt

#alpha diversity
qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny rooted-tree.qza --p-max-depth 50000 --m-metadata-file metadata.txt --o-visualization alpha-rarefaction.qzv
qiime diversity alpha-group-significance   --i-alpha-diversity core-metrics-results/faith_pd_vector.qza   --m-metadata-file metadata.txt   --o-visualization core-metrics-results/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance   --i-alpha-diversity core-metrics-results/evenness_vector.qza   --m-metadata-file metadata.txt   --o-visualization core-metrics-results/evenness-group-significance.qzv

#beta diversity significance
for f in core-metrics-results/*distance_matrix.qza; do
		qiime diversity beta-group-significance --i-distance-matrix $f --m-metadata-file metadata.txt --m-metadata-column Location --o-visualization ${f%.*}-significance.qzv --p-pairwise &
done


#taxonomy assignment and plots
qiime feature-classifier classify-sklearn --i-classifier ../Silva_classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.txt --o-visualization taxa-bar-plots.qzv


#get differential ASVs
qiime composition add-pseudocount --i-table table.qza --o-composition-table comp-table.qza
qiime composition ancom --i-table comp-table.qza --m-metadata-file metadata.txt --m-metadata-column Location --o-visualization ancom-asv.qzv

#get differential taxa
for t in 1 2 3 4 5 6; do
	qiime taxa collapse --i-table table.qza --i-taxonomy taxonomy.qza --p-level $t --o-collapsed-table table-d${t}.qza
	qiime composition add-pseudocount --i-table table-d${t}.qza --o-composition-table comp-table-d${t}.qza
	qiime composition ancom --i-table comp-table-d${t}.qza --m-metadata-file metadata.txt --m-metadata-column Location --o-visualization ancom-d${t}.qzv
done



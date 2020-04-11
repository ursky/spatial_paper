#!/usr/bin/env bash

source activate qiime2

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
DEPTH=5183

#run diversity metrics
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza  --i-table table.qza --p-sampling-depth $DEPTH --output-dir core-metrics-results --m-metadata-file metadata.txt


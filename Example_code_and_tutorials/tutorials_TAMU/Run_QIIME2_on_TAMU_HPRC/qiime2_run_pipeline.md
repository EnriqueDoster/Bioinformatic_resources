

# Classification

```
qiime feature-classifier classify-sklearn \
--i-classifier /scratch/group/vero_research/QIIME_Classifiers/qiime2-2023.2/SILVA138.1/341f-785r/silva-138.1-ssu-nr99-341f-785r-classifier.qza \
--i-reads merged-rep-seqs.qza --o-classification taxonomy.qza

qiime taxa filter-table --i-table merged-dada-table.qza --i-taxonomy taxonomy.qza \
--p-exclude mitochondria,chloroplast --o-filtered-table no-chloro-no-mito/merged-dada-table.qza

qiime taxa filter-seqs --i-sequences merged-rep-seqs.qza --i-taxonomy taxonomy.qza \
--p-exclude mitochondria,chloroplast \
--o-filtered-sequences no-chloro-no-mito/merged-rep-seqs.qza

cd no-mito-no-chloro/

qiime alignment mafft --i-sequences merged-rep-seqs.qza \
--o-alignment aligned-merged-rep-seqs.qza

qiime alignment mask --i-alignment aligned-merged-rep-seqs.qza \
--o-masked-alignment masked-aligned-merged-rep-seqs.qza

qiime phylogeny fasttree --i-alignment masked-aligned-merged-rep-seqs.qza \
--o-tree unrooted-tree.qza

qiime phylogeny midpoint-root --i-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

mkdir exported 

qiime tools export --input-path merged-rep-seqs.qza --output-path exported/                                 
qiime tools export --input-path taxonomy.qza --output-path exported/
qiime tools export --input-path rooted-tree.qza --output-path exported/
qiime tools export --input-path table.qza --output-path exported/


```




```

wget   -O "85_otus.fasta"   "https://data.qiime2.org/2020.2/tutorials/training-feature-classifiers/85_otus.fasta"

wget \
  -O "85_otu_taxonomy.txt" \
  "https://data.qiime2.org/2020.2/tutorials/training-feature-classifiers/85_otu_taxonomy.txt"
  
# add any sequences to the 85_otus.fasta file

# modify taxonomy file to include header name of new sequences and classification on the second column
  

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path spikein_85_otus.fasta \
  --output-path spikein_85_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path spikein_85_otu_taxonomy.txt \
  --output-path spikein_ref-taxonomy.qza

  qiime feature-classifier extract-reads \
  --i-sequences spikein_85_otus.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 120 \
  --p-min-length 100 \
  --p-max-length 400 \
  --o-reads spikein_ref-seqs.qza
  
  qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads spikein_ref-seqs.qza \
  --i-reference-taxonomy spikein_ref-taxonomy.qza \
  --o-classifier spikein_85gg_classifier_June2020.qza

```

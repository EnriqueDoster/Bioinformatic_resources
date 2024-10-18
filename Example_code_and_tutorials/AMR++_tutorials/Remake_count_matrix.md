# How to run the python script to merge results and create a new count matrix



From the AMRplusplus repository, you can run this command to make a new count matrix. 

```
python3 bin/amr_long_to_wide.py -i Grayson_TE_results/ResistomeAnalysis/ResistomeCounts/*.dedup_AMR.gene.tsv -o deduped_AMR_analytic_matrix.csv
```
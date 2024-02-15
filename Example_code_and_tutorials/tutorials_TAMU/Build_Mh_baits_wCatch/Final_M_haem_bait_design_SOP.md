# Bait enrichment for metagenomic sequencing of *Mannheimia haemolytica*


# Table of contents
* [Reference genomes](#reference-genomes)
* [Bait design](#bait-design-with-catch)
* [Bait design refinement](#bait-design-refinement)
  * [Increasing specificity with kraken2](#increasing-probe-specificity-with-kraken2)
  * [Informative loci](#informative-loci)
  * [Integrative and conjugative elements (ICE)](#ice-elements)
* [Final bait design overview](#final-bait-design-overview)


# Reference genomes
* 69 *M. haemolytica* genomes (CP017484 – CP017552)
  * Accession: PRJNA340884  
* CP006957

# Bait design with CATCH
Bait design software
* We employed the [CATCH software](https://github.com/broadinstitute/catch) to create the starting bait design.
   * Metsky HC and Siddle KJ et al. Capturing sequence diversity in metagenomes with comprehensive and scalable probe design. Nature Biotechnology, 37(2), 160–168 (2019). doi: 10.1038/s41587-018-0006-x
* We ran the "design.py" script on the 70 reference genomes using the following flags:
   * `-pl 120` = “Probe Length”, number of nucleotides for each probe
   * `-ps 120` = “Probe stride”, distance between probes
   * `-m 20` = “mismatches”, number of mismatches that can be tolerated
     * This was selected because the Agilent baits allow up to a 20 nucleotide mismatch from the target sequence for successful hybridization.
   * `-e 50` = “cover extension”, number of nucleotides on either side of the probe target that is assumed to be captured with the same DNA fragment.
     * This flag value is recommended by the CATCH developers.
   * `--cover-groupings-separately --cluster-from-fragments 10000 --max-num-processes 120` = flags to improve the speed of bait design (days vs weeks)

```
design.py all_70_combined_genomes.fasta --cluster-and-design-separately 0.5 -e 50 -pl 120 -ps 120 -m 20 -o M_haem_70genome_probes_pl120_ps120_e50_m20.fasta --verbose --cover-groupings-separately --cluster-from-fragments 10000 --max-num-processes 120 --write-analysis-to-tsv analysis_M_haem_70genome_probes_pl120_ps120_e50_m20 --write-sliding-window-coverage sliding_coverage_M_haem_70genome_probes_pl120_ps120_e50_m20
```
* **51,612 probes created** from the 70 reference genomes

# Bait design refinement

## Increasing probe specificity with Kraken2
We used kraken2 to taxonomically classify all the probes as if they were metagenomic reads. Then, only reads classified as *M. haemolytica* (or it's taxonomic lineage) were filtered out for inclusion in the final bait design.

```
# Run kraken2 with nt kraken
kraken2 --db /home/noyes046/shared/databases/kraken2_databases/nt_kraken_db_July2020 M_haem_70genome_probes_pl120_ps120_e50_m20.fasta --threads 20 --report kraken_report_M_haem_70genome_probes_pl120_ps120_e50_m20 > kraken_raw_M_haem_70genome_probes_pl120_ps120_e50_m20

# extract kraken-specific reads
python extract_kraken_reads.py -k kraken_raw_M_haem_70genome_probes_pl120_ps120_e50_m20 -s M_haem_70genome_probes_pl120_ps120_e50_m20.fasta -o all70_M_haem_blast_sp_probes.fasta -t 75985 --include-children --include-parents -r kraken_report_M_haem_70genome_probes_pl120_ps120_e50_m20

```
* **47,011 probes extracted from the kraken results (nt database)**



---
## Informative loci
To improve coverage of "informative loci" in the *M. haemolytica* genome, Dr. Mike Clawson helped us identify 17 loci for consideration.
* Unfortunately, we only found 13/17 loci
  * Could not find loci D650_20320, or BG559_12965 in the reference genome
  * Could not find 2 “intergenic” regions due to differences in coordinates

To improve the coverage of these regions, we used probes created on the 70 *M. haemolytica* reference genomes without any additional flags (n=103,024) and aligned them to the 13 informative loci using the GUI software, geneious.
* All 13 loci had a mean coverage > 3 baits at each nucleotide position.
* **Probes that aligned to these informative loci (n=1,361) were extracted and added to the bait design, now at 48,372 probes**


## ICE elements

Dr. Mike Clawson shared the sequences and flanking regions for 3 integrative ice elements (ICE) associated with *M. haemolytica*. To increase our coverage of any further natural sequence variation in ICE elements, Dr. Ilya Slizovskiy probed the corresponding annotations for genes associated with the ICE machinery.
* Then, these genes were analyzed using "blast" and NCBI's "nt" database to identify other sequence variants.
* Following the removal of any sequences matching the reference genes with < 80% sequence identity, **399 novel sequence variants** were extracted.
* The "design.py" script was used with the command below to create probes for the 3 ICE elements and 399 novel gene variants.

```
# Run blast to identify sequence variants
blastn -db nt -query M_haemolytica_ICE_machinery.fasta -max_target_seqs 20 -out top20_blast_M_haem_ICE_classified -outfmt "6 qseqid qlen sscinames sseq evalue bitscore score qstart qend length pident mismatch" -remote

# All three ICE sequences with added variants
design.py all_ICE_and_80per_extra_variants.fasta -pl 120 -ps 120 -o all_ICE_and_80per_extra_variants.fasta.probes --max-num-processes 6 --verbose
```
* **A total of 2,131 probes were created and added to the bait design, now at 50,503 total probes.**
  * Next step is to reduce the number of redundant probes by extracting baits from "M_haem_specific_probes_70genomes.fasta" mapping to the CP006957 reference genome, which contributed a disproportionate number of baits (~20k) to the design.



  # Final bait design overview

  Final changes to bait design:
  * Redundant probe sequences were removed using "CD-hit", resulting in a **bait design with 50,011 probes**
  * The bait design with 50,011 probes equates with a 6MB bait design, which is just above the threshold of 5.9MB for a lower pricing tier with Agilent.
  * Therefore, we used CD-hit and clustered probe sequences by 99% and 98% sequence identity to slightly reduce the bait design.

  ```
  # Remove redundant probes
  cd-hit-est -i all_combined_baits.fasta -o Final_M_haem_bait_design.fasta -c 1.0 -AS 0 -AL 0 -aL 1.0 -aS 1.0

  # Remove probes with 99% sequence similarity to get bait design under the 5.9MB limit.
  # Resulted in 49,542 baits which was still 5.94MB
  cd-hit-est -i Final_M_haem_bait_design.fasta -o 99_Final_M_haem_bait_design.fasta -c 0.99

  # Next, cluster by 98% cluster identity.
  # Resulted in 48,394 probes (5.8MB) which also creates space for the GC boosting that Agilent's software will perform.
  cd-hit-est -i Final_M_haem_bait_design.fasta -o 98_Final_M_haem_bait_design.fasta -c 0.98

  #46932
  cd-hit-est -i Final_M_haem_bait_design.fasta -o 95_Final_M_haem_bait_design.fasta -c 0.95

  # 43209
  cd-hit-est -i Final_M_haem_bait_design.fasta -o 90_Final_M_haem_bait_design.fasta -c 0.90

  # 40322
  cd-hit-est -i Final_M_haem_bait_design.fasta -o 85_Final_M_haem_bait_design.fasta -c 0.85

  # 37539
  cd-hit-est -i Final_M_haem_bait_design.fasta -o 80_Final_M_haem_bait_design.fasta -c 0.80

  #
  cd-hit-est -i Final_M_haem_bait_design.fasta -o 75_Final_M_haem_bait_design.fasta -c 0.75

  ```
  * We will go forward with the **final, reduced bait design consisting of 48,394 baits**

  We then tested the **final bait design of 48,394 baits** by mapping it to:
  * The refSeq genome for *M. haemolytica*, "NC_021082.fasta"
    * coverage of whole genome 2,732,423 bases:
      * Mean: 2, Std Dev: 1
      * Minimum: 0, Maximum: 16
  * 13/17 informative loci
    * 100% of informative loci sequence covered by, on average > 5 baits (min = 3, max = 8) per nucleotide.
  * ICE elements
    * On average ICE elements were covered by 1.1, 2.7, and 3.5 baits per nucleotide.

### Updated on 11/2/2020

The bait design with 48,394 baits was determined to be too big and would fit in a large tier with Agilent. CD-HIT is limited to 80% sequence identity, so we tried other methods for reducing the bait design and ultimately decided on using Geneious to map the baits to a reference genome and extracting baits to reduce the overall coverage.

```


Alignment of M_haem_specific_probes_70genomes.fasta to CP006957.fasta with geneious

42,124 baits aligned (4,888 not aligned)

Nucleotide Statistics:
Length: 2,525,685 bp
Sequences: 42,124
Identical Sites: 2,345,213 (98.7%)
Pairwise Identity: 99.1%

Coverage of 2,525,685 bases:
  Mean: 2.0       Std Dev: 0.9
  Minimum: 0    Maximum: 9
  Forward: ?   Reverse: ?
  Ref-Seq: 94.1% (2,375,421 of 2,524,779)

Ungapped lengths of 42,123 reads:
  Mean: 120.0       Std Dev: 0.0
  Minimum: 120    Maximum: 120


    Freq       % of non-gaps
A:  2,222,332  29.3%
C:  1,573,231  20.8%
G:  1,543,169  20.4%
T:  2,240,807  29.6%

GC: 3,116,400  41.1%
All:7,579,539  100.0%
-:  6,445       0.1% (of any)

# Highlighted first 4 baits per position, 41,060 probes extracted
# Highlighted first 3 baits per position, 36,076 probes extracted
```


Combine final components of bait design:
* M. haemolytica specific baits (from kraken)
  * Reduced bait number by setting a 3x coverage threshold for baits aligning to CP006957 and including the baits that didn't map to the reference (n=4888).
* Baits for ICE elements
* Baits to 13 informative loci



```

# concatenate files for temp_final_design
# 44456 final probes
cat all_ICE_and_80per_extra_variants.fasta.probes M_haem_specific_probes_70genomes_bowtie2_to_CP006957.2_UnusedReadsUnpaired.fasta 3X_M_haem_specific_probes_70genomes_bowtie2_to_CP006957.2_extractionUnpaired.fasta probes_to_extracted_regionsUnpairedUnpaired.fasta > Oct4_temp_final_design.fasta

# Remove redundant probes, 44033 final probes
cd-hit-est -i Oct4_temp_final_design.fasta -o nonred_Final_M_haem_bait_design.fasta -c 1

# Cluster by 90% sequence identity, 39478 probes
cd-hit-est -i nonred_Final_M_haem_bait_design.fasta -o Oct4_90_nonred_Final_M_haem_bait_design.fasta -c 0.90

# Cluster by 80% sequence identity, 39478 probes
cd-hit-est -i nonred_Final_M_haem_bait_design.fasta -o Oct4_80_nonred_Final_M_haem_bait_design.fasta -c 0.8

#
## Random sample alignment of baits to M haem reference genome and add that back to the rest of the design.
#
# 3X extracted reads - 35576 baits


cat all_ICE_and_80per_extra_variants.fasta.probes M_haem_specific_probes_70genomes_bowtie2_to_CP006957.2_UnusedReadsUnpaired.fasta 500subtracted_3X_Final_M_haem_bait_design_extractionUnpairedUnpaired.fasta probes_to_extracted_regionsUnpairedUnpaired.fasta > temp_M_haem_reduced_500baits_final_bait_design.fasta

# 43541 clusters
cd-hit-est -i temp_M_haem_reduced_500baits_final_bait_design.fasta -o M_haem_reduced_500baits_final_bait_design.fasta -c 1



#
##
### Bigger bait design
##
#

# 3X extracted reads - 35576 baits

cat M20_E50_M_haem_bait_design.fasta all_ICE_and_80per_extra_variants.fasta.probes probes_to_extracted_regionsUnpairedUnpaired.fasta > temp_Mhaem_sp_full_design.fasta

cd-hit-est -i temp_Mhaem_sp_full_design.fasta -o Mhaem_sp_full_design_M20_E50.fasta -c 1

# Non M haem specific
cat M_haem_70genome_probes_pl120_ps120_e50_m20.fasta all_ICE_and_80per_extra_variants.fasta.probes probes_to_extracted_regionsUnpairedUnpaired.fasta > temp_Mhaem_nonsp_full_design_M20_E50.fasta

cd-hit-est -i temp_Mhaem_nonsp_full_design_M20_E50.fasta -o Mhaem_nonsp_full_design_M20_E50.fasta -c 1

#
##
# Relaxed flags
##
#

M_haem_70genome_probes_pl120_ps120.fasta

cat M_haem_70genome_probes_pl120_ps120.fasta all_ICE_and_80per_extra_variants.fasta.probes probes_to_extracted_regionsUnpairedUnpaired.fasta > temp_Mhaem_nonsp_full_design.fasta

# 104284
cd-hit-est -i temp_Mhaem_nonsp_full_design.fasta -o Mhaem_nonsp_full_design.fasta -c 1

#86781
cd-hit-est -i Mhaem_nonsp_full_design.fasta -o 95_Mhaem_nonsp_full_design.fasta -c .95

# 72240
cd-hit-est -i Mhaem_nonsp_full_design.fasta -o 90_Mhaem_nonsp_full_design.fasta -c .90

# 61168
cd-hit-est -i Mhaem_nonsp_full_design.fasta -o 85_Mhaem_nonsp_full_design.fasta -c .85

# 53563
cd-hit-est -i Mhaem_nonsp_full_design.fasta -o 80_Mhaem_nonsp_full_design.fasta -c .80


```

# Creating a Manifest File for 16S Sequencing Data

This guide outlines a streamlined process for creating a manifest file from 16S sequencing data.

## Prerequisites

Ensure you replace `USERNAME` with your actual username in the paths provided below.

## Steps

### Step 1 & 2: List Absolute Paths

Generate `R1.txt` and `R2.txt` with absolute paths to forward and reverse reads. Change the regular expression to match the suffix for your files (e.g. *R1_001.fastq.gz)

```bash
ls /scratch/users/USERNAME/path/to/reads/*R1.fastq.gz > R1.txt
ls /scratch/users/USERNAME/path/to/reads/*R2.fastq.gz > R2.txt
```

### Step 3: Extract Sample IDs

Extract sample IDs and write them to `sampleIDs.txt`. Remember to replace the "_R1.fastq.gz" to match the suffix in your files that you want to remove for the sampleIDs.

```bash
awk -F'/' '{print $NF}' R1.txt | sed 's/_R1.fastq.gz//' > sampleIDs.txt
```

### Step 4: Combine Files into Manifest

Combine `sampleIDs.txt`, `R1.txt`, and `R2.txt` into `manifest.tsv`.

```bash
paste sampleIDs.txt R1.txt R2.txt > manifest.tsv
```

### Step 5: Add Headers to Manifest

Prepend headers to the manifest file.

```bash
(echo -e "sample-id"$'\t'"forward-absolute-filepath"$'\t'"reverse-absolute-filepath" && cat manifest.tsv) > manifest_with_header.tsv && mv manifest_with_header.tsv manifest.tsv
```

## Conclusion

Following these steps, you'll have a ready-to-use manifest file, crucial for bioinformatics pipelines involving sequencing data.

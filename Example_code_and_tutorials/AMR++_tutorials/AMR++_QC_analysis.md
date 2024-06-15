# QC analysis

In AMR++, quality control (QC) analysis of metagenomic sequencing data occurs in two subworkflows or pipelines. The pipeline is used "eval_qc" to make QC reports and select trimming parameters, which is then paired with "trim_qc" to run Trimmomatic and output QC trimmed reads.


# Contents
* [##first ]

# First step - Evaluate QC

## Running the "eval_qc" pipeline

Input to the "eval_qc" pipeline:
* "--reads" flag pointing to paired end reads
* "--output" flag to store results

The "eval_qc" pipeline runs two tools, fastqc and multiqc, to help inform two major questions:
1. Was the sequencing performance adequate?
2. What trimming parameters should I choose for QC trimming?

FastQC analyzes your sample reads and multiQC summarizes the data in a useful html report.

An example command to run the pipeline looks like this:
```
# From inside the AMRplusplus repository on your server
nextflow run main_AMR++.nf -profile local --pipeline eval_qc --reads "/PATH/TO/YOUR/READS/_R{1,2}_001.fastq.gz" --output AMR++_results
```

You can find more information about how to make the regular expression matching your read location in the [Frequently Asked Questions page].

## eval_qc output

Now, let's look at the output folder "AMR++_results" new directories that were created inside, "AMR++_results/


## Example reports
For an example of some reports, check out these links:
* Example of good QC scores
* Example of bad QC scores
* Example of weird QC scores

# Second step - Trim reads
Next, we will use Trimmomatic help improve the overall QC scores and accuracy of your downstream analysis. Based on the QC reports for your samples, you'll select adequate trimming parameters and also be able to identify any failures in sequencing for specific samples.

## Choose trimming parameters

The parameters we have to change are:
```
--leading = 3
--trailing = 3
--slidingwindow = "4:15"
--minlen = 36
```

Remember, we can change these parameters by adding them to the command itself, or by modifying the "params.config" file in the main directory. 

## Running the "trim_qc" pipeline

To run this next pipeline, we have to 

Notice, we'll still be pointing at the same raw reads as above so the "--reads" flag won't change. 

An example command to run the pipeline looks like this:
```
# From inside the AMRplusplus repository on your server
nextflow run main_AMR++.nf -profile local --pipeline trim_qc --reads "/PATH/TO/YOUR/READS/_R{1,2}_001.fastq.gz" --output AMR++_results
```


# Next step, remove host DNA from your QC trimmed reads.
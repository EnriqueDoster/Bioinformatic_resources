# Running AMR++ in pieces

## Installing AMR++

* Anaconda
* Singularity + nextflow
* Install all tools on your server (manually or with modules)

## Understanding how to change how AMR++ runs
If you run AMR++ using the default settings like this:

```
nextflow run AMR++_main.nf 
```

AMR++ knows what pipeline to run and with which samples by looking at the file "params.config" and using the default settings for those variables. We can change those values directly in the params.config file, or add any of those variables to the command itself to override what's in the 


## Understanding the use of "-profile"

AMR++ uses many different types of software and custom scripts to analyze your data. For AMR++ to know where to find these tools, we can specify different options to the "-profile" flag. Notice that only 1 dash "-", this detail is significant to the underlying working of the pipeline, but for now just keep track of whether each flag has 1 or 2 dashes.

By default, AMR++ assumes all software dependencies are in your \$PATH, as in the "local" profile. Here are the other options:

        - local: Assumes all sofware is already in your \$PATH
        - local_slurm: Local installation and adds control over slurm job submission.
        - conda: Uses "conda" to install a conda environment. 
        - conda_slurm: Uses "conda" and adds control over slurm job submission.
        - singularity: Uses a "singularity" image container.
        - singularity_slurm: Singularity image and adds control over slurm job submission.
        - docker: Uses a docker image container.


## AMR++ and the --pipeline flag

AMR++ can run the entire pipeline including QC evaluation, QC trimming, host removal, resistome and microbiome analysis. However, running the entire pipeline on large datasets can be computationally prohibitive. The main reason for this is that nextflow creates many temporary files while running the pipeline. The temporary files can then be used to "resume" an ongoing pipeline if it's interrupted for whatever reason. Unfortunately, these files can often take 1-2X the amount of storage as your input files. Therefore, we provide the option to run each component of the AMR++ pipeline individually by changing the "--pipeline" flag in your command. 


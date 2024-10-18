# Host removal


# The --rm_host ##

## Host removal parameters

The parameters we have to change are:

```
--reads = "/path/to/qc_trimmed/reads/*R{1,2}.fastq.gz"
--host = "${baseDir}/data/host/chr21.fasta.gz"
--host_index = "${baseDir}/data/host/chr21.fasta.gz*"
```

If you don't have a host_index yet, you can leave that 
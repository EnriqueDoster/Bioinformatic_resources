# TAMU server use SOP

# Table of contents
* [Computing servers](#computing-servers)
* [Example AMR++ run](#example-amr-run)


# Computing servers
* Terra
  * newest, 304 computing nodes, some GPUs
  * still being developed and new tools are being added
* Ada
  * main computing nodes 
  * 15 nodes are 1TB and 2TB
  * loaded with many bioinformatic tools
* Curie
  * 50 nodes (256GB each)
  * loaded with many bioinformatic tools
  * Curie is being retired August 31, 2020
* Lonestar – joint effort with UT Austin
  * computing hours available


# Example AMR++ run
1. Sign-in to a server, ada.
  * Log on to VPN (Anyconnect) if you are not on TAMU’s network, connect.tamu.edu
2. Enter username and password
```
ssh enriquedoster@ada.tamu.edu
```
3. Navigate the server directories
/scratch/user/$USER

4. Prepare directory for running AMR++
module load Singularity/2.5.2-GCC-6.4.0-2.28

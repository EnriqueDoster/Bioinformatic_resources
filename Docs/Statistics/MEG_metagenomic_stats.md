# MEG statistics overview
## **partial document in progress**


# Table of contents
* [Normalization](#normalization)
* [Characterization](#characterization)
* [Alpha diversity](#alpha-diversity)


# Normalization
Often needed to account for differences in sequencing depth. What method you use depends on the analysis you want to run. Not all analyses require normalization.
* Rarefying - not used by our group, but commonly accepted. 
* Total sum scaling - basic way of normalizing based on total counts in each sample.
* Cumulative sum scaling - currantly our main method for normalization


# Characterization
Performed on normalized counts.

* Relative abundance plots
* Upset plots

# Alpha diversity
Performed on raw counts at the lowest taxonomic level (ASVs, gene accessions).
* Richness (observed features)
* Eveness 
  * Typically we use "Shannon's index"
  * Other options are "Simpson's index" and "Inverse Simpson", 


# Beta diversity


# Differential abundance testing


# Correlation/
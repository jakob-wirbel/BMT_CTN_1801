# BMT_CTN_1801 Microbiome

Analysis code for the BMT CTN 1801 microbiome project.

Please refer to our publication for more information:
> [Authors. Title. bioRxiv (2025)](https://www.biorxiv.org/)

## Repository organization

### Data

Some of the data for the presented analyses are available via 
[Zenodo](https://doi.org/10.5281/zenodo.17127362). 
Please refer to the BMT CTN Data Access Committee for information
on how the metadata for all samples and participants can be
obtained.  
Unfortunately, this repository is not a full plug-and-play
reproduction of results, since the metadata and participant
outcome data are still under embargo and can therefore not
be provided with other data on Zenodo.

### Raw sequencing data

The raw sequencing data for all metagenomic stool samples in the project
have been uploaded to ENA and are available under the identifier
[PRJEB97686](https://www.ebi.ac.uk/ena/browser/view/PRJEB97686). 
The derived data tables (sequencing statistics, motus and metaphlan
profiles) are available on Zenodo as well.

### Scripts

Here is a quick description of the included scripts, in the order that they
were executed to produce the results presented in the paper:

1. `absolute_abundance_prediction`  
This script uses the relationship between DNA concentration and 16S copy number
to train a machine learning model for the prediction of 16S copy number. Please
see the [ancillary publication](https://doi.org/10.1016/j.crmeth.2025.101030) 
for more detail. A second batch of samples were measured by ddPCR, as we had 
selected a subset of samples with relatively high concentration in the 
first batch.

	input:
	- meta_samples.tsv
	- sequencing_stats.tsv
	- ddPCR_results_*
	- ddPCR_layout_*
	- metaphlan_all.tsv

	output:
	- copies_16S.tsv
	- figures/ddPCR/*

2. `clean_tables`  
This script collates all data and metadata and saves the resulting tables in
an `.RData` object for easier use in subsequent scripts. Additionally, some
quality checks and comparisons (between `next_day` and `same_day` sample 
processing) are performed.

input:
- meta_samples.tsv
- copies_16S.tsv
- sequencing_stats.tsv
- motus_all.tsv
- metaphlan_all.tsv
- meta_participants.tsv

output:
- figures/sample_type/*
- figures/general/*
- all_data.RData
   - feat.motus     # motus table
   - feat.metaphlan # metaphlan table
   - df.meta.clean  # cleaned metadata table for each sample
   - df.response    # cleaned metadata table for each participant

3. `diversity`  
This script explores the trajectory of microbial alpha diversity and absolute
abundance over time in both treatment arms.

input:
- all_data.RData

output:
- figures/alpha/*


4. `abx_analyses`  
This script  

5. `domination`  
This script 

6. `associations`  
This script 

7. `humann`  
This script 

#### Primary and secondary clinical trial endpoints

The primary and secondary clinical trial analyses were performed by Michael J.
Martens at the Division of Biostatistics, Medical College of Wisconsin, WI, USA,
based on the pre-specified statistical plan. The official BMT CTN 1801 analysis
report is confidential/available upon request (? actually, not sure...
will have to ask Mike ?).

## Contact

If you have any questions/issues, please feel free to open an 
issue in this repo or reach out to
_wirbel[at]stanford.edu_.

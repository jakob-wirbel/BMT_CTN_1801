# BMT_CTN_1801 Microbiome

Analysis code for the BMT CTN 1801 microbiome project.

Please refer to our publication for more information:
> [Wirbel, Saber, Martens et al. Differential effects of two common GVHD prophylaxis regimens on the gut microbiome: Results from the BMT CTN 1801 study. bioRxiv (2026)](https://www.biorxiv.org/content/10.64898/2026.02.19.706769v1)

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
profiles) are available on Zenodo.

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
This script investigates the influence of different antibiotics on the absolute
microbial load. Treatment arms are compared by their exposure to antibiotics
and each sample is classified as being exposed to antibiotics or not. Additionally,
the influence of prophylactic vs non-prophylactic antibiotics is explored, 
showing that the difference between GVHD prophylaxis arms in absolute microbial
abundance is not due to differences in antibiotics exposure.

	input:
	- all_data.Rdata
	- antibiotics_data.tsv

	output:
	- figures/abx/*
	

5. `domination`  
This script deals with domination events, when a single microbial species is
present at more than 30% relative abundance in a sample. Differences between 
study arms are explored and what type of species lead to domination.
Interestingly, not all commonly dominating species are pathogens, but also 
species with comparably high relative abundance in not-dominated samples, 
arguing for a nuanced view on microbial domination. Associations with outcome
are also explored.

	input:
	- all_data.Rdata

	output:
	- figures/domination/*

6. `composition`  
This script investigates differences in microbiome composition between treatment
arms. Overall relative abundance differences are computed and additional 
analyses zoom in on _Clostridium scindes_ as top hit in the differential 
abundance testing.

	input:
	- all_data.Rdata
	- files/abx_exposure.RData
	- files/mOTUs_3.0.0_GTDB_tax.tsv

	output:
	- figures/composition/*
	- files/enrichment_results.tsv
	
7. `humann`  
This script explores the functional microbiome profiles generated with the
Humann package, mapped to metacyc pathways. Overall differences between 
treatment arms are explored and, similar to compositional differences,
additional analyses zoom in on bile acid modification pathways as the top hit
in the differential abundance analysis.

	input:
	- all_data.Rdata
	- data/humann.tsv

	output:
	- figures/humann/*

## Primary and secondary clinical trial endpoints

The primary and secondary clinical trial analyses were performed by Michael J.
Martens at the Division of Biostatistics, Medical College of Wisconsin, WI, USA,
based on the pre-specified statistical plan. The official BMT CTN 1801 analysis
results are included as Supplementary Information in the manuscript (see
link above).

## Contact

If you have any questions/issues, please feel free to open an 
issue in this repo or reach out to
_jakob.wirbel[at]helmholtz-hzi.de.

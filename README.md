# S-EPS
R code of a novel Site-based mutation dynamics framework to profile the source-sink dynamics of key amino acid mutations in major SARS-CoV-2 variants. The framework is demonstrated using a randomly sampled dataset from the period when the Delta Variant of Concern (VOC) became prevalent (B = 1, June 2020 < t < August 2021, mutation k ∈ Delta VOC characteristic mutations as described in the Methods section).

Input data:

First Stage Sampling: Genetic sequences of the Spike gene of SARS-CoV-2 collected from thirteen geographic regions between June 2020 and August 2021. The sample size for each month-region stratum is up to 1,705 sequences. (Folder: first_stage_sampling)
Second Stage Sampling: Genetic sequences of the Spike gene of SARS-CoV-2 collected from thirteen geographic regions between June 2020 and August 2021. The sample size for each month-region stratum is up to 163 sequences. (Folder: second_stage_sampling)
Reference sequence: The Spike gene of the original SARS-CoV-2 strain, hCoV-19/Wuhan/Hu-1/2019 (File: wuhan_reference_spike_aa.csv)
Sample Size Data: The size of sequences for each month-region stratum collected from GISAID. (File: sample_size.csv)
List of Delta-characteristic mutations on the Spike gene (File: delta_mutations.csv)

Installation and setup:

Step 0: Install dependent packages R (4.1.3) or newer version is required. To download R, please see the website https://www.r-project.org/ for detailed instructions. The R packages prettyR (2.2-3) and scales (1.2.1) are pre-requisites, which can be installed from the following websites.

prettyR (2.2-3): https://www.rdocumentation.org/packages/prettyR/versions/2.2-3
scales (1.2.1): https://www.rdocumentation.org/packages/scales/versions/1.2.1

Analysis Pipeline:

Step 1: Identify Delta-characteristic key mutations.
 
Step 1.1: Fitness criteria.
Calculate the proportion of mutations at each amino acid site across the Spike gene, and identify all amino acid mutations with a proportion of no less than 50% in at least one month-region stratum. Output the list of mutations that meet the Fitness criteria of key mutations as described in the Methods section. Run time: ~1 min.

Step 1.2: Uniqueness criteria 
Filter the potential key mutations identified in Step 1.1 using the Uniqueness criteria described in the Methods section to identify Delta-characteristic mutations. Output the list of mutations that meet the Uniqueness and Fitnes criteria of key mutations as described in the Methods section

Step 1.3: Coverage criteria
Among the potential key mutations from Step 1.2, identify mutations that have been observed with a prevalence exceeding 20% in more than half of the 13 regions. Output the list of key mutations that meet all the key mutation criteria as described in the Methods section. Run time: ~1 min.

Step 2: Determine occurrence ranking using Equal Power Sampling (EPS) Framework
Run time: ~1 min

Step 2.1: Identification of source region
Use the first stage sampling data to determine the source region of each Delta-characteristic key mutation.

Step 2.2: Occurrence ranking of the remaining regions
Use the second stage sampling data to determine the occurrence ranking of the remaining regions for each Delta-characteristic key mutation.

Step 2.3: Normalization
Normalize the ranking to a scale of 1-13 to facilitate joint analysis of all mutations. Output the normalized detection rank of each region along the transmission pathway of each key mutation

Step 3: Calculate rank statistics probability
Use rank statistics to infer the relative positions of geographic regions along the viral transmission pathway. Output the matrix of probability of each region being the rth order along the key mutation transmission among the selected regions.

Data Attribution

The genetic sequences provided in the folders first_stage_sampling and second_stage_sampling include the virus names of SARS-CoV-2 sequences used in this study. These viral sequence data were downloaded from the Global Initiative on Sharing All Influenza Data (GISAID) at http://platform.gisaid.org/. We thank the contributions of all the health care workers and scientists, the GISAID team, and the submitting and the originating laboratories.

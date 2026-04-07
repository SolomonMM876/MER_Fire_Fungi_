HMSC Fire Response Analysis
==============================================================================

OVERVIEW
------------------------------------------------------------------------------
This repository contains the data and scripts used to analyze the response of 
mycorrhizal fungal communities to fire regimes in Australian systems. The 
analysis involves:
1.  Calculating Alpha Diversity and testing fire effects using LMMs.
2.  Using Hierarchical Modelling of Species Communities (HMSC) to estimate 
    species-level responses to fire (beta_fire).
3.  Analyzing how beta_fire varies by site, fire traits, and fungal functional 
    traits.


DIRECTORY STRUCTURE
------------------------------------------------------------------------------
1. data/
   Contains all raw and clean CSV inputs required to run the analysis.

2. scripts/
   Contains the R analysis scripts.

3. plots/
   Contains all plots generated 
   
4. processed_data/
   Contains summary outputs of all analyses
   
5. HMSC_Output/
   Contains output of HMSC analysis in script 1


DATA FILES (data/)
------------------------------------------------------------------------------
1. Input Data (Used for HMSC Fitting & Alpha Diversity)
   
   a. input_site_metadata.csv
      Contains geospatial and treatment information for every plot.
      - Site, Plot, Lat, Long, Fire_Treatment.

   b. input_soil_community.csv
      The raw DNA sequence count matrix (OTU table) for soil samples.
      - Plot: Unique identifier matching the metadata.
      - ITSall_...: Sequence counts for fungal OTUs.

   c. input_host_frequency.csv
      Vegetation covariate data.
      - freq_AM, freq_ECM: Frequency of mycorrhizal host plants.


2. Analysis Data (Used for Beta Coefficient Analysis)

   a. beta_estimates.csv
      Posterior mean estimates of response to fire from HMSC model.

   b. site_metadata_analysis.csv
      Detailed site metadata (Aridity, TSF, Vegetation Type).

   c. taxonomy.csv
      Consensus taxonomy and foraging strategies for fungal OTUs.

   d. soil_abundance.csv & hyph_abundance.csv
      Sequence count data used for Relative Abundance calculations.

   e. traits_repo_strat.csv
      Fungal reproductive strategies (Epigeous, Hypogeous, etc.).

   f. traits_amf_guilds.csv
      AMF functional guild classifications (Weber et al. 2019).


ANALYSIS WORKFLOW
------------------------------------------------------------------------------
1. Run 'Script_01_HMSC_Fitting.R'
   Fits the HMSC model using the inputs in 'data/' and saves the model object.

2. Run 'Script_02_Beta_Analysis.R'
   Reads CSVs from 'data/', performs linear models and ANOVAs on the beta 
   coefficients, and generates publication plots for Beta responses.

3. Run 'Script_03_Alpha_Diversity.R'
   Reads community data from 'data/', calculates diversity metrics (Shannon, 
   Simpson, Chao1), and performs LMMs to test for fire treatment effects.


OUTPUTS (project_script/)
------------------------------------------------------------------------------
* Plots (project_script/plots/):
    - Effect_Size_site_Plot.png
    - Effect_Size_Fire_interval.png
    - Effect_Size_aridity.png
    - Effect_Size_repo_strat_Plot.png
    - Effect_Size_Guild_Plot.png
    - Effect_Size_explo_Plot.png
    - Effect_Size_phylum_Plot.png
    - RA_OTU_Effect_Size_Plot.png

* Data (project_script/processed_data/):
    - Alpha_Diversity_ANOVA.csv: ANOVA results for alpha diversity metrics.
    - Master_ANOVA_Results.csv: Combined ANOVA results for Beta analyses.
    - Master_PostHoc_Results.csv: Combined pairwise comparison results.
##################################################

OUTPUT FILE: ITS_output_clean.tsv
Tab-delimited text file containing OTU counts per sample and OTU metadata. 

VARIABLE NAME	---	DESCRIPTION
barcode	---	sample index used to process DNA sequencing data
sample	--- lab/field sample code (if available)
OTU --- operational taxonomic unit identifier
count	---	number of reads observed
resampled_count	---	number of reads following random resampling to constant depth
resampled_date	---	date when resampling performed (for tracking versions)
kingdom:species	---	taxonomic assignment (if NA, no close match to BLAST database)
n	---	number of close BLAST database matches 
score	--- proportion of close BLAST database matches belonging to the assigned taxon (assignment made if one taxon makes up more than 50% of close matches)
starts with 'SH'	---	additional metadata from the top match in UNITE
identity	---	proportion of bases matching in alignment between query sequence and top hit
coverage	---	proportion overlap between query sequence and top hit
trophicMode:citationSource	---	FUNGuild assignment (if NA, no match)
unite_db_version	---	UNITE database version used for taxonomic assignments (for tracking versions)
date_assigned	---	date when taxonomic assignments performed (for tracking versions)

##################################################

OUTPUT FILE: ITS_output_summarised.tsv
Tab-delimited text file containing community structure estimates for each sample.

VARIABLE NAME	---	DESCRIPTION
barcode	---	sample index used to process DNA sequencing data
sample	--- lab/field sample code (if available)
n_reads	---	number of DNA sequence reads assigned to Fungi at the phylum level or lower
n_reads_resampled	---	number of DNA sequence reads used to calculate estimates requiring resampling
obs_richness_fungi	---	observed number of fungal OTUs following random resampling to constant sequencing depth per sample
chao1_richness_fungi	---	estimated number of fungal OTUs using the Chao1 estimator and OTU counts before resampling
obs_shannon_fungi	---	observed Shannon diversity of fungal OTUs observed following random resampling to constant sequencing depth per sample
obs_richness_fungi_ArbuscularMycorrhizal	---	observed number of OTUs assigned to arbuscular mycorrhizal fungal taxa following random resampling to constant sequencing depth per sample
obs_richness_fungi_Ectomycorrhizal	---	observed number of OTUs assigned to ectomycorrhizal fungal taxa following random resampling to constant sequencing depth per sample
obs_richness_fungi_OrchidMycorrhizal	---	observed number of OTUs assigned to orchid mycorrhizal fungal taxa following random resampling to constant sequencing depth per sample
obs_richness_fungi_PlantPathogen	---	observed number of OTUs assigned to plant pathogenic fungal taxa following random resampling to constant sequencing depth per sample
obs_richness_fungi_WoodSaprotroph	---	observed number of OTUs assigned to wood saprotrophic fungal taxa following random resampling to constant sequencing depth per sample
chao1_richness_fungi_ArbuscularMycorrhizal	---	estimated number of OTUs assigned to arbuscular mycorrhizal fungal taxa using the Chao1 estimator and OTU counts before resampling
chao1_richness_fungi_Ectomycorrhizal	---	estimated number of OTUs assigned to ectomycorrhizal fungal taxa using the Chao1 estimator and OTU counts before resampling
chao1_richness_fungi_OrchidMycorrhizal	---	estimated number of OTUs assigned to orchid mycorrhizal fungal taxa using the Chao1 estimator and OTU counts before resampling
chao1_richness_fungi_PlantPathogen	---	estimated number of OTUs assigned to plant pathogenic fungal taxa using the Chao1 estimator and OTU counts before resampling
chao1_richness_fungi_WoodSaprotroph	---	estimated number of OTUs assigned to wood saprotrophic fungal taxa using the Chao1 estimator and OTU counts before resampling
obs_shannon_fungi_ArbuscularMycorrhizal	---	observed Shannon diversity of OTUs assigned to arbuscular mycorrhizal fungal taxa following random resampling to constant sequencing depth per sample
obs_shannon_fungi_Ectomycorrhizal	---	observed Shannon diversity of OTUs assigned to ectomycorrhizal fungal taxa following random resampling to constant sequencing depth per sample
obs_shannon_fungi_OrchidMycorrhizal	---	observed Shannon diversity of OTUs assigned to orchid mycorrhizal fungal taxa following random resampling to constant sequencing depth per sample
obs_shannon_fungi_PlantPathogen	---	observed Shannon diversity of OTUs assigned to plant pathogenic fungal taxa following random resampling to constant sequencing depth per sample
obs_shannon_fungi_WoodSaprotroph	---	observed Shannon diversity of OTUs assigned to wood saprotrophic fungal taxa following random resampling to constant sequencing depth per sample
relabund_fungi_ArbuscularMycorrhizal	---	proportion of ITS2 sequence reads ('reads.fungi') assigned to arbuscular mycorrhizal fungal taxa
relabund_fungi_Ectomycorrhizal	---	proportion of ITS2 sequence reads ('reads.fungi') assigned to ectomycorrhizal fungal taxa
relabund_fungi_OrchidMycorrhizal	---	proportion of ITS2 sequence reads ('reads.fungi') assigned to orchid mycorrhizal fungal taxa
relabund_fungi_PlantPathogen	---	proportion of ITS2 sequence reads ('reads.fungi') assigned to plant pathogenic fungal taxa
relabund_fungi_WoodSaprotroph	---	proportion of ITS2 sequence reads ('reads.fungi') assigned to wood saprotrophic fungal taxa

##################################################

OUTPUT FILE: ITS_data_quality_summary.tsv
Tab-delimited text file containing counts of the number of DNA sequence reads remaining after each filtering step, plus other sample metadata.

VARIABLE NAME	---	DESCRIPTION
barcode	---	sample index used to process DNA sequencing data
sample	--- lab/field sample code (if available)
starts with 'dna'	---	estimates of DNA concentration and quality before normalising
n_R1	---	number of DNA sequence reads in the R1 file
n_R2	--- number of DNA sequence reads in the R2 file
n_merged	--- number of DNA sequence reads following successful merging of R1 and R2 reads
n_passed	--- number of DNA sequence reads passing the quality threshold
n_assigned_to_fungi	--- number of DNA sequence reads assigned to kingdom Fungi
n_assigned_to_fungi_and_phylum	--- number of DNA sequence reads assigned to a phylum within kingdom Fungi
n_after_removing_rare_otus	--- number of DNA sequence reads assigned to a phylum within kingdom Fungi and remaining after removing questionable observations (n < 10 in a sample)

##################################################

OUTPUT FILE: ITS_output_repseqs.fasta
Fasta file containing the representative sequence for each detected OTU

##################################################

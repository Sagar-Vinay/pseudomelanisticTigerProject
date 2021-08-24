# pseudomelanisticTigerProject
This repository contains scripts for the sequence data analysis, poulation genetic analysis and simulations done for the psueomelanistic tiger genetics project.
Publication - High frequency of an otherwise rare phenotype in a small and isolate tiger population, Sagar et al., 2021.

List of files in this repository:

1. genomeAnalysisPipeline.sh - Script with commands to analyse raw fastq sequence data and obtain the variant calling format (VCF) file for Taqpep genomic region.
2. singleReadmpcrDataAnalysisPipeline.sh - Script with commands to process raw single-read fastq files obtained from mPCR-MiSeq/mPCR-HiSeq experiment and obtain VCF file.
3. pairedReadmpcrDataAnalysisPipeline.sh - Script with commands to process raw paired-read fastq files obtained from mPCR-MiSeq/mPCR-HiSeq experiment and obtain VCF file.
4. mpcrVCFfiltering.sh - Script with commands to filter the VCF file obtained from mPCR-MiSeq/mPCR-HiSeq data and obtain pairwise individual relatedness values.
5. driftProbPseudomelanismFreq0.5founder1het.R - R script for isolated population, discrete generation - post bottleneck drift simulaitons to estimate probability of pseudomelanistic allele frequency reaching 0.5 or above.
6. timeToFixation6popModels.R - R script for estimating the time to fixation of one allele over another in different population growth models with or without migration.
7. taqpepGenotypedIndividuals.csv - A CSV file detailing the 428 tigers genotyped at Taqpep p.H454Y site with biosample and bioproject accession wherever applicable. For some bengal tigers only taqpep genomic region bam file is uploaded to SRA database (mentioned as 'TaqpepRegionBamFile' in the 'Data Type' coloumn) becuase this data was generated for other projects that are still undergoing. 
8. samplesUsedForPopGenAnalysis.csv - A CSV file detailing the 137 tigers used for different population genetics analyses in the paper with Biosample and Bioproject accesions.  

# In case of querries/concerns, email  - vinays@ncbs.res.in

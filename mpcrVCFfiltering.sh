#Additional supplemental material for -
#High frequence of an otherwise rare phenotype in a small and isolated tiger population#
#Goes with Methods section - Species identification of the noninvasive samples and individual recaptures
#Commands used for filtering the VCF file prepared by processing the raw fastq files obtained after the multiplex PCR and MiSeq sequencing; script modified from Natesh et al, 2019 - Empowering conservation practice with efficient and economical genotyping from poor quality samples.

#This script includes steps of filtering the mpcr data VCF before individual identification is done using pi-hat relatedness values.

#Remove indels and uncompress the vcf - requires vcftools
vcftools --gzvcf mpcr_bcftools_consensusCaller.vcf.gz --remove-indels --out mpcr_bcftools_consensusCaller_noIndel --recode #assumes the input file name as mpcr_bcftools_consensusCaller

#Filtering out genotypes based on genotype quality and depth - requires GATK
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R referenceGenome.fa --setFilteredGtToNocall --genotypeFilterExpression "GQ < 10" --genotypeFilterName "GQbelow10" --variant mpcr_bcftools_consensusCaller_noIndel.recode.vcf -o mpcr_bcftools_consensusCaller_noIndel_GQ10.vcf

java -jar GenomeAnalysisTK.jar -T VariantFiltration -R referenceGenome.fa --setFilteredGtToNocall --genotypeFilterExpression "DP < 10" --genotypeFilterName "DPbelow10" --variant mpcr_bcftools_consensusCaller_noIndel_GQ10.vcf -o mpcr_bcftools_consensusCaller_noIndel_GQ10_DP10.vcf

#Calculate missing data  per site - requires vcftools
vcftools --vcf mpcr_bcftools_consensusCaller_noIndel_GQ10_DP10.vcf --missing-site --out mpcr_bcftools_consensusCaller_noIndel_GQ10_DP10

#following this, loci with missing data for more than 10% individuals and minor allele count 1 were filtered out - requires vcftools
vcftools --vcf mpcr_bcftools_consensusCaller_noIndel_GQ10_DP10.vcf --max-missing 0.9 --mac 1 --out mpcr_bcftools_consensusCaller_noIndel_GQ10_DP10_mac1_maxMiss0.9 --recode

#caldulate missing data per sample - requires vcftools
vcftools --vcf mpcr_bcftools_consensusCaller_noIndel_GQ10_DP10_mac1_maxMiss0.9.recode.vcf --missing-indv --out mpcr_bcftools_consensusCaller_noIndel_GQ10_DP10_mac1_maxMiss0.9

#Filtered out samples with genotypes called at less than 50 sites (out of remaining 102) - requires vcftools
awk '$4>52' mpcr_bcftools_consensusCaller_noIndel_GQ10_DP10_mac1_maxMiss0.9.imiss | cut -f1 > indWithMissingDataOn52sites.txt #from 102 sites, samples with genotypes called at minimum 50 sites were kept
vcftools --vcf mpcr_bcftools_consensusCaller_noIndel_GQ10_DP10_mac1_maxMiss0.9.recode.vcf --remove indWithMissingDataOn52sites.txt --out mpcr_bcftools_consensusCaller_noIndel_GQ10_DP10_mac1_maxMiss0.9_ind50snps --recode

#Estimating pi-hat relatedness - requires plink
plink-1.9/plink --vcf mpcr_bcftools_consensusCaller_noIndel_GQ10_DP10_mac1_maxMiss0.9_ind50snps.recode.vcf --genome --aec --const-fid --out mpcr_bcftools_consensusCaller_noIndel_GQ10_DP10_mac1_maxMiss0.9_ind50snps

#From the .genome file obtained the pi-hat values for each pair of samples and identified the individuals as described in the methods section - Species identification of the noninvasive samples and individual recaptures.

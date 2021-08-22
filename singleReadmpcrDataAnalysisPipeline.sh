#Additional supplemental material for -
#High frequence of an otherwise rare phenotype in a small and isolated tiger population#
#Goes with Methods section - Species identification of the noninvasive samples and individual recaptures
#Commands used for obtaining and filtering SNPs from the raw fastq files obtained after the multiplex PCR and MiSeq sequencing; script modified from Natesh et al, 2019 - Empowering conservation practice with efficient and economical genotyping from poor quality samples.
#This script is a pipeline to obtain a variant calling format file from raw fastq files (single read)

#Make a separate directory for each analysis run, and within this directory create subdirectories for raw data, trimmed files, bam files etc#
##Move into this directory and define path to the analysis directory
DATA_PATH='/home/uramakri/vinays/trialRun'

#copy all the raw fastq files in $DATA_PATH}/rawData

#Making subdirectories within this directory
mkdir trimmedFiles ##trimmed files to be stored here##
mkdir bamFiles ##bamFiles to be stored here##
mkdir vcfFiltering ##vcf file to be created and filtered here##

#Trimming the fastq files - requires trim_galore

for file in ${DATA_PATH}/rawData/*.gz
do
trim_galore --quality 30 --phred33 --fastqc --illumina --stringency 5 --length 5 --output_dir ${DATA_PATH}/trimmedFiles ${file}
done

#Mapping the trimmed reads to reference genome - requires bwa, samtools

for file in ${DATA_PATH}/trimmedFiles/*.fq.gz
do
base=$(basename ${file} _L001_R1_001_trimmed.fq.gz)#assumes trimmed files have extension '_L001_R1_001_trimmed.fq.gz'
readGroupString=@RG\\tID:${base}\\tSM:${base}\\tLB:${base}\\tPL:illumina
bwa-0.7.17/bwa mem -M -R ${readGroupString} -t 30 -B 3 referenceGenome.fa ${file} | gzip -3 | samtools view -hu - | samtools sort -O bam - > ${DATA_PATH}/bamFiles/${base}_bwa_mem.bam
samtools index ${DATA_PATH}/bamFiles/${base}_bwa_mem.bam
done

#Variant calling at specific positions - requires bcftools
#first need to create list of bam files
cd ~/
ls ${DATA_PATH}/bamFiles/*.bam > ${DATA_PATH}/bamFiles/listBams #this will create the file name list along with the path in bamFiles directory#
cd ${DATA_PATH}/
mkdir ${DATA_PATH}/temp #create a temporary directory for bcftools to store temperary files#
#now beginning the variant calling at specific postions (SNP panel)#
bcftools-1.6/bcftools mpileup -Ou --max-depth 10000 -f referenceGenome.fa --bam-list ${DATA_PATH}/bamFiles/listBams -R SNP-positions-file --annotate FORMAT/AD,FORMAT/DP,INFO/AD | bcftools-1.6/bcftools call -Ou -c -f GQ | bcftools-1.6/bcftools sort --temp-dir ${DATA_PATH}/temp -Oz -o ${DATA_PATH}/vcfFiltering/output-file.vcf.gz

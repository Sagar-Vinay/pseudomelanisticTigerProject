#Additional supplemental material for -
#High frequence of an otherwise rare phenotype in a small and isolated tiger population#
#Goes with Methods section - identifying the causal mutation
#Commands used for obtaining variants in the taqpep genomic region from raw fastq files (whole genome sequences of captive tigers)

#Trimming the raw sequence reads, executed for each sample - requires trim_galore
trim_galore --quality 30 --phred33 --fastqc --paired --stringency 3 --length 30 --max_n 1 --retain_unpaired --output_dir path/to/output/directory sample_R1.fastq.gz sample_R2.fastq.gz

#Mapping trimmed reads to reference genome, executed on each file - requires bwa, samtools
gunzip <sample_R1_val_1.fq.gz> sample_R1_val_1.fq #assumes the trimmed file has extension _R1_val_1.fq.gz
gunzip <sample_R2_val_2.fq.gz> sample_R2_val_2.fq

bwa-0.7.17/bwa mem -M -t30 referenceGenome.fa sample_R1_val_1.fq sample_R2_val_2.fq | samtools-1.9/samtools view -hu - | samtools-1.9/samtools sort -O bam - > sample_bwa.mem.bam

samtools-1.9/samtools index inputFile_bwa.mem.bam

#merging bam files of the same sample sequence in separate lanes, executed for samples with data from multiple lanes - requires samtools
samtools-1.9/samtools merge sample_merged_bwa.mem.bam sample_a_bwa.mem.bam sample_b_bwa.mem.bam #here sample_a and sample_b are two bams for the same sample sequenced on different lanes or in different sequencing run and sample_merged_bwa.mem.bam is the merged bam file

#Marking duplicate reads, executed on each bam file - requires picardtools
java -jar picard.jar MarkDuplicates I=sample_bwa.mem.bam O=sample_mdup.bam M=sample_marked_dup_matrix.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

#Adding read groups to bam files - requires picardtools
for file in /path/to/bam/files/*.bam
do
base=$(basename ${file} _mdup.bam)#change extension according to input file name#
echo ${base}
java -jar picard.jar AddOrReplaceReadGroups INPUT=${file} OUTPUT=output/directory/${base}_mdupRG.bam RGID=${base} RGLB=library1 RGPL=illumina RGPU=unit1 RGSM=${base}
done

#testing the resultant bam file for integrity and truncation, executed on each bam file - requires picardtools
java -jar picard.jar ValidateSamFile I=sample_mdupRG.bam MODE=SUMMARY MAX_OPEN_TEMP_FILES=1000

#subsampling the data for Taqpep genomic region, executed on each bam file - requires samtools
samtools-1.9/samtools view -uh sample_mdupRG.bam chrA1:94,850,000-94,950,000 > sample_taqpep.bam #the coordinates are for 100kb region containing Taqpep gDNA sequence in reference genome felCat8.0

#Calling variants - requires freebayes
freebayes -f referenceGenome.fa -L list/of/bam/files -F 0.2 -C 1 --haplotype-length 0 > output.vcf

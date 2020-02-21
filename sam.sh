#!/bin/bash

#Variables for the reference fasta, bam files, samtools software,  
samtools="/Users/chin-hua-huang/samtools-1.9/samtools"
bam_file21="/Users/chin-hua-huang/project_inputs/NA12878.chrom21.CGI.teramap.CEU.high_coverage.20120417.bam"
#bam_file22="/Users/chin-hua-huang/Desktop/Sequence_Bioinformatics/Assignments/Projekt"
sam_out="/Users/chin-hua-huang/project_outputs/output_chr21.sam"
ref21="/Users/chin-hua-huang/project_inputs/chrom.fasta"
ref21_fai="/Users/chin-hua-huang/project_inputs/chrom.fasta.fai"

#Variables for the simulated bam:

#-----------------------------------------------------
#Convert bam file to sorted bam file, then index it:(This process is mandatory(on my Mac-PC...) for the softwares except Varscan and samtools. )
output_sorted_bam="/Users/chin-hua-huang/project_outputs/chr21.sorted.bam"
#$samtools sort $bam_file21 -o $output_sorted_bam -m 10G

#View the sorted bam file:
#$samtools view -h /Users/chin-hua-huang/project_outputs/chr21.sorted.bam.tmp.0000.bam | head

#Index the sorted bam file:
#$samtools index $bam_file21 


#-----------------------------------------------------------------------------------
#Preprocessing to bcf file then to vcf file(samtools-bcftools pipeline)

output_bcf="/Users/chin-hua-huang/project_outputs/output_chr21_10x.bcf"
#$samtools mpileup -g -f $ref21 $bam_file21 > $output_bcf

#Code for bcf tools to generate vcf file.
#Variables:
bcf=/Users/chin-hua-huang/bcftools-1.10/bcftools
output_vcf="/Users/chin-hua-huang/project_outputs/output_chr21_10x.vcf"

#Commandline:
#$bcf call -c -v $output_bcf > $output_vcf 

#----------
#Code for VarScan
#Variables:
output_mpileup="/Users/chin-hua-huang/project_outputs/output_chr21_10x.mpileup"
output_vcf_varscan="/Users/chin-hua-huang/project_outputs/output_chr21_10x_varscan.vcf"
output_vcf_varscan_indel="/Users/chin-hua-huang/project_outputs/output_chr21_10x_varscan_indel.vcf"
varscan="/Users/chin-hua-huang/VarScan.v2.3.9.jar"

#Commandlines: 1. line is for the preparation of mpileup file for the 2. and 3. commandlines:
#$samtools mpileup -f $ref21 $bam_file21 > $output_mpileup
#java -jar $varscan mpileup2snp $output_mpileup --output-vcf 1 > $output_vcf_varscan
#java -jar $varscan mpileup2indel $output_mpileup --output-vcf 1 > $output_vcf_varscan_indel

#-------------------
#Add read groups to bam file for GATK:

#Variables
picard="/Users/chin-hua-huang/picard/build/libs/picard.jar"
bam_rg="/Users/chin-hua-huang/project_inputs/bam_rg.bam"
bam_tmp="/Users/chin-hua-huang/project_outputs/chr21.sorted.bam.tmp.0000.bam"

#Two alternative Commandlines
#java -jar $picard AddOrReplaceReadGroups CREATE_INDEX=true I=$bam_tmp O=$bam_rg RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
#java -jar $picard AddOrReplaceReadGroups I=$bam_tmp O=$bam_rg RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

#------------
#Code for GATK

#Variables:
output_gatk_vcf="/Users/chin-hua-huang/project_outputs/output21_gatk.vcf"
gatk="/Users/chin-hua-huang/gatk-4.1.4.1/gatk"

#Commandline:
#$gatk HaplotypeCaller -R $ref21 -I $bam_tmp -O $output_gatk_vcf

#------------
#Code for Platypus

#Variables
bam_file21_sorted="/Users/chin-hua-huang/project_outputs/chr21.sorted.bam.tmp.0000.bam"
output_platypus_vcf="/Users/chin-hua-huang/project_outputs/output_platypus.vcf"
platypus="/Users/chin-hua-huang/Platypus/bin/Platypus.py"

#Commandline:
#python $platypus callVariants --bamFiles=$bam_rg --refFile=$ref21 --output=$output_platypus_vcf

#----------------
#Code for FreeBayes

#Variables
output_freebayes_vcf="/Users/chin-hua-huang/project_outputs/output_freebayes.vcf"
freebayes="/Users/chin-hua-huang/freebayes/bin/freebayes"

#Commandline
#$freebayes -f $ref21 $bam_file21 > $output_freebayes_vcf



#-----------------
#Code for vcftools to compare the different vcf-files:
#variables
vcftools="/Users/chin-hua-huang/vcftools_0.1.13/bin/vcftools"
ground_truth="/Users/chin-hua-huang/project_outputs/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_all.vcf"
varscan_indel_vcf="/Users/chin-hua-huang/project_outputs/output_chr21_10x_varscan_indel.vcf"
varscan_vcf="/Users/chin-hua-huang/project_outputs/output_chr21_10x_varscan.vcf"
samtools_vcf="/Users/chin-hua-huang/project_outputs/output_chr21_10x.vcf"
compare_varscan_indel="/Users/chin-hua-huang/project_outputs/compare_varscan_indel"
out_vcf_indel="/Users/chin-hua-huang/project_outputs/picard_indel_vcf"

ground_truth_sim="/Users/chin-hua-huang/project_inputs/truevar.txt"
freebayes_sim="/Users/chin-hua-huang/project_inputs/freebayes.vcf"
gatk_vcf="/Users/chin-hua-huang/project_outputs/gatk_sim.vcf"


#Parallel commands:
#$vcftools --vcf $ground_truth --diff $varscan_indel_vcf --diff-site --out $compare_varscan_indel --chr chr21
#$vcftools --vcf $ground_truth --diff $varscan_vcf --diff-site --out $compare_varscan_indel --chr chr21
#$vcftools --vcf $ground_truth --diff $samtools_vcf --diff-site --out $compare_varscan_indel --chr chr21

#$vcftools --vcf $gatk_vcf --diff $freebayes_sim --diff-site --out $compare_varscan_indel --chr simref


#-------------------
sample="/Users/chin-hua-huang/project_inputs/example.sam"
sample_bam="/Users/chin-hua-huang/project_inputs/example.bam"
ref="/Users/chin-hua-huang/project_inputs/ref.fa"
vcf_varscan_sim="/Users/chin-hua-huang/project_outputs/varscan_sim.vcf"
vcf_varscan_sim_indel="/Users/chin-hua-huang/project_outputs/varscan_sim_indel.vcf"
output_dict="/Users/chin-hua-huang/project_outputs/output.dict"

#$samtools view -S -b $sample > $sample_bam

#$samtools mpileup -f $ref $sample_bam > $output_mpileup
#java -jar $varscan mpileup2snp $output_mpileup --output-vcf 1 > $vcf_varscan_sim
#java -jar $varscan mpileup2indel $output_mpileup --output-vcf 1 > $vcf_varscan_sim_indel

#$vcftools --vcf $vcf_varscan_sim --diff $freebayes_sim --diff-site --out $compare_varscan_indel --chr simref
#$gatk HaplotypeCaller -R $ref21 -I $bam_tmp -O $output_gatk_vcf


#java -jar $picard CreateSequenceDictionary REFERENCE=$ref OUTPUT=$output_dict

#--------------

#Preprocessing to bcf file then to vcf file(samtools-bcftools pipeline)

sample_bcf="/Users/chin-hua-huang/project_outputs/sample.bcf"
$samtools mpileup -g -f $ref $sample_bam > $sample_bcf

#Code for bcf tools to generate vcf file.
#Variables:
bcf=/Users/chin-hua-huang/bcftools-1.10/bcftools
output_vcf="/Users/chin-hua-huang/Desktop/sample.vcf"

#Commandline:
$bcf call -c -v $sample_bcf > $output_vcf

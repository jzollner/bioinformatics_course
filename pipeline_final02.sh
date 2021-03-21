#### Pipeline for NGS ###
##Advanced bioinformatics course
##Create the suitable environment with the appropriate folders
##See read.me with the description of folders and files
##Prealignment QC
##1.0 FASTQC on raw data
##please ensure that the files are unzipped using zcat

fastqc -t 4 ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/WES01_chr22m_R1.fastq.gz ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/WES01_chr22m_R2.fastq.gz
mkdir ~/ngs_course/dnaseq_pipeline/results/fastqc_untrimmed_reads
mv ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/*fastqc* ~/ngs_course/dnaseq_pipeline/results/fastqc_untrimmed_reads
#2.0 Trimmomatics
trimmomatic PE \
-threads 4 \
-phred33 \
/home/ubuntu/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/WES01_chr22m_R1.fastq.gz \
/home/ubuntu/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/WES01_chr22m_R2.fastq.gz \
-baseout /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/WES01_chr22m_trimmed_R \
ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 \
TRAILING:25 MINLEN:50

##3.0 FASTQC on trimmed data
fastqc -t 4 /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/WES01_chr22m_trimmed_R_1P \
/home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/WES01_chr22m_trimmed_R_2P

mkdir /home/ubuntu/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads

mv /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/*fastqc* /home/ubuntu/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads/

cd ~/ngs_course/dnaseq_pipeline/data/reference


#reference genome indexing ===================
#==============================================
#This i did in a different step - to allow this 
#indexed files are within /reference folder
#the usual fa, bwt, sa, etc
#===========================
#5.0 running BWA mem
# using the RG information

mkdir ~/ngs_course/dnaseq_pipeline/data/aligned_data

bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.WES01\tSM:WES01\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50  ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/WES01_chr22m_trimmed_R_1P ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/WES01_chr22m_trimmed_R_2P > ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m.sam

cd ~/ngs_course/dnaseq_pipeline/data/aligned_data

##using samtools
#6.0

samtools view -h -b WES01_chr22m.sam > WES01_chr22m.bam
samtools sort WES01_chr22m.bam > WES01_chr22m_sorted.bam
samtools index WES01_chr22m_sorted.bam

##post alignment filtering and QC
#7.0

picard MarkDuplicates -I WES01_chr22m_sorted.bam -O WES01_chr22m_sorted_marked.bam -M marked_dup_metrics.txt
samtools index WES01_chr22m_sorted_marked.bam

#Filter BAM based on mapping quality and bitwise flags using samtools
#8.0

samtools view -F 1796  -q 20 -o WES01_chr22m_sorted_filtered.bam WES01_chr22m_sorted_marked.bam
samtools index WES01_chr22m_sorted_filtered.bam

#variant calling with freebayes
#9.0

zcat ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa
samtools faidx ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa
freebayes --bam ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa --vcf ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf
bgzip ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf
tabix -p vcf ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf.gz

##10 Filtering VCF
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf.gz > ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered.vcf

##use bedtools now
bedtools intersect -header -wa -a ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered.vcf -b ~/ngs_course/dnaseq_pipeline/data/chr22.genes.hg19.bed  > ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf

bgzip ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf

tabix -p vcf ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf.gz

##11.0 ANNOVAR

cd ~/

cd annovar

##VCF to ANNOVAR input
./convert2annovar.pl -format vcf4 /home/ubuntu/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf.gz > /home/ubuntu/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.avinput

#Make a table format
#final step
#12.0

./table_annovar.pl /home/ubuntu/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.avinput humandb/ -buildver hg19 -out /home/ubuntu/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22 -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout


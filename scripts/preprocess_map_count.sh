#!/bin/bash -login
#PBS -l nodes=1:ppn=10,walltime=02:00:00:00,mem=20gb
#PBS -j oe
set -e
set -u
set -o pipefail

# working directory:
cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/RNA_DE

# Look at your data in FastQC and generate a report
module load FastQC/0.11.3

find ./data/${Var}/ -name '*_L00*_R*.fastq.gz' | \
xargs -i sh -c "fastqc {} -o ./analysis/fastqc/"

# Trim the data with trimmomatic
module load Trimmomatic/0.33

find ./data/${Var}/ -name '*_L00*_R*.fastq.gz' | \
xargs -i basename {} | \
sed 's/_R[1-2].fastq.gz//' | \
sort | uniq | \
xargs -i sh -c "java -jar -Xmx20g $TRIM/trimmomatic PE -threads 8 ./data/${Var}/{}_R1.fastq.gz ./data/${Var}/{}_R2.fastq.gz ./data/${Var}/{}_R1_paired.fastq.gz ./data/${Var}/{}_R1_unpaired.fastq.gz ./data/${Var}/{}_R2_paired.fastq.gz ./data/${Var}/{}_R2_unpaired.fastq.gz ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE.fa:2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:25"

# This is equivalent to:

#java -jar $TRIM/trimmomatic \
#PE \
#./data/File1_L001_R1.fastq.gz \
#./data/File1_L001_R2.fastq.gz \
#./data/File1_L001_R1_paired.fastq.gz \
#./data/File1_L001_R1_unpaired.fastq.gz \
#./data/File1_L001_R2_paired.fastq.gz \
#./data/File1_L001_R2_unpaired.fastq.gz \
#ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE.fa:2:30:10 \
#LEADING:2 \
#TRAILING:2 \
#SLIDINGWINDOW:4:2 \
#MINLEN:25

# Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
# ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
# Remove Illumina adapters provided in the TruSeq3-PE.fa file (provided). Initially Trimmomatic will look for seed matches (16 bases) allowing maximally 2 mismatches. 
# These seeds will be extended and clipped if in the case of paired end reads a score of 30 is reached (about 50 bases), 
# or in the case of single ended reads a score of 10, (about 17 bases).
# Remove leading low quality or N bases (below quality 2) (LEADING:2)
# Remove trailing low quality or N bases (below quality 2) (TRAILING:2)
# Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 2 (SLIDINGWINDOW:4:2)
# Choose a phed cutoff of two based on MacManes, 2014.
# Drop reads below the 25 bases long (MINLEN:25)

# Remove the trimmed unpaired reads
find ./data/${Var}/ -name '*_L00[1-8]_R[1-2]_unpaired.fastq.gz' | \
xargs -i sh -c "rm {}"

# Remove the original raw data files that were previously copied
find ./data/${Var}/ -name '*_L00*_R[1-2].fastq.gz' | \
xargs -i sh -c "rm {}"

# Look at paired trimmed reads in FastQC
find ./data/${Var}/ -name '*_L00[1-8]_R[1-2]_paired.fastq.gz' | \
xargs -i sh -c "fastqc {} -o ./analysis/fastqc/"

# Obtained chicken specific rRNA from GenBank accession number KT445934
# file: /mnt/scratch/steepale/birdman/MDV_project/RNA_dif_express/data/ref/KT445934.2/KT445934.2_Gallus_rRNA_gg4.fasta

# Examine the rRNA data to make sure it is in proper format for bbmap
/mnt/home/steepale/Apps/bbmap/testformat.sh \
./data/ref/KT445934.2/KT445934.2_Gallus_rRNA_gg4.fasta
#sanger	fasta	raw	single-ended

# Remove the rRNA sequences specific to chicken
find ./data/${Var}/ -name '*_L00[1-8]_R[1-2]_paired.fastq.gz' | \
xargs -i basename {} | \
sed 's/_R[1-2]_paired.fastq.gz//' | \
sort | uniq | 
xargs -i sh -c "/mnt/home/steepale/Apps/bbmap/bbduk.sh -Xmx18g in1=./data/${Var}/{}_R1_paired.fastq.gz in2=./data/${Var}/{}_R2_paired.fastq.gz out1=./data/${Var}/{}_R1_paired_norRNA.fastq.gz out2=./data/${Var}/{}_R2_paired_norRNA.fastq.gz ref=./data/ref/KT445934.2/KT445934.2_Gallus_rRNA_gg4.fasta"

# This is equivalent to:
#/mnt/home/steepale/Apps/bbmap/bbduk.sh \
#in1=./data/Sample_lane_R1_paired.fastq \
#in2=./data/Sample_lane_R2_paired.fastq \
#out1=./data/Sample_lane_R1_paired_norRNA.fastq \
#out2=./data/Sample_lane_R2_paired_norRNA.fastq \
#ref=./data/ref/KT445934.2/KT445934.2_Gallus_rRNA_gg4.fasta

# Remove the redundent paired trimmed reads that were output from Trimmomatic
find ./data/${Var}/ -name '*_L00[1-8]_R[1-2]_paired.fastq.gz' | \
xargs -i sh -c "rm {}"

# Look at paired trimmed reads minus rRNA in FastQC
find ./data/${Var}/ -name '*_L00[1-8]_R[1-2]_paired_norRNA.fastq.gz' | \
xargs -i sh -c "fastqc {} -o ./analysis/fastqc/"

# Map the reads to a reference genome
module load TopHat/2.1.0
module load bowtie2/2.2.6

# Need to map by certain paired lanes (1&2, 3&4, 5&6, 7) based on how samples were organized during sequencing run
#1&2
find ./data/${Var}/ -name '*_L00[1-2]_R[1-2]_paired_norRNA.fastq.gz' | \
xargs -i basename {} | \
sed 's/_L00[1-8]_R[1-2]_paired_norRNA.fastq.gz//' | \
sort | uniq | \
xargs -i sh -c "tophat -p 10 --transcriptome-index=./data/ref/GCF_000002315.4_Gallus_gallus-5.0_genomic -o ./data/${Var} ./data/ref/galgal5 ./data/${Var}/{}_L001_R1_paired_norRNA.fastq.gz,./data/${Var}/{}_L002_R1_paired_norRNA.fastq.gz ./data/${Var}/{}_L001_R2_paired_norRNA.fastq.gz,./data/${Var}/{}_L002_R2_paired_norRNA.fastq.gz"
#3&4
find ./data/${Var}/ -name '*_L00[3-4]_R[1-2]_paired_norRNA.fastq.gz' | \
xargs -i basename {} | \
sed 's/_L00[1-8]_R[1-2]_paired_norRNA.fastq.gz//' | \
sort | uniq | \
xargs -i sh -c "tophat -p 10 --transcriptome-index=./data/ref/GCF_000002315.4_Gallus_gallus-5.0_genomic -o ./data/${Var} ./data/ref/galgal5 ./data/${Var}/{}_L003_R1_paired_norRNA.fastq.gz,./data/${Var}/{}_L004_R1_paired_norRNA.fastq.gz ./data/${Var}/{}_L003_R2_paired_norRNA.fastq.gz,./data/${Var}/{}_L004_R2_paired_norRNA.fastq.gz"
#5&6
find ./data/${Var}/ -name '*_L00[5-6]_R[1-2]_paired_norRNA.fastq.gz' | \
xargs -i basename {} | \
sed 's/_L00[1-8]_R[1-2]_paired_norRNA.fastq.gz//' | \
sort | uniq | \
xargs -i sh -c "tophat -p 10 --transcriptome-index=./data/ref/GCF_000002315.4_Gallus_gallus-5.0_genomic -o ./data/${Var} ./data/ref/galgal5 ./data/${Var}/{}_L005_R1_paired_norRNA.fastq.gz,./data/${Var}/{}_L006_R1_paired_norRNA.fastq.gz ./data/${Var}/{}_L005_R2_paired_norRNA.fastq.gz,./data/${Var}/{}_L006_R2_paired_norRNA.fastq.gz"
#7
find ./data/${Var}/ -name '*_L007_R[1-2]_paired_norRNA.fastq.gz' | \
xargs -i basename {} | \
sed 's/_L00[1-8]_R[1-2]_paired_norRNA.fastq.gz//' | \
sort | uniq | \
xargs -i sh -c "tophat -p 10 --transcriptome-index=./data/ref/GCF_000002315.4_Gallus_gallus-5.0_genomic -o ./data/${Var} ./data/ref/galgal5 ./data/${Var}/{}_L007_R1_paired_norRNA.fastq.gz ./data/${Var}/{}_L007_R2_paired_norRNA.fastq.gz"

# This is equivalent to:
#tophat \
#-p 4 \
#--transcriptome-index=./data/ref/GCF_000002315.4_Gallus_gallus-5.0_genomic \
#-o ./data/sample_barcode \
#./data/ref/galgal5 \
#./data/sample_barcode_firstlane_R1_paired_norRNA.fastq.gz,./data/sample_barcode_secondlane_R1_paired_norRNA.fastq.gz \
#./data/sample_barcode_firstlane_R2_paired_norRNA.fastq.gz,./data/sample_barcode_secondlane_R2_paired_norRNA.fastq.gz

# Count the number of mapped reads with samtools
module load SAMTools/1.2

find ./data/${Var}/ -name '*_L00[1-8]_R[1-2]_paired_norRNA.fastq.gz' | \
xargs -i basename {} | \
sed 's/_L00[1-8]_R[1-2]_paired_norRNA.fastq.gz//' | \
sort | uniq | \
xargs -i sh -c "echo {}; samtools view -c -F 4 ./data/${Var}/accepted_hits.bam"

# Look at the mapping summary files
find ./data/${Var}/ -name '*_L00[1-8]_R[1-2]_paired_norRNA.fastq.gz' | \
xargs -i basename {} | \
sed 's/_L00[1-8]_R[1-2]_paired_norRNA.fastq.gz//' | \
sort | uniq | \
xargs -i sh -c "echo {}; cat ./data/${Var}/align_summary.txt"

# 60% of reads mapped when all said and done... 60% of 50,000,000 is 30,000,000... literally just enough

# Sort the BAMs with Samtools
find ./data/${Var}/ -name '*_L00[1-8]_R[1-2]_paired_norRNA.fastq.gz' | \
xargs -i basename {} | \
sed 's/_L00[1-8]_R[1-2]_paired_norRNA.fastq.gz//' | \
sort | uniq | \
xargs -i sh -c "samtools sort -n ./data/${Var}/accepted_hits.bam ./data/${Var}/{}_paired_norRNA_sorted"

# Remove the tophat output bam file
rm ./data/${Var}/accepted_hits.bam
rm ./data/${Var}/unmapped.bam

# Count the mapped RNA read coverage on respective genes
module load PySAM/0.6
module load HTSeq/0.6.1

find ./data/${Var}/ -name '*_paired_norRNA_sorted.bam' | \
xargs -i basename {} | \
sed 's/_paired_norRNA_sorted.bam//' | \
sort | uniq | \
xargs -i sh -c "htseq-count --format=bam --stranded=yes --order=name ./data/${Var}/{}_paired_norRNA_sorted.bam ./data/ref/ensemble/Gallus_gallus.Gallus_gallus-5.0.86.gtf > ./data/gene_counts/{}_ensembl_gene_counts.txt"

# Remove the paired trimmed reads minus rRNA
find ./data/${Var}/ -name '*_L00[1-8]_R[1-2]_paired_norRNA.fastq.gz' | \
xargs -i sh -c "rm {}"

qstat -f ${PBS_JOBID}

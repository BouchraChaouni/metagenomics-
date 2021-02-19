# metagenomics-
coding
# Illumina Miseq paired-end reads with an insert size of 300bp and 2.000.000 reads per sample
# Kraken: Taxonomic classifier 
# Build Kraken Database: 
module load kraken/0.10.5-beta
kraken-build --standard --db Kraken
# Run Kraken:
kraken --db ~/Kraken --threads 1 --fastq-input osd24_2014_SM_1.fastq osd24_2014_SM_2.fastq --output osd24_kraken.out.txt
kraken --db ~/Kraken --threads 1 --fastq-input osd91_2014_SM_1.fastq osd91_2014_SM_2.fastq --output osd91_kraken.out.txt
# Getting taxonomic name associated with input sequences: 
kraken-translate --db ~/Kraken osd24_kraken.out.txt > osd24_kraken_labels.txt
kraken-translate --db ~/Kraken osd91_kraken.out.txt > osd91_kraken_labels.txt

# Run FastQC: quality control tool for high throughput sequence data
module load fastqc/0.11.4
fastqc -o fastqc *fastq
# Since I only got a small fraction of reads classified, and most of them assigned to phages. We doubted it might be some sort of contamination, but the fastqc showed a clean output. We will run Assembly next.  

# Run Assembly: 
module load spades/3.11.0
/bioware/spades-3.11.0/bin/metaspades.py -1 /users/bchaouni/CGSB_BC/sm_raw_data/osd24_2014_SM_1.fastq -2 /users/bchaouni/CGSB_BC/sm_raw_data/osd24_2014_SM_2.fastq -t 24 -m 96 -o /users/bchaouni/CGSB_BC/sm_raw_data/osd24_metaspades
/bioware/spades-3.11.0/bin/metaspades.py -1 /users/bchaouni/CGSB_BC/sm_raw_data/osd91_2014_SM_1.fastq -2 /users/bchaouni/CGSB_BC/sm_raw_data/osd91_2014_SM_2.fastq -t 24 -m 96 -o /users/bchaouni/CGSB_BC/sm_raw_data/osd91_metaspades
# Blasting the first 10 contigs with an average coverage of 10 didn't give a 100% match, and they were assigned to synechococcus phages. Thus, we will use bwa for assembly validation and Trimmomatic then run metaspades one more time 

# Run bwa: 
module load bwa/0.7.15
# build aligner specific reference index
bwa index contigs_10.fasta
# align reads to reference 
bwa mem contigs_10.fasta /users/bchaouni/CGSB_BC/sm_raw_data/osd24_2014_SM_1.fastq /users/bchaouni/CGSB_BC/sm_raw_data/osd24_2014_SM_2.fastq > osd24_contigs_10.sam 
# convert SAM to BAM
module load samtools/1.3.1 
# Verify if the header starts with @
samtools view -bS osd24_contigs_10.sam > osd24_contigs_10.bam
# sort BAM by position 
samtools sort osd24_contigs_10.bam -o osd24_contigs_10.sorted.bam
# index BAM
samtools index osd24_contigs_10.sorted.bam
# display BAM content 
samtools view osd24_contigs_10.sorted.bam | head
# alignment metrics and QC
samtools flagstat osd24_contigs_10.sorted.bam
samtools idxstats osd24_contigs_10.sorted.bam

# Run Trimmomatic:
java -jar /bioware/trimmomatic-0.36/trimmomatic-0.36.jar PE osd24_2014_SM_1.fastq osd24_2014_SM_2.fastq osd24_trimmomatic_1P.fq osd24_trimmomatic_1U.fq osd24_trimmomatic_2P.fq osd24_trimmomatic_2U.fq ILLUMINACLIP:/users/bchaouni/CGSB_BC/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar /bioware/trimmomatic-0.36/trimmomatic-0.36.jar PE osd91_2014_SM_1.fastq osd91_2014_SM_2.fastq osd91_trimmomatic_1P.fq osd91_trimmomatic_1U.fq osd91_trimmomatic_2P.fq osd91_trimmomatic_2U.fq ILLUMINACLIP:/users/bchaouni/CGSB_BC/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
wc osd*q -l
# Run metaspades with only Trimmomatic paired files:
/bioware/spades-3.11.0/bin/metaspades.py  --pe1-1 /users/bchaouni/CGSB_BC/sm_raw_data/trimmomatic_rerun/osd24_trimmomatic_1P.fq   --pe1-2 /users/bchaouni/CGSB_BC/sm_raw_data/trimmomatic_rerun/osd24_trimmomatic_2P.fq   -t 24 -m 96 -o      /users/bchaouni/CGSB_BC/sm_raw_data/trimmomatic_rerun/osd24_trimmomatic_metaspades
# Didn't work with unpaired !!!

# Run QUAST 
module load quast/intel/4.4
quast.py /scratch/bc2609/CGSB_BC/sm_raw_data/osd24_metaspades/
contigs_10.fasta -1 /scratch/bc2609/CGSB_BC/sm_raw_data/osd24_2014_SM_1.fastq -2 /scratch/bc2609/CGSB_BC/sm_raw_data/osd24_2014_SM_2.fastq —bam /scratch/bc2609/CGSB_BC/sm_raw_data/osd24_metaspades/
osd24_contigs_10.sorted.bam --bam  /scratch/bc2609/CGSB_BC/sm_raw_data/osd24_metaspades/
osd24_contigs_10.sorted.bam.bai -o quast_output_osd24
# View the html report
module load firefox/55.0
firefox report.html

# Run IGV
srun --x11 --mem=62GB --time=03:00:00 --cpus-per-task=20 --pty /bin/bash
module load igv/2.3.90
igv

# Run GOTTCHA
module load gottcha/intel/20171102
gottcha.pl --outdir gottcha_osd24 --input osd24_2014_SM_1.fastq,osd24_2014_SM_2.fastq --mode all --database /scratch/at120/bouchra-metagenomics/gottcha-db/database/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species

# Run Krona
module load kronatools/2.7
share/apps/kronatools/2.7/scripts/ImportTaxonomy.pl /scratch/at120/bouchra-metagenomics/gottcha_osd24/osd24_2014_SM_merged_temp/osd24_2014_SM_merged.family.tsv -o osd24.krona.html
# Open html report with firefox

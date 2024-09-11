#### Run stats for each HQ MAG ####

``for i in <HQ_MAGs_directory>/*.fasta; do stats.sh in=${i} format=5 > ${i}_stats.txt; done &`` 

#### Coverage of HQ Nelore MAGs (required by NCBI) ####

``cat <HQ_MAGs_directory>/*.fasta > HQ_MAGs.fasta``

``bowtie2-build HQ_MAGs.fasta HQ_MAGs --threads 24``

``bowtie2 --threads 24 -x HQ_MAGs.fasta -1 <reads_R1>.fastq -2 <reads_R2>.fastq --no-unal | samtools view -@ 24 -F 4 -bS > <sample>_HQ_mags.bam``

``samtools sort -@ 24 <sample>_HQ_mags.bam -o <sample>_HQ_mags_sorted.bam``

``coverm genome --bam-files <sample>_HQ_mags_sorted.bam -m mean --min-read-aligned-percent 0.75 --min-read-percent-identity 0.95 --min-covered-fraction 0 --genome-fasta-directory <HQ_MAGs_directory> --genome-fasta-extension fasta -t 24 --sharded > <sample>_HQ_mags_coverm_mean.txt``

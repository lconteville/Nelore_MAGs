#### Metagenomic Assembly ####

``megahit/bin/megahit -1 <reads_R1>.fastq -2 <reads_R2>.fastq --kmin-1pass -m 60e+10 --k-list 27,37,47,57,67,77,87 --min-contig-len 1000 -t 16 -o <megahit_output>``

#### Create index for contigs fasta ####

``bwa index <megahit_output>/<contigs>.fa``

#### Mapping Reads against Contigs ####

``bwa mem -t 24 <megahit_output>/<contigs>.fa <reads_R1>.fastq <reads_R2>.fastq | samtools view -@ 4 -bS - > <mapping>.bam``

#### Sort Bam Files ####

``samtools sort -@ 24 -m 2G -o <sorted>.bam <mapping>.bam``

#### Index Sorted Bam Files ####

``samtools index <sorted>.bam``

#### Generate Coverage Files for Metabat and MaxBin2####

``jgi_summarize_bam_contig_depths --outputDepth <depth_file>.txt <sorted>.bam``

``cut -f1,4 <depth_file>.txt > <depth_file_maxbin>.txt``

#### Generate Coverage Files for Concoct ####

``cut_up_fasta.py <contigs>.fa -c 10000 -o 0 --merge_last -b <contigs_10K>.bed > <contigs_10K>.fa``

``concoct_coverage_table.py <contigs_10K>.bed <sorted>.bam > <coverage_table>.tsv``

#### Run Metabat2 ####

``metabat2 -i <contigs>.fa -a <depth_file>.txt -o <metabat_output/bins> -m 2000 --minContigDepth 2``

#### Run Maxbin2 ####

``run_MaxBin.pl -thread 24 -min_contig_length 1500 -contig <contigs>.fa -abund <depth_file_maxbin>.txt -out <maxbin_output>``

#### Run Concoct ####

``concoct --composition_file <contigs_10K>.fa --coverage_file <coverage_table>.tsv --threads 24 -b <concoct_output>``


#### Generate inputs for DAS Tool ####

``Fasta_to_Contig2Bin.sh -i <metabat_output> -e fa > <metabat.contigs2bin>.tsv``

``Fasta_to_Contig2Bin.sh -i <maxbin_output> -e fasta > <maxbin.contigs2bin>.tsv``

``cp clustering_gt1000.csv > <concoct.contigs2bin>.tsv``

#### Run DAS Tool ####

``DAS_Tool -i <metabat.contigs2bin>.tsv,<maxbin.contigs2bin>.tsv,<concoct.contigs2bin>.tsv -l metabat,maxbin,concoct -c <contigs>.fa -t 24 --search_engine diamond --write_bin_evals --write_bins -o <DASTool_output>``

#### Run dRep ####

``dRep dereplicate -p 32 -comp 50 -con 10 -pa 0.95 -sa 0.99 <dRep_output> -g <DASTool_output>/*fa``

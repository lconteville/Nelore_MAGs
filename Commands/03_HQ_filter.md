#### Generate table with MAGsID, Completeness >= 80, Contamination <= 10 and Quality Score >=50 ####

``awk -F',' '{print $1"\t"$3"\t"$4"\t"$3-5*$4}' <dRep_output>/data_tables/Widb.csv | awk -F$'\t' '$2 >= 80 && $3 <= 10 && $4>=50' > <dRep_output>/data_tables/Widb_80_10_50.tsv``

#### Create a directory exclusively with the high-quality MAGs ####

``mkdir <HQ_MAGs_directory>``

``for k in `awk '{print $1}' <dRep_output>/data_tables/Widb_80_10_50.tsv`; do
	cp <dRep_output>/dereplicated_genomes/${k} <HQ_MAGs_directory>;
done
``

#### Run GTDB tk ####

``export GTDBTK_DATA_PATH=~/release207_v2``

``gtdbtk classify_wf --genome_dir <HQ_MAGs_directory> -x fa --out_dir <GTDBtk_output> --skip_ani_screen --full_tree --cpus 32 --pplacer_cpus 12 ``

#### Generate a smaller bacterial tree for figure ####

``gtdbtk infer --msa_file <GTDBtk_output>/align/gtdbtk.bac120.user_msa.fasta --out_dir <GTDBtk_output> --prefix gtdbtk.smaller_tree --cpus 32``

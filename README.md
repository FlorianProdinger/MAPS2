Welcome to MAPS2!

The main script "MAPS2.tcsh" is the MAPS2 pipeline.
It calls in succession these steps:


 1) The "R_dada_script" this script does most of the work
 2) blastx on the ASVs to check if they are viral sequences AND blast returns the amino acid sequence
 3) a bash command "sed 's/^/>/' {$MAIN_DIR}$BLASTX_OUT..." to make a fasta file from the blast output
 4) mafft to add viral sequences to a reference file
 5) pplacer to place the amplpicon files in endo 2020's reference tree
 6) The "R_tree_script" this script returns a filtered ASV table and a statistics table to show were reads were lost
    (no singeltons, no ASVs placed outside assigned outside the Mimi branch of the tree) 
 7) The "trim_MIMI_ASV_to_commom_region_v4.py" to trimm the ASV in a common region (common region hardcoded in this pipeline) 
 8) cd-hit to cluster at 100% (and other percentages)
 9) The "R_script_OTU_tables" script makes ASV / OTU tables from the cd-hit output.

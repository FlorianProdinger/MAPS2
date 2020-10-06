#!/bin/tcsh

#This script is the pipeline it calls in succession the
# 1) The "R_dada_script" this script does most of the work
# 2) blastx on the ASVs to check if they are viral sequences AND blast returns the amino acid sequence
# 3) a bash command "sed 's/^/>/' {$MAIN_DIR}$BLASTX_OUT..." to make a fasta file from the blast output
# 4) mafft to add viral sequences to a reference file
# 5) pplacer to place the amplpicon files in endo 2020's reference tree
# 6) The "R_tree_script" this script returns a filtered ASV table and a statistics table to show were reads were lost
#    (no singeltons, no ASVs placed outside assigned outside the Mimi branch of the tree) 
# 7) The "trim_MIMI_ASV_to_commom_region_v4.py" to trimm the ASV in a common region (common region hardcoded in this pipeline) 
# 8) cd-hit to cluster at 100% (and other percentages)
# 9) The "R_script_OTU_tables" script makes ASV / OTU tables from the cd-hit output.

#cdb or APC
#PBS -q cdb 
#PBS -o qsub_MAPS2_OUT.txt
#PBS -e qsub_MAPS2_ERROR.txt
#PBS -N MP_dada2_pipe
#PBS -l nice=0
#PBS -l mem=150gb
#PBS -l ncpus=16
#PBS -m ea
#PBS -M florian.prodinger@gmx.net

#eval `/usr/bin/modulecmd tcsh load blast+/2.9.0`

#creating environment for modules to work
#check which command works
#eval `/usr/bin/modulecmd tcsh load blast+/2.9.0 pplacer/1.1.alpha19 mafft/7.453`
source /etc/profile.d/modules.sh
module purge
module load blast+/2.9.0
module load pplacer/1.1.alpha19
module load mafft/7.453
module load R/3.6.1
module load Python/3.7.5
module load cd-hit/4.6.1


####################
# define variables # 
####################

######################
# User input section #
######################

set THREADS_LIMIT="16"  #match with the PBS command!

#these are the input file and directory for this pipeline:
set MAIN_DIR="/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/20201005_test_xia_data/"
#"/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/20200611_pipeline_test/"
set R_path_to_raw="/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/20201005_test_xia_data/0_RAW_TEST_3_FILES/"
#"/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/test_qsub_20200527/0_RAW/"
#"/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/20200609_pipeline/"

######################
# User input section #
#        END         #
######################


#directory of the pipeline
set MAPS2_DIR="/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/20201005_MAPS2_pipeline/MAPS2/"
set R_dada_script="${MAPS2_DIR}MAPS2_1_dada2.R" #MP_1_dada2_ASvs_for_qsub.R"
set R_tree_script="${MAPS2_DIR}MAPS2_2_filter_ASV_table.R"  #MP_4_dada2_filter_jplace_for_qsub_v2.R"
set Python_script_name="${MAPS2_DIR}MAPS2_3_trim_MIMI_ASV_to_commom_region.py" #trim_MIMI_ASV_to_commom_region_v4.py"
set python_add_OTUs_to_table_script="${MAPS2_DIR}MAPS2_4_add_OTUs_to_final_table.py" #add_OTUs_to_final_table_v1.py"
set R_script_OTU_tables="${MAPS2_DIR}MAPS2_5_generate_final_ASV_table.R" #generate_ASV_and_OTU_tables_v1.R"

#names of the generated files
set R_dada_script_out_name="20201005_MAPS2_dada2_out"
set QUERY_IN="${R_dada_script_out_name}_DNA_sequences.fasta"
set AA_FASTA_FILE={$QUERY_IN}_AA.fasta
set BLASTX_OUT={$QUERY_IN}_BLASTX_out.txt  


#directory of the references 
set REFERENCES_DIR="/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/20201005_MAPS2_pipeline/MAPS2/reference_MAPS2/"
#"/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/dada2_references/"
set REF_PPLACER_ALN=${REFERENCES_DIR}"20200601_MAPS2_MIMI_TARA_EUK_BAC_ref_seqs_alingment.fasta"
set REF_PPLACER_TREE=${REFERENCES_DIR}RAxML_bestTree.20200601_MAPS2_MIMI_TARA_EUK_BAC_ref_seqs_pplacer
set REF_PPLACER_INFO=${REFERENCES_DIR}RAxML_info.20200601_MAPS2_MIMI_TARA_EUK_BAC_ref_seqs_pplacer
set REF_POLB_HOMOLOGY_SEARCH=${REFERENCES_DIR}PolB_homology_search.faa

#set REF_POLB_HOMOLOGY_SEARCH=${REFERENCES_DIR}MAPS_references/PolB_homology_search.faa
#
#set REF_PPLACER_TREE=${REFERENCES_DIR}20200601_MAPS2_MIMI_TARA_EUK_BAC_ref_seqs_nwk.tre
#set REF_PPLACER_ALN=${REFERENCES_DIR}20200601_MAPS2_MIMI_TARA_EUK_BAC_ref_seqs_alingment.phylip




#####################
# starting pipeline #
#####################

#calling dada
#script takes
#   output_directory  
#   output name of files > these files will be used by all other scripts!!
#   path to the raw files
#   how many threads can be used 

#echo "Rscript - dada2 pipeline: import, filterAndTrim, learnErrors, ASV formation (dada2), chimera filter, output: table and fasta file"
#Rscript \
#   $R_dada_script \
#   $MAIN_DIR \
#   $R_dada_script_out_name \
#   $R_path_to_raw \
#   $THREADS_LIMIT 
#if ($? != 0) then
# echo " command exit status after previous command is $?"
# exit
#endif
#
###calls blastx on the ASV sequences that were generated by dada2 in R
###needs more than 50gb in memory
###$QUERY_IN is generated by dada2 R script
##making database:
##makeblastdb -dbtype prot \
## -in $REF_POLB_HOMOLOGY_SEARCH -out $REF_POLB_HOMOLOGY_SEARCH \
## -logfile ${REF_POLB_HOMOLOGY_SEARCH}_logfile -input_type fasta -parse_seqids
#echo "blastx - search the ASV against a custom PolB aminoacid sequence database (Li&Hingamp et al 2018), saving aminoacid sequence of ASV"
#blastx \
# -query {$MAIN_DIR}$QUERY_IN \
# -db $REF_POLB_HOMOLOGY_SEARCH \
# -out {$MAIN_DIR}$BLASTX_OUT \
# -outfmt '10 qseqid qseq evalue pident qframe sseqid sallseqid' \
# -evalue 1e-5 \
# -max_target_seqs 1 \
# -num_threads $THREADS_LIMIT
#if ($? != 0) then
# echo " command exit status after previous command is $?"
# exit
#endif
#
###converting the blast output to a amino acid fasta file
###cat test_20200526_5_samples_250ASVs_BLASTX_out.txt | cut -d, -f1,2
#echo "extracting translated Mimiviridae sequences"
#sed 's/^/>/' {$MAIN_DIR}$BLASTX_OUT | cut -d, -f1,2 | tr , '\n' > {$MAIN_DIR}$AA_FASTA_FILE 
#if ($? != 0) then
# echo " command exit status after previous command is $?"
# exit
#endif
#
#
#echo "mafft - adding the amino acid sequence ASV to a reference alignment (Li&Hingamp et al 2018)"
##adding the translated Megaprimer reads to the reference
#mafft \
#  --thread $THREADS_LIMIT \
#  --quiet \
#  --6merpair\
#  --addfragments \
#  ${MAIN_DIR}$AA_FASTA_FILE $REF_PPLACER_ALN > ${MAIN_DIR}${AA_FASTA_FILE}_alingment.fasta
#if ($? != 0) then
# echo " command exit status after previous command is $?"
# exit
#endif
#
#echo "pplacer - placing the translated Megaprimer ASV in the reference tree (Endo et al 2020)"
##placing the translated Megaprimer reads in the reference tree
##using more than 1 core doesn't work with qsub. Seems like a bug
#pplacer -j 1 --timing\
# --verbosity 0 \
# ${MAIN_DIR}${AA_FASTA_FILE}_alingment.fasta \
# -t ${REF_PPLACER_TREE} \
# -s ${REF_PPLACER_INFO} \
# -o {$MAIN_DIR}${AA_FASTA_FILE}_pplacer.jplace 
#if ($? != 0) then
# echo " command exit status after previous command is $?"
# exit
#endif


echo "Rscript - removes sequences filtered by blast & pplacer"
##input: main directory, jplace file, dada2 ASV table from first R script 
Rscript \
 $R_tree_script \
  $MAIN_DIR \
  {$MAIN_DIR}{$AA_FASTA_FILE}_pplacer.jplace  \
  ${R_dada_script_out_name}_ASV_table.tsv \
  {$MAIN_DIR}$BLASTX_OUT 
if ($? != 0) then
 echo " command exit status after previous command is $?"
 exit
endif


echo "python script - trimming the ASVs to a common region"
##python ${Python_script_name} 20200528_dada2_out_DNA_sequences.fasta  20200528_dada2_out_DNA_sequences.fasta_AA.fasta_alingment.fasta  20200528_dada2_out_DNA_sequences.fasta_BLASTX_out.txt  23350  23970
python ${Python_script_name} \
 ${MAIN_DIR}{$QUERY_IN} \
 ${MAIN_DIR}{$AA_FASTA_FILE}_alingment.fasta \
 ${MAIN_DIR}{$BLASTX_OUT} \
 23350 \
 23970
if ($? != 0) then
 echo " command exit status after previous command is $?"
 exit
endif



echo "cd-hit - generate OTUs at several percentages using cd-hit"
/bin/cp ${MAIN_DIR}${R_dada_script_out_name}_final_ASV_table_noS.tsv ${MAIN_DIR}${R_dada_script_out_name}_final_ASV_table_noS_save.tsv
#foreach PERCENT ( 1 0.99 0.97 0.95 0.90 0.85 0.80  )
foreach PERCENT ( 1  0.99 0.97 )
set FILE_NAME=${QUERY_IN}_trimmed_DNA.fasta 
  cd-hit-est -i ${MAIN_DIR}${QUERY_IN}_trimmed_DNA.fasta \
  -o ${MAIN_DIR}${QUERY_IN}_trimmed_DNA.fasta_${PERCENT} \
  -d 0 \
  -c $PERCENT \
  -T $THREADS_LIMIT
 if ($? != 0) then
  echo " command exit status after previous command is $?"
  exit
 endif
   
  #add the cluster name of different OTU % to the ASV table
  python $python_add_OTUs_to_table_script \
  ${MAIN_DIR}${R_dada_script_out_name}_final_ASV_table_noS.tsv \
  ${MAIN_DIR}${QUERY_IN}_trimmed_DNA.fasta_${PERCENT}.clstr\
  $PERCENT 
 if ($? != 0) then
  echo " command exit status after previous command is $?"
  exit
 endif
end


#takes an ASV table with the OTUs added and splits them
#this script will fail if the python script failed
echo "making the finalt ASV table (no singeltons, no duplicate sequences, attached clades) in R"
Rscript $R_script_OTU_tables \
${MAIN_DIR} \
${R_dada_script_out_name}_final_ASV_table_noS.tsv
if ($? != 0) then
 echo " command exit status after previous command is $?"
 exit
endif


echo "pipeline finished"

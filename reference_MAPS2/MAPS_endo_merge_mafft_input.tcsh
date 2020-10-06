#!/bin/tcsh
#PBS -q cdb
#PBS -o qsub_mafft_o.txt
#PBS -e qsub_mafft_e.txt
#PBS -N create_ref_tree_MP
#PBS -l mem=50gb
#PBS -l ncpus=8
#PBS -m ea
#PBS -M florian.prodinger@gmx.net




source /etc/profile.d/modules.sh
module load mafft/7.453
module load FastTree/2.1.11 
module load pplacer/1.1.alpha19
module load Python/3.7.5
module load raxml/7.2.8.PTHREADS

#module load raxml/8.2.12.AVX2.PTHREADS


set MAIN_DIR="/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/dada2_references/"
set MAPS_DIR="MAPS_references/"
set NEW_SEQS_TO_ADD="20200601_MAPS2_MIMI_TARA_EUK_BAC_ref_seqs.fasta" #"ncldv_6012.m700_10new_mimi.faa"
set MAPS_OLD_TREE="Pplacer.aln"

set OUT_FILE="20200601_MAPS2_MIMI_TARA_EUK_BAC_ref_seqs"


#merges two of endo's files without using seqs present in MAPS reference alignment
#cd $MAIN_DIR
#python merge_MAPS_reference_and_endo_tree.py

#echo "mafft fasta..."
#mafft --quiet --thread 8 --6merpair --addfragments \
# {$MAIN_DIR}{$NEW_SEQS_TO_ADD}  {$MAIN_DIR}{$MAPS_DIR}{$MAPS_OLD_TREE} > {$MAIN_DIR}{$OUT_FILE}_alingment.fasta
#
#exit
#
#echo "mafft phylip..."
#mafft --quiet --thread 8 --6merpair  --phylipout --namelength 75  --addfragments \
# {$MAIN_DIR}{$NEW_SEQS_TO_ADD}  {$MAIN_DIR}{$MAPS_DIR}{$MAPS_OLD_TREE} > {$MAIN_DIR}{$OUT_FILE}_alingment.phylip
#
#echo "FastTree"
#FastTree {$MAIN_DIR}{$OUT_FILE}_alingment.phylip > {$MAIN_DIR}{$OUT_FILE}_nwk.tre
#

echo "raxml"
cd $MAIN_DIR
raxmlHPC -s {$MAIN_DIR}{$OUT_FILE}_alingment.phylip \
 -T 24 \
 -p 1 \
 -m PROTGAMMAWAG \
 -g {$MAIN_DIR}{$OUT_FILE}_nwk.tre \
 -n {$OUT_FILE}_pplacer

# -s = input file
# -m = Protin matrix  -m PROTGAMMAWAG
# -g = tree 
# -n = output file

#mihara san's input for her tree in MAPS: 
#/user1/scl1/mihara/src/RAxML-7.2.7-ALPHA/raxmlHPC -s PolB_Refs_Cells_Viruses.phylip -n PolB_Refs -m PROTGAMMAWAG -g /user1/scl1/mihara/result/rpob/pi2/dpo/PolB_Refs_Cells_Viruses_rooted.nwk 





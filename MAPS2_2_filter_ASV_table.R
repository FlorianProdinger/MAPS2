#!/bin/R/3.6.1


#check user input
#read  out user input
print("[R] script for reading jplace and the blastx output") 

if (length(commandArgs(trailingOnly = T)) == 0){
 print("please enter a vaild directory after the R script name")
 quit()
} else if ( dir.exists( commandArgs(trailingOnly = T)[1] )){
 input_u  <- commandArgs(trailingOnly = T)
 out_dir <-input_u[1]
 #directory needs to be set to write tables
 setwd(out_dir)
 tree_file <- input_u[2]           #"test_20200526_5_samples_DNA_sequences.fasta_AA.fasta_pplacer.jplace"
 ASV_table_file <- input_u[3]      #"test_20200526_5_samples_ASV_table.tsv"
 fasta_output_file <- input_u[4]   #"test_20200526_5_samples_DNA_sequences.fasta_BLASTX_out.txt"
} else {
 print("[Rscript] directory not found")
 quit()}



require(ggtree, quietly = T)
require(treeio, quietly = T)
#tutorial to jplace trees in R
#https://bioc.ism.ac.jp/packages/3.2/bioc/vignettes/ggtree/inst/doc/treeImport.html#s4-classes

tree_MP <- read.jplace(tree_file)

#makes a dataframe with all the tree data (no placement)
tree_MP_df <-  fortify(tree_MP)



#define a clade which is considred Mimiviridae from two tips:
MIMI_node1 <- "MIMI_POV"
MIMI_node2 <- "MIMI_megavirus_bus"

MIMI_MRCA_node <- MRCA(get.tree(tree_MP), c( MIMI_node1, MIMI_node2))
#get all the "offspring" nodes of the most recent ancestor of Mimiviridae
MIMI_nodes <-  offspring(get.tree(tree_MP), MIMI_MRCA_node)


#returns a data frame with all the placed ASVs and where they were placed
tree_MP_placement_df <- as.data.frame(get.placements(tree_MP, by="best"))

#add a column to the dataframe showing if ASV was placed in Mimiviridae (bool)
tree_MP_placement_df$IS_MIMI <- tree_MP_placement_df$node %in% MIMI_nodes




########################

#load the OTU table from previous R script
ASV_tab <- read.table( ASV_table_file , header=1, row.names=1, sep = "\t" )
samples <- colnames( ASV_tab )



#was OTU placed by pplacer??
ASV_tab$pplaced_T_F <- rownames(ASV_tab) %in% tree_MP_placement_df$name



##############################
# best hit to none MIMI ASVs #
##############################
tab_fasta <- read.table(fasta_output_file, header=F, sep=",")
colnames(tab_fasta) <- c("ASVid", "AA_seq", "e_value", "percent_identity", "reading_frame", "best_hit", "other_hit")
not_MIMI_ASV <- tab_fasta$ASVid[ !grepl( "MEGA", tab_fasta$best_hit)]

ASV_tab_2 <- ASV_tab[ASV_tab$pplaced_T_F,]

##############################


###################################################
# check ASVs that were assigned to multiple nodes #
##################################################


asvID_to_number <-function( ASV_ID ){ as.numeric(gsub("ASV_","", ASV_ID))}
all_T <- function( T_F_vector){ length(T_F_vector)==sum(T_F_vector) }

MIMI_ASV_aggre <- aggregate(tree_MP_placement_df$IS_MIMI  ,  by=list( tree_MP_placement_df$name), FUN=all_T)
MIMI_ASV_aggre <- MIMI_ASV_aggre[ order(unlist(lapply( MIMI_ASV_aggre$Group.1, FUN=asvID_to_number))) , ]

#print( head(tree_MP_placement_df) )
#print( head( MIMI_ASV_aggre ))
#print( head( ASV_tab_2 ))
#print( sum( MIMI_ASV_aggre$Group.1 == rownames(ASV_tab_2) ))

if( !(nrow(ASV_tab_2)  ==  sum( MIMI_ASV_aggre$Group.1 == rownames(ASV_tab_2) ))){
print("ERROR with the ASV IDs")
quit()}

ASV_tab_2$IS_MIMI <- MIMI_ASV_aggre$x



ASV_tab_final <- ASV_tab_2[ ASV_tab_2$IS_MIMI , samples  ] 
#print( head( ASV_tab_final ))


#####################
# singelton removal #
#####################
ASV_tab_final_no_singelton <- ASV_tab_final[ colSums( ASV_tab_final  ) != 1 ]



#output the loss in reads, filter statistics
read_loss_df <- data.frame(  "dada2" = colSums(ASV_tab[,samples] ),
            "pplaced_ASVs" = colSums(ASV_tab_2[,samples] ),
            "MIMI_ASVs" =    colSums(ASV_tab_final[,samples] ))
read_loss_df$sample_name <- rownames(read_loss_df)
print("[R] created statistics file...")
print( read_loss_df ) 



#####################
# outputting tables #
#####################

ASV_table_file_new <- gsub( "ASV_table.tsv", "final_ASV_table_noS.tsv", ASV_table_file, ) 
write.table( file= ASV_table_file_new, ASV_tab_final_no_singelton, sep="\t", quote = F )

#########################
# outputting statistics #
#########################
 
previous_stat_file_name <- gsub("ASV_table.tsv", "import_merge_statititics_table.tsv", ASV_table_file)
ASV_table_stats_file <- gsub( "ASV_table.tsv", "filter_statistics.tsv", ASV_table_file )

#read stats of previous R script output (qual filter, merging, chimera checking)
filter_stats_1 <- read.table(previous_stat_file_name , header=1, row.names=1,sep="\t")

#combine the two stats tables:
ASV_table_stats <- cbind( filter_stats_1, read_loss_df)

write.table( file= ASV_table_stats_file , ASV_table_stats , sep="\t", quote = F )

save.image( file="dada2_pipeline_Rscript_tree")
#save.image( file=".dada2_pipeline_Rscript_tree")
print("[R] finished ASV_table filtering.")


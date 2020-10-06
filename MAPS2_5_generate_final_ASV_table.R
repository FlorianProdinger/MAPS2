#takes an 1) current dir and 2) OTU table that has clstr information in columns
#splits it in several tables

user_input <- commandArgs( trailingOnly=T)

if (length( user_input ) != 2) {
 print("[R] enter OTU table")
 print( length( user_input ))
 quit()
}


work_dir <- user_input[1] 
setwd( work_dir )

tab_file_name <-  user_input[2] # "step_1_dada2_derep_out_final_ASV_table_noS.tsv"

tab <- read.table( tab_file_name, header=1, row.names=1, sep="\t")

#tab with only reads
tab_no_clstr <- tab[ ,  !grepl("OTU", colnames(tab)) ] 

#at which percentages were OTUs clustered
OTU_per <- colnames(tab)[ grepl("OTU", colnames(tab)) ]

#make 1 table for each OTU percentage
for (per in OTU_per){
 print( paste0( "[R] working on ", per , " table") )
 #print(per)
 tab_test <- as.data.frame(aggregate( tab_no_clstr, FUN=sum, by =list( tab[ ,per  ] ) ))
 rownames( tab_test ) <- tab_test$Group.1
 tab_test$Group.1 <- NULL
 tab_test <- tab_test[ order( rowSums(tab_test), decreasing=T),]
 tab_test_noSingel <- tab_test[ rowSums(tab_test) >1,]
 print(paste0( "[R] ", sum( rowSums(tab_test) < 2 )," singeltons removed in: ",  per ))
 tab_file_name_out <- paste0( gsub(".tsv","_", tab_file_name) , per, ".tsv")
 print( paste0( "[R] writing ", per , " table") )
 write.table( file=tab_file_name_out, sep="\t", quote=F, tab_test_noSingel)
}

print("[R] generated several OTU tables.")


#This script uses dada to read raw fasta files
#and output an OTU table and a fasta file with cleaned reads
#=====
#use "module load Bioconductor/3.11" to run this script 
#tutorial here: https://benjjneb.github.io/dada2/tutorial.html
# most of the dada pipeline is from this page


if (length(commandArgs(trailingOnly = T)) != 4){ 
 print("[R] please enter: out dir, out file name,  the raw file directory and how many threads should be used")
 quit()
} else if ( dir.exists( commandArgs(trailingOnly = T)[1] )){
input_u <- commandArgs(trailingOnly = T)
 out_dir <- input_u[1]
 out_name <- input_u[2]  #"20200528_dada2_out"
 path_to_raw <- input_u[3]
 THREADS_LIMIT <- as.integer(input_u[4] )
} else {
 print("[R] directory not found")
 quit()}

#if user input is correct load moduls
library(dada2)
library(seqinr)
library(RcppParallel)

#make sure that tables are in the specified dir
print("[R] setting out directory:")
print( out_dir )
setwd(out_dir)

#path for raw files
path <-  path_to_raw
#<- "/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/test_qsub_20200527/0_RAW"


#limiting the CPU usage depending on user input
print(paste( "[R] using", THREADS_LIMIT , "threads"))
setThreadOptions(numThreads = THREADS_LIMIT )


 
#"/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/test"
#path <- "/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/test"
list.files( path )

##########################
### IMPORT FASTA FILES ###
##########################
#getting fastq files
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

#sample names 
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
out_dir_no_slash <- substring( out_dir  ,1,nchar( out_dir  )-1)

#make new directory for filtered sequences
filtFs <- file.path( out_dir_no_slash  , "filtered_dada2", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path( out_dir_no_slash  , "filtered_dada2", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

print("[R] performing trim and filter, saving files in this directory:")
print( filtFs[1] )


##############
### FILTER ###
##############
#default maxEE=Inf, tutorial: c(2,2)
#higher values gives more reads, hence I chose c(5,5)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(240,240),
                     maxN=0,
                     maxEE=c(2,2),
                     truncQ=2,
                     rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,
                     trimLeft=40, trimRight=38) 

print("[R] filtered and trimmed sequences: saving workspace")
save.image("R_script_dada_workspace")


##############
### ERRORs ###
##############

#learn error rates dor dada2 ASV generation
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
print("[R] learned error rates: saving workspace")
save.image("R_script_dada_workspace")


#visualize the error rates
out_name_pdf <- paste0( out_name , ".pdf")
pdf(out_name_pdf)
plotErrors(errF, nominalQ=TRUE)
dev.off()
print("[R] generated error rates plot:")
print( out_name_pdf)


#####################
### DEREPLICATION ###
#####################

#dereplicating sequences
print("[R] dereplicating sequences")
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
#Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names



#############
### DADA  ###
#############
#core sample interference 
print("[R] performing dada2" )
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
print("[R] performed dada2 ASV generation: saving workspace")
save.image("R_script_dada_workspace")



################
### MERGING  ###
################
#The mergers object is a list of data.frames from each sample. Each data.frame contains the merged $sequence, its $abundance, and the indices of the $forward and $reverse sequence variants that were merged. Paired reads that did not exactly overlap were removed by mergePairs, further reducing spurious output.
print("[R] merging reads")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, maxMismatch = 0, minOverlap = 12, )
print("[R] merged ASVs: saving workspace")
save.image("R_script_dada_workspace")



#######################
### CHIMERA FILTER  ###
#######################
#save some of the steps for later output
seqtab <- makeSequenceTable(mergers)
#dim(seqtab)
#table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
print("[R] chimera filter: saving workspace")
save.image("R_script_dada_workspace")

#dim(seqtab.nochim)
print("[R] chimera filter:")
print(sum(seqtab.nochim)/sum(seqtab)) #use this function to save reads per step


##################################
####  not from documentation  ####
#### these are my own scripts ####
##################################

#seqtab.nochim => rownames is DNA seq table is read count
#change the seqtab.nochim table to an ASV table and a fasta file
ASV_tab <- as.data.frame(seqtab.nochim)

# sites are rows columns are OTUs. Different from usual!!
# colnames(ASV_tab)[1]
#"ATTCTGTCTTCTATCGAAGAAACGATACGTTGGTATGTTATACGAAACGGACCCTAATAAATGTAAGCAAAAAAGCATGGGAATCGTTCTGAAACGACGTGATAATGCACCAATTGTAAAAGATATTTATGGGGGGATTATAGATATACTGATGAAGGAAGGGAATATTGTGAAAGCTGTGAACTTCTTACAAGACTGCTTACAAAAGATGATGGATGAACAATATCCTTTAGAAAAGCTCATTATCACAAAATCATTACGATCGAACTATAAGAACCCCAAGCAA"

#generating ASV_IDs and saving sequences
ASV_seq_df <- data.frame(
             "ASV_ID" = paste("ASV", c(1: length(colnames(ASV_tab))  ), sep="_"),
             "sequence" = colnames(ASV_tab)
             )


#replacing the column names with ASV_IDs
colnames(ASV_tab) <-ASV_seq_df$ASV_ID
out_ASV_tab <- t( ASV_tab ) #more natural orientation to me, rows are ASVs, columns are sites

#keep singeltons, might not be singeltons after trimming
#print("[R] singeltons per sample (not removed):")
#print( sum( rowSums(out_ASV_tab) == 1 ))

print("[R] script has finished, outputting tables and saving workspace")
save.image("R_script_dada_workspace")

#####################
### OUTPUT TABLES ###
###   AND FASTA   ###
#####################

write.table( file= paste(out_name, "ASV_table.tsv",sep="_"), out_ASV_tab, sep="\t", quote = F, row.names = T)

write.fasta( file.out= paste(out_name, "DNA_sequences.fasta",sep="_"),
             sequences= as.list(ASV_seq_df$sequence),
             names = ASV_seq_df$ASV_ID,
             as.string = T) 


##########################
### OUTPUT STATISTITCS ###
##########################
stat_file_name <- paste(out_name, "import_merge_statititics_table.tsv",sep="_")

stat_table <- cbind( out,
                     "merge" = (rowSums(unlist(makeSequenceTable(mergers)))),
                     "chimera_removed" = (rowSums(unlist( seqtab.nochim )))
                      )

write.table( file= stat_file_name, stat_table, sep="\t", quote = F, row.names = T)

print("[R] wrote table and fasta files. script finished succesfully.")


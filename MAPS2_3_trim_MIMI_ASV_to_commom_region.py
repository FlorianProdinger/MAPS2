import os
import sys
from Bio import SeqIO
import random
import re #for searching in strings

#use user input here
script_id = "python - trimming script"

user_input = sys.argv
if len( user_input ) != 6:
 print( "python - trimming script" )
 print("please input:\n\tfasta_seqs\n\tthe AA alingment file\n\tthe blast output filename\n\tcut left border\n\tcut right border.")
 print( len( user_input ))
 quit()
else:
 dada2_read_file_name = user_input[1] #"20200528_dada2_out_DNA_sequences.fasta"
 aln_file_name = user_input[2] #"20200528_dada2_out_DNA_sequences.fasta_AA.fasta_alingment_3.fasta"
 blast_out_file = user_input[3] #"20200528_dada2_out_DNA_sequences.fasta_BLASTX_out.txt"
 left_border = int(user_input[4])  #23350
 right_border= int(user_input[5]) #23970


#if files are not found, quit
for file_name in user_input[1:4]:
 if not os.path.isfile(file_name): # not in os.listdir():
  quit( f"{file_name} not found, quitting" )



#############
# functions #
#############

def first_AA( seq_string ):
 #returns place of first not "-" character
 return re.search( r"[^-]", str( seq_string  )).start() 

def last_AA( seq_string ):
 #returns place of last not "-" character
 str_len = len( str( seq_string ))
 return str_len - re.search( r"[^-]", str( seq_string  )[::-1] ).start()  

def get_polB_region( seq_string ):
 return str( seq_string )[ first_AA(seq_string) : last_AA(seq_string) ]


######################
# parse the DNA seqs #
######################

print("[python] parsing DNA fasta")
ASV_DNAseq_dic = {}
for record in SeqIO.parse( dada2_read_file_name  , "fasta"):
 ASV_DNAseq_dic[record.id] = str(record.seq) 

#####################
# parse the AA seqs #
#####################

print("[python] parsing AA fasta")
ASV_info_dic = {}
min_seq_len = 50
#cut_left_side = 0
#cut_right_side = 30000 #just a large placeholder number


#parse the alignment fasta file, ASV_aln_dic has only ASVs
for record in SeqIO.parse( aln_file_name  , "fasta"):
 if record.id.startswith("ASV") and len(str(record.seq).replace("-","")) > min_seq_len:
  ASV_info_dic[record.id] = [str(record.seq)]
  first_AA_aln = first_AA( record.seq )
  last_AA_aln  = last_AA( record.seq )
  ASV_info_dic[record.id].append( first_AA_aln  )
  ASV_info_dic[record.id].append( last_AA_aln   )

#adding the reading frame shift and the unsaligned AA seq to the information list
with open(blast_out_file) as blast_handle:
 for line in blast_handle:
  line_list = line.split(",")
  ASV_ = line_list[0]
  blast_seq = line_list[1] 
  reading_frame = line_list[4] 
  if ASV_ in ASV_info_dic and len(ASV_info_dic[ASV_]) != 5 :
   ASV_info_dic[ ASV_ ].append( reading_frame )
   ASV_info_dic[ ASV_ ].append( blast_seq ) 

#make a list of removed OTUs due to including stop codons
STOP_codon_ASV = []
for ASV_, info_list in ASV_info_dic.items(): 
 if "*" in info_list[4]:
  STOP_codon_ASV.append(ASV_)


#for ASV_, info_list in ASV_info_dic.items():
# seq_1 = str(info_list[0]).replace("-", "")
# seq_2 = info_list[4].replace("-","")
# if "*" in seq_2: continue
# if seq_1 != seq_2:
#  print(ASV_, info_list, "\n", seq_1, "\n", seq_2)


###########################
## Cutting the sequences ##
###########################

#>>> [input(f'{k},{v}') for k,v in ASV_info_dic.items()]
#looked at R plot to choose borders
# user input


#adding the cut sequence to the info_list and adding the cut DNA seq:
print("[python] cutting the DNA and AA sequence...")
ASV_cut_AA_DNA_dic = {}

for ASV_, info_list in ASV_info_dic.items():
 alinged_AA = str( info_list[0] )
 if "*" in alinged_AA: continue #don't include if stop codon
 frame_shift = int( info_list[3] )
 alinged_AA_noGap = alinged_AA.replace("-","")
 trimed_AA = alinged_AA[ left_border : right_border ] #here the AA sequence is cut
 trimed_AA_noGap = trimed_AA.replace("-","")
 if len(trimed_AA_noGap) == 0: continue #if the sequence was mapped to a not G domain region don't include
 cut_AA_left = alinged_AA_noGap.split( trimed_AA_noGap )[0] #how many AA were removed by cutting
 #cut_AA_left, cut_AA_right = alinged_AA_noGap.split( trimed_AA_noGap )
 #print(cut_AA_left, trimed_AA_noGap ,cut_AA_right)
 #print( len(cut_AA_left), len(trimed_AA_noGap) , len(cut_AA_right))
 DNA_seq = str( ASV_DNAseq_dic[ ASV_ ])
 #print( ASV_ in ASV_DNAseq_dic )
 cut_DNA_left = len(cut_AA_left)*3 + frame_shift-1 #cut the DNA sequence
 cut_DNA_right = cut_DNA_left + len(trimed_AA_noGap)*3 + frame_shift-1
 DNA_seq_trimmed = DNA_seq[ cut_DNA_left : cut_DNA_right ]
 ASV_cut_AA_DNA_dic[ ASV_ ] = [ DNA_seq_trimmed, trimed_AA_noGap]

 #ASV_cut_AA_DNA_dic[ ASV_ ] = {"DNA": DNA_seq_trimmed, "AA": trimed_AA_noGap}

########################
## write output fasta ##
########################
print("[python] outputting DNA and AA fasta with trimmed seqs")

file_name_out_AA  = dada2_read_file_name + "_trimmed_AA.fasta"
file_name_out_DNA = dada2_read_file_name + "_trimmed_DNA.fasta" 

with open( file_name_out_AA, "w" ) as AA_handle:
 with open( file_name_out_DNA, "w" ) as DNA_handle:
  for ASV_, DNA_AA_list in ASV_cut_AA_DNA_dic.items():
   DNA_seq = DNA_AA_list[0]
   AA_seq = DNA_AA_list[1]
   DNA_handle.write( f">{ASV_}\n{DNA_seq}\n")
   AA_handle.write( f">{ASV_}\n{AA_seq}\n")

#########
## END ##
#########
quit()







####################################################
            ### D E B U G I N G ### 
####################################################

quit()
for ASV_, info_list in ASV_info_dic.items():
  print(ASV_, info_list[1:])
  input()

#check the length ditribution of the polb gene files
len_seq_dic = {}
for info_list in ASV_info_dic.values():
 seq = info_list[0]
 len_polb = len(get_polB_region(seq).replace("-",""))
 if len_polb in len_seq_dic:
    len_seq_dic[ len_polb ] += 1
 else:
    len_seq_dic[ len_polb ] = 1

#save the length ditribution of the polb genes to plot in R (i dont like pandas :/ )
polb_border_dist_file = "polb_border_dist_file.txt"
with  open( polb_border_dist_file , "w") as border_out_handle:
 border_out_handle.write("ASVs,leftBorder,rightBorder\n")
 [border_out_handle.write( f"{ASV},{info[1]},{info[2]}\n") for ASV, info in ASV_info_dic.items()]


print("read lenth distribution (AA)\nlength,ASVs")
[print( f"{x},{len_seq_dic[x]}") for x in sorted(list(len_seq_dic.keys()))]

polb_len_dist_file = "polb_len_dist_file.txt"
with  open( polb_len_dist_file, "w") as len_out_handle:
 len_out_handle.write("length,ASVs\n")
 [len_out_handle.write( f"{x},{len_seq_dic[x]}\n") for x in sorted(list(len_seq_dic.keys()))]



quit()




#check length after trimming
trim_len_dic = {} 
for info_list in ASV_info_dic.values():
 seq = str(info_list[0])
 trim_seq = seq[ left_border : right_border ]
 trim_seq_len = len( trim_seq.replace("-","" ))
 #if trim_seq_len<40: continue 98% of all OTUs are longer than 40 AA 
 if trim_seq_len in trim_len_dic:
  trim_len_dic[ trim_seq_len ] += 1
 else:
  trim_len_dic[ trim_seq_len ] = 1


#make a practice dictionary with ~1000 entries 


ASV_info_dic_2 = {}
while len( ASV_info_dic_2 ) < 40:
 ASV = random.choice(list(  ASV_info_dic.keys() ))
 ASV_info_dic_2[ASV] = ASV_info_dic[ASV]


ASV_aln_dic_2 = {}
while len( ASV_aln_dic_2 ) < 100:
 ASV = random.choice(list(  ASV_aln_dic.keys() ))
 ASV_aln_dic_2[ASV] = ASV_aln_dic[ASV]


seqlen_list = [int( len(get_polB_region(seq).replace("-",""))  ) for seq in ASV_aln_dic.values()]
sorted( set(seqlen_list ))


#check the length ditribution of the polb gene files
len_seq_dic = {}
for seq in ASV_aln_dic.values():
 len_polb = len(get_polB_region(seq).replace("-",""))
 if len_polb in len_seq_dic:
    len_seq_dic[ len_polb ] += 1
 else:
    len_seq_dic[ len_polb ] = 1

polb_len_dist_file = "polb_len_dist_file.txt"
with  open( polb_len_dist_file, "w") as len_out_handle:
 [len_out_handle.write( f"{x},{len_seq_dic[x]}") for x in sorted(list(len_seq_dic.keys()))] 







#checking if functions work
#[print(first_AA(seq), last_AA(seq)) for seq in ASV_aln_dic_2.values()]
#[print( len(get_polB_region(seq).replace("-",""))  ) for seq in ASV_aln_dic_2.values()]




quit()
##### R plots #### 
#library(ggplot2)
#library(reshape2)
#
#tab_1 <- read.csv("polb_len_dist_file.txt")
#pdf("polb_len_dist_file.pdf")
#ggplot(tab_1, aes(x=length,y=ASVs)) + geom_col()
#dev.off()
#
#
#tab2 <- read.csv("polb_border_dist_file.txt")
#tab2_ <- melt(tab2 )
#> table(tab2$rightBorder)[table(tab2$rightBorder)>500]
#23999 24000 24197 24217 24220 24262 24710 24742 24963 24995 24996 
#  560  1000  1832   543 11626  1113  1055  2216   545  2432  9662 
#> table(tab2$leftBorder)[table(tab2$leftBorder)>500]
#23069 23178 23195 23206 23214 23251 23291 23294 23324 23345 
#  699 10512  1155  7489  3168   507  5446  1343   991  1118 
#
#pdf("polb_border_dist_file.pdf")
#ggplot(tab2_, aes(x=value, col=variable)) + geom_bar(stat="count") + xlim(13000, 29200)
#dev.off()

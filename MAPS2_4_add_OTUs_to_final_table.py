#this script takes an 
# 1) OTU table
# 2) a cd hit cltr file of the sequences of the OTU table
# 3) a name for the a new column of the output OUT table
# the script outputs the same OTU table it takes, but adds another column, in which the clstr of the OTU is shown. If the OTU is not present in the clstr file it will add "NA" insted of the clstr number 


import sys

if len(sys.argv) != 4:
 print(len( sys.argv ))
 quit("enter an:\n\tOTU table\n\ta clstr file\n\tnew col name.")
else:
 table_file_name = sys.argv[1]
 #"step_1_dada2_derep_out_final_ASV_table_noS.tsv"

 clstr_file_name = sys.argv[2]
 #"step_1_dada2_derep_out_DNA_sequences.fasta_trimmed_DNA.fasta_1.clstr"
  
 new_col_name = "OTU_" + str(sys.argv[3])
 #"OTU_100" 

 new_table_file_name = table_file_name # + new_col_name + ".tsv"

#this function returns a dictionary
def make_dic_from_file_name( file_name ):
 ASV_clstr_dic = {}
 with open(file_name) as clstr_handle:
  for line in clstr_handle:
   line = line.strip("\n")
   if line.startswith(">Cluster"):
    clstr = line.strip(">")
   else:
#    print(line)
    ASV_1 = line.split(">")[1]
    ASV_ = ASV_1.split("...")[0]
    if ASV_ not in ASV_clstr_dic:
     ASV_clstr_dic[ ASV_  ] =  clstr 
    else: 
     print(f"double ASV {ASV_}")
     break
 #the dictionary has ASVs as key and its OTU as value
 return ASV_clstr_dic

#checking dic
#[input( f"{k}: {v}") for k,v in clstr_ASV_dic.items()]


ASV_clst_dic = make_dic_from_file_name( clstr_file_name )


def make_new_table( table_file_name, new_col_name):
 table_line_list = []
 with open( table_file_name ) as tab_handle:
  for line in tab_handle:
   line = line.strip("\n")
   line_split = line.split("\t")
   ASV_ = line_split[0]
   if not line.startswith("ASV_"):
    row_append = new_col_name 
   elif ASV_ in ASV_clst_dic:
    row_append =  ASV_clst_dic[line_split[0]]
   else:
    row_append = "NA" 
    #print( f"{ASV_} not found :/")
   line_split.append( row_append  )
   line_2 =  "\t".join(line_split)
   #input(line_2)
   table_line_list.append( line_2 )
 return table_line_list


file_out_lines = make_new_table( table_file_name, new_col_name)



with open(new_table_file_name, "w") as out_tab_handle:
 for line in file_out_lines:
  out_tab_handle.write(line)
  out_tab_handle.write("\n")







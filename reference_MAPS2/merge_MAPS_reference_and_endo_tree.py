#reads out some fasta files and writes unique entries in an output file

MAPS_ali_file = "Pplacer.aln"
endo_tara_file = "ncldv_6012.m700.faa"
endo_ref_file = "PolB_Refs_NCLDVs20200102.faa"

endo_path = "/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/dada2_references/endo_references"
MAPS_path = "/lustre1/aptmp/florian/megaviridae/all_experiments/Uranouchi_MP10/20191021_43_samples/MP_dada2/dada2_references/MAPS_references"

from Bio import SeqIO
import os

id_seq_dic = {}
MAPS_merge_list = ["PITH_pithovirus", "POX_Amsacta_moorei_entomopoxvirus_L_", "PHYCO_Bathycoccus_sp_RCC1105_virus_BpV1", "PHYCO_Micromonas_sp_RCC1109_virus_MpV1"]



for record in SeqIO.parse( os.path.join( MAPS_path , MAPS_ali_file)  , "fasta") :
 MAPS_merge_list.append( record.id )
  #print (record.seq)

for record in SeqIO.parse( os.path.join(endo_path, endo_ref_file)  , "fasta") :
  if record.id in MAPS_merge_list:
   continue
  else:
   id_seq_dic[ record.id ] = record.seq

  
for record in SeqIO.parse( os.path.join(endo_path, endo_tara_file )  , "fasta") :
  if record.id in MAPS_merge_list:
   continue
  else:
   id_seq_dic[ record.id ] = record.seq

out_file = "20200601_MAPS2_MIMI_TARA_EUK_BAC_ref_seqs.fasta"

with open( out_file, "w") as out_handle:
 for ID, seq in id_seq_dic.items():
  out_handle.write( f">{ID}\n")
  out_handle.write( f"{seq}\n")
  



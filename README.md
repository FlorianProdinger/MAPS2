# Welcome to MAPS2!

MAPS is short for "*M*EGAPRIMER *A*mplicon *P*rocessing *S*ystem" and is a analysis pipeline for Mimiviridae *polB* gene amplicon analysis. 

The first version (i.e. MAPS) was published by Li & Hingamp et al. (2018) and was mostly based on bash and python. The current version was built on R, bash and python and incorporates dada2 (Callahan et al. 2016).

The main script "MAPS2.tcsh" is the MAPS2 pipeline.
It calls in succession these steps:

## Strucutre of the pipeline
 1. The script called at first uses dada2 in R. This script does most of the work.
 2. blastx on the ASVs to check if they are viral sequences. (saves translated amino acid sequences).
 3. a bash command (sed...) to make a fasta file from the blast output.
 4. mafft adds *Mimiviridae* ASVs sequences to a reference file.
 5. pplacer places the ASVs in a reference tree (Endo et al. 2020).
 6. Another R script removes ASV that were not placed within the Mimiviridae branch (saves a statistics table).
 7. All ASVs are trimmed in a common region (common region hardcoded in this pipeline).
 8. Clustering trimmed ASVs at 100% (99%, 97%) nucleotide indentity
 9. Creating a final ASV table from clustered ASVs.


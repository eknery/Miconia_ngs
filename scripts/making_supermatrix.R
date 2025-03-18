### library
library("seqinr")
### read fastas
fasta1 = read.fasta("1_target_data/4_filtered_sequences/all_loci_filtered.fas")  
fasta2 = read.fasta("1_skimm_data/4_filtered_sequences/all_loci_filtered.fas")  
### species in each aligment
spp1 = names(fasta1)
spp2 = names(fasta2)
### length of each aligment
length1 = length(fasta1[[1]])
length2 = length(fasta2[[2]])
### insert empty lines for absent species
for(sp_name in spp1){
  boll = sp_name %in% spp2
  if(boll == FALSE){
    fasta2[[sp_name]] = rep("-", length2)
  }
}
### ordering species
fasta2 = fasta2[spp1]
### super matrix
smatrix = fasta1
for(sp_name in spp1){
  smatrix[[sp_name]] = c(smatrix[[sp_name]], fasta2[[sp_name]])
}
### output directory 
dir_out = "2_supermatrix/full_matrix/"
### exporting sequences 
write.fasta(
  sequences = smatrix, 
  as.string = F, 
  names = spp1,
  file.out = paste0(dir_out, "supermatrix.fasta"),
  nbchar = 100
)

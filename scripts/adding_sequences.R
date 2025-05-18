### library
library("seqinr")
### read fastas
fasta1 = read.fasta("1_ITS_data/ITS_merged.fas")  
fasta2 = read.fasta("1_skimm_data/4_filtered_sequences/all_loci_filtered.fas")
### species in each aligment
spp1 = names(fasta1)
spp2 = names(fasta2)
### length of each aligment
length1 = length(fasta1[[1]])
length2 = length(fasta2[[2]])
### more complete aligment
add_alignment = fasta1
### size of less complete aligment
n_patterns = length(fasta2[[1]])
### adding novel data and missing data
for(sp_name in spp1){
  boll = sp_name %in% spp2
  if(boll == TRUE){
    add_alignment[[sp_name]] = c(add_alignment[[sp_name]], fasta2[[sp_name]] )
  }
  if(boll == FALSE){
    add_alignment[[sp_name]] = c(add_alignment[[sp_name]], rep("-", n_patterns) )
  }
}

### export
dir_out = "2_supermatrix/"
write.fasta(
  sequences = add_alignment, 
  as.string = F, 
  names = names(add_alignment),
  file.out = paste0(dir_out, "ITS_plastid.fas"),
  nbchar = 100
)

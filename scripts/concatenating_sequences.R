### library
library("seqinr")
### all species names
spp_names = read.table("spp_names.txt", h=F)
spp_names = spp_names$V1
### choose data type
dtype = "1_skimm_data/"
### all species names
n_spp = read.table(paste0(dtype, "n_spp.txt"), h=F)
n_spp = n_spp$V1
### list file names
loci_names = list.files(path = paste0(dtype, "2_completed_aligned_sequences"), pattern = ".FNA")
### loci and number of sequences
loci_nspp = as.data.frame(cbind(loci_names,  c(44,n_spp)))
### pick the first aligment
concatenation = read.fasta(paste0(dtype,"2_completed_aligned_sequences/", loci_names[1]))
### minimum number of species to consider a locus
min_nspp = 42
### loop
for(i in 2:length(loci_names) ){
  ### pick other locus
  other_locus = read.fasta(paste0(dtype,"2_completed_aligned_sequences/", 
                                  loci_names[i]))
  ### retrieve number of species with data
  nspp_data = as.integer(loci_nspp$V2[loci_nspp$loci_names == loci_names[i]])
  ### check if aligment has enough data
  if(nspp_data >= min_nspp){
    ### add sequences from other locus
    for(sp_name in spp_names){
      concatenation[[sp_name]] = c(concatenation[[sp_name]], other_locus[[sp_name]])
    }
  }
}

### export
dir_out = "3_concatenated_sequences/"
write.fasta(
  sequences = concatenation, 
  as.string = F, 
  names = spp_names,
  file.out = paste0(dtype, dir_out, "all_loci.fasta"),
  nbchar = 100
)

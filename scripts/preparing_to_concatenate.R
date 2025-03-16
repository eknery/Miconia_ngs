library("seqinr")

### all species names
all_names = read.table("all_names.txt", h=F)
all_names = all_names$V1

### list file names
locus_names = list.files(path = "1_aligned_sequences", pattern = ".FNA")

### number of species in each aligment
n_spp = c()

### pick one fasta aligment
one_locus = read.fasta(paste0("1_aligned_sequences/", locus_names[1]))
### update n_spp
n_spp = c(n_spp, length(one_locus) )
### including missing species
n_patterns = length(one_locus[[1]])
for(name in all_names){
  boll = name %in% names(one_locus)
  if(boll == FALSE){
    one_locus[[name]] = as.integer(rep("-", n_patterns))
  }
}
### ordering sepcies
one_locus = one_locus[all_names]

### exporting sequences 
dir_out = "3_final_sequences/"
for(i in 1:length(my_fasta_list)){
  fasta_name = paste0(names(my_fasta_list)[i], ".fasta")
  write.phyDat(x = my_fasta_list[[i]], 
               file = paste0(dir_out,fasta_name), 
               format = "fasta", 
               colsep = "", 
               nbcol =100
  )
}


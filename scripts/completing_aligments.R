### library
library("seqinr")
### all species names
spp_names = read.table("spp_names.txt", h=F)
spp_names = spp_names$V1
### choose data type
dtype = "1_target_data/"
### list file names
loci_names = list.files(path = paste0(dtype,"1_aligned_sequences"), pattern = ".FNA")
### number of species in each alignment
n_spp = c()
### loop 
for(i in 1:length(loci_names)){
  tryCatch(
    {
      ### pick one fasta aligment
      one_locus = read.fasta(paste0(dtype,"1_aligned_sequences/", loci_names[i]))  
      ### update n_spp
      n_spp = c(n_spp, length(one_locus) )
      ### including missing species
      n_patterns = length(one_locus[[1]])
      for(name in spp_names){
        boll = name %in% names(one_locus)
        if(boll == FALSE){
          one_locus[[name]] = rep("-", n_patterns)
        }
      }
      ### ordering species
      one_locus = one_locus[spp_names]
      ### output directory 
      dir_out = "2_completed_aligned_sequences/"
      ### exporting sequences 
      write.fasta(
        sequences = one_locus, 
        as.string = F, 
        names = spp_names,
        file.out = paste0(dtype, dir_out, loci_names[i]),
        nbchar = 100
      )
    },
    error = function(e) {
      n_spp = c(n_spp, 0)
      cat("Skipping: ", loci_names[i],"; ")
      return(NULL)  # Return NULL to indicate failure
    }
  )
}

write.table(n_spp,paste0(dtype,"n_spp.txt"), row.names = F)

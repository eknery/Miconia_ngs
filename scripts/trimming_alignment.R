### library
library("seqinr")
library("ape")

### choose data type
dtype = "1_target_data/"

### chose alignmnet directory
dir_align = "3_concatenated_sequences/"

### load alignment
align = read.fasta(paste0(dtype, dir_align, "201_loci.fasta"))

### species names e number of sites
spp_names = names(align)
n_sites = length(align[[1]])

### converting to a matrix
mtx_align = matrix(unlist(align), ncol = n_sites, byrow = T)

### trimming
trim_align = del.colgapsonly(x = mtx_align, 
                                 threshold = 0.1,
                                 freq.only = FALSE)

### convert back to matrix
mtx = as.matrix(as.character(trim_align))

### convert back to list
list_trim = list()
for(i in 1:nrow(mtx)){
  list_trim[[i]] = paste0(mtx[i,], collapse = "")
}
names(list_trim) = spp_names

### export
write.fasta(
  sequences = list_trim, 
  as.string = T, 
  names = spp_names,
  file.out = paste0(dtype, "4_trimmed_sequences", "trim_201_loci.fasta"),
  nbchar = 100
)


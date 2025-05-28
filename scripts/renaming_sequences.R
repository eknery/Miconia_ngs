### library
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### choose input directory
dir_input = "1_aligned_sequences/"

### list FASTA names
loci_names = list.files(path = paste0(dir_input), pattern = ".FNA")

### trimming loci in loop
for(i in 1:length(loci_names) ){
  ### name of one locus
  locus_name = loci_names[i]
  tryCatch({
    ### load locus
    one_locus = read.fasta(paste0(dir_input, locus_name))
    ### remove old 'vatlheri'
    if("valtheri" %in% names(one_locus)){
      ### select species
      one_locus = one_locus[!names(one_locus) %in% "valtheri"]
    }
    ### name new 'valtheri'
    if("buddlejoides" %in% names(one_locus)){
      names(one_locus)[names(one_locus) == "buddlejoides"] = "valtheri"
    }
    ### export
    write.fasta(
      sequences = one_locus, 
      as.string = F, 
      names = names(one_locus),
      file.out = paste0("1_new_aligned_sequences/", locus_name),
      nbchar = 100
    )
  },
  error = function(e) {
    print(paste0("Could not select: ", locus_name))
    return(NULL)  # Return NULL to indicate failure
  })
} 

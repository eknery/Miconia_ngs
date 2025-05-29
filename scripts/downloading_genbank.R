if(!require("ape")) install.packages("ape"); library("ape")

### input directory
dir_input = "0_sanger_data/"
### read acessions table
access = read.csv(paste0(dir_input, "0_sanger_acessions.csv") )

### loci names
loci_names = colnames(access)[colnames(access) != "species"]

for(locus_name in loci_names){
  ### lines with information
  access_clean = access[!is.na(access[,locus_name]),]
  ### accession numbers and names
  access_num = access_clean[,locus_name]
  access_names = access_clean$species
  ### downloading sequences
  one_locus = read.GenBank( access.nb = access_num )
  ### naming sequences
  names(one_locus) = access_names
  ### exrpoting directory
  dir_out = dir_input
  ### export
  write.dna(
    one_locus, 
    file = paste0(dir_out, locus_name, "_raw.fas"), 
    format = 'fasta' 
  )
  ### check
  print(paste0("Download complete: ", locus_name))
}


if(!require("ape")) install.packages("ape"); library("ape")

### input directory
dir_input = "0_sanger_data/"
### read acessions table
access = read.csv(paste0(dir_input, "sanger_acessions.csv") )

### pick one locus
locus_name = "ETS"

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
  file = paste0(dir_out, locus_name, "_sanger.FNA"), 
  format = 'fasta' 
  )


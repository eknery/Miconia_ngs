### library
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### choose directory with FASTA data
dir_input = "0_sanger_data/"

### data from NGS - inn
skimm = read.fasta(paste0(dir_input, "ITS_skimm.FNA"))

### data from sanger - out
sanger = read.fasta(paste0(dir_input, "ITS_sanger.FNA"))

### get skimm data for species without sanger data
skimm_plus = skimm[!names(skimm) %in% names(sanger)]

### merging sanger with complementary skimm data
merged = sanger
for(i in 1:length(skimm_plus)){
  sp_name = names(skimm_plus)[i]
  merged[[sp_name]] = skimm_plus[[sp_name]] 
}

### exporting directory
dir_out = dir_input 

### export
write.fasta(
  sequences = merged, 
  as.string = F, 
  names = names(merged),
  file.out = paste0(dir_out, "ITS_merged.fas"),
  nbchar = 100
)

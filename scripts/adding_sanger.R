### library
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### choose data type
dtype = "1_target_data/"

### choose directory with FASTA data
dir_input = "3_selected_sequences/"

### list file names
loci_names = list.files(path = paste0(dtype,dir_input), 
                        pattern = ".FNA")

### choose sanger data
outtype = "1_ITS_data/"

### data from NGS - inn
inn_locus = read.fasta(paste0(outtype, "ITS_skimm.FNA"))

### data from sanger - out
out_locus = read.fasta(paste0(outtype, "ITS_sanger.FNA"))

### getting all species names across loci
all_spp_names = c()
for(i in 1:length(loci_names)){
  locus_name = loci_names[i]
  one_locus = read.fasta(paste0(dtype, dir_input, locus_name))  
  spp_names = names(one_locus)
  all_spp_names = sort(unique(c(all_spp_names, spp_names)))
}

### keep data from sampled species
sampled = inn_locus[names(inn_locus) %in% all_spp_names]

### get data from unsampled species
unsampled = out_locus[!names(out_locus) %in% names(sampled) ]

### merging sampled and unsampled
merged = sampled
for(i in 1:length(unsampled)){
  sp_name = names(unsampled)[i]
  merged[[sp_name]] = unsampled[[sp_name]] 
}

### remove species
remove_spp = c(
  "amoena",
  "budlejoides",
  "burchellii",
  "castaneiflora",
  "cinerascens",
  "collatata",
  "dorsaliporosa",
  "elata",
  "eriodonta",
  "fallax",
  "hyemalis",
  "lanata",
  "macuxi",
  "mayarae",
  "pennipilis",
  "petroniana",
  "ruficalyx",
  "ruschiana",
  "schwackei",
  "sclerophylla",
  "triplinervis",
  "valtheri"
)

merged = merged[!names(merged) %in% remove_spp]

### export
write.fasta(
  sequences = merged, 
  as.string = F, 
  names = names(merged),
  file.out = paste0(dtype, "3_selected_sequences/", "ITS.FNA"),
  nbchar = 100
)

### library
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### choose data type
dtype = "1_target_data/"

### chose directory with FASTA data
dir_input = "2_trimmed_sequences/"

### list FASTA names
loci_names = list.files(path = paste0(dtype, dir_input), pattern = ".FNA")

### species to remove
spp_remove = c(
  "albicans",
  "buddlejoides",
  "burchellii",
  "castaneiflora",
  "collatata",
  "cyathanthera",
  "dorsaliporosa",
  "elata",
  "eriodonta",
  "fallax",
  "hyemalis",
  "lanata",
  "pennipilis",
  "petroniana",
  "ruficalyx",
  "schwackei",
  "sclerophylla",
  "triplinervis",
  "valtheri"
)

### minimum number of species to maintain a locus
min_nspp = 60

### trimming loci in loop
for(i in 1:length(loci_names) ){
  ### name of one locus
  locus_name = loci_names[i]
  ### load locus
  one_locus = read.fasta(paste0(dtype, dir_input, locus_name))
  ### select species
  slc_locus = one_locus[!names(one_locus) %in% spp_remove]
  ### number of remaining species
  nspp = length(slc_locus)
  ### export if enough species remained
  if(nspp >= min_nspp){
    ### export
    write.fasta(
      sequences = slc_locus, 
      as.string = F, 
      names = names(slc_locus),
      file.out = paste0(dtype, "3_selected_sequences/", locus_name),
      nbchar = 100
    )
    ### check!
    print(paste0("Selection done: ", locus_name)) 
  } else {
    print(paste0("Discarding: ", locus_name))
  }
}

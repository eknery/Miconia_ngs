### library
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### choose input directory
dir_input = "1_aligned_sequences/"

### list FASTA names
loci_names = list.files(path = paste0(dir_input), pattern = ".FNA")

### species to remove
spp_remove = c(
  "albicans",
  "castaneiflora",
  "collatata",
  "cyathanthera",
  "dorsaliporosa",
  "elata",
  "eriodonta",
  "fallax",
  "hyemalis",
  "lanata",
  "longispicata",
  "pennipilis",
  "ruficalyx",
  "sclerophylla",
  "triplinervis"
)

### minimum number of species to maintain a locus
min_nspp = 40

### trimming loci in loop
for(i in 1:length(loci_names) ){
  ### name of one locus
  locus_name = loci_names[i]
  tryCatch({
    ### load locus
    one_locus = read.fasta(paste0(dir_input, locus_name))
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
        file.out = paste0("2_selected_sequences/", locus_name),
        nbchar = 100
      )
      ### check!
      print(paste0("Selection done: ", locus_name)) 
    } else {
      print(paste0("Discarding: ", locus_name))
    }
  },
  error = function(e) {
    print(paste0("Could not select: ", locus_name))
    return(NULL)  # Return NULL to indicate failure
  })
}
############################## SELECTING SANGER ###############################

## choose input directory
dir_input = "0_sanger_data/"

### list FASTA names
loci_names = list.files(path = paste0(dir_input), pattern = "_aligned.fas")

for(locus_name in loci_names){
  ### list FASTA names
  sanger = read.fasta(paste0(dir_input, locus_name))
  ### species to remove
  spp_remove = c(
    "collatata",
    "dorsaliporosa",
    "elata",
    "hyemalis",
    "lanata",
    "longispicata",
    "pennipilis",
    "triplinervis"
  )
  ### remvoving species
  slc_sanger = sanger[!names(sanger) %in% spp_remove]
  ### clean name
  clean_name = gsub("_aligned.fas","", locus_name)
  ### export
  write.fasta(
    sequences = slc_sanger, 
    as.string = F, 
    names = names(slc_sanger),
    file.out = paste0("2_selected_sequences/", clean_name, ".FNA"),
    nbchar = 100
  )
  ### check!
  print(paste0("Selection done: ", locus_name)) 
}



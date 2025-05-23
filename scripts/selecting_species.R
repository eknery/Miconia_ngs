### library
library("seqinr")
library("ape")

### choose data type
dtype = "1_target_data/"

### chose directory with FASTA data
fdir = "2_trimmed_sequences/"

### list FASTA names
loci_names = list.files(path = paste0(dtype, fdir), pattern = ".FNA")

### species to remove
spp_remove = c(
  "albicans",
  "buddlejoides",
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

### trimming loci in loop
for(i in 1:length(loci_names) ){
  ### name of one locus
  locus_name = loci_names[i]
  tryCatch({
    ### load alignment
    one_locus = read.fasta(paste0(dtype, fdir, locus_name))
    ### select species
    slc_locus = one_locus[!names(one_locus) %in% spp_remove]
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
  },
  error = function(e) {
    print(paste0("Skipping: ", locus_name))
    return(NULL)  # Return NULL to indicate failure
  })
}

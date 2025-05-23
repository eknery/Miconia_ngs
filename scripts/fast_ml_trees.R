### load libraries
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("stringr")) install.packages("stringr"); library("stringr")

### choose data type
dtype = "1_target_data/"

### chose directory with FASTA data
dir_input = "3_selected_sequences/"

### select output dir
dir_out = "4_fast_ml_trees/"

### list FASTA names
loci_names = list.files(path = paste0(dtype, dir_input), pattern = ".FNA")

### infer ML tree per locus
for(i in 1:length(loci_names)){ 
  ### name of one locus
  locus_name = loci_names[i]
  ### load alignment
  one_locus = read.phyDat(paste0(dtype, dir_input,  locus_name),
              format = "fasta",
              type = "DNA"
  )
  ### initial NJ tree
  nj_tree = NJ(dist.ml(one_locus))
  ### optimized ML tree
  ml_opt = optim.pml(
    pml(tree = nj_tree, data = one_locus), 
    model = "GTR"
  )
  ### ML tree name
  tree_name = str_remove(string = locus_name, pattern = ".FNA")
  ### export ML tree
  write.tree(phy = ml_opt$tree,
             file = paste0(dtype, dir_out, tree_name,".tree" ))
  ### check
  print(paste0("ML tree done: ", locus_name))
}

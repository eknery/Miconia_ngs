### load libraries
if(!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if(!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("TreeTools")) install.packages("TreeTools"); library("TreeTools")

### choose data type
dtype = "1_target_data/"

### chose directory with tree data
dir_input = "4_fast_ml_trees/"

### file names
tree_names = list.files(paste0(dtype, dir_input) )

### loading data
tree_list = list()
for(i in 1:length(tree_names) ){
  tree_name = tree_names[i]
  tree_list[[i]] = read.tree(file = paste0(dtype,dir_input, tree_name))
  names(tree_list)[i] =  str_remove(string = tree_name, pattern = ".tree")
}

################################# FIND COMMON SPECIES ###########################

### getting species per locus
all_names = c()
for(i in 1:length(tree_list)){
  some_names = tree_list[[i]]$tip.label
  all_names = sort(unique(c(all_names, some_names)))
}

### get names per locus
names_loci =  all_names
for(i in 1:length(tree_list)){
  boll_names = all_names %in% tree_list[[i]]$tip.label
  names_loci = cbind(names_loci, boll_names)
}

### transform to tibble
names_loci = as_tibble(names_loci)
### get species with all loci sequenced 
common_names = names_loci %>% 
  filter_at(vars(-names_loci), all_vars(. == TRUE) ) %>% 
  select(names_loci) %>% 
  pull()

################################ PROCESSING TREES #############################

# pruning trees to sampled species
pruned_trees_list = tree_list
for (i in 1:length(tree_list) ){
  pruned_tree = keep.tip(phy = tree_list[[i]], 
                        tip = common_names)
  pruned_trees_list[[i]] = pruned_tree
}

### transforming in vector
pruned_trees_vec = c()
for(i in 1:length(pruned_trees_list)){
  pruned_trees_vec = c(pruned_trees_vec, pruned_trees_list[[i]] )
}

################################ CONGRUENCE ANALYSES ##########################

### conver to multiphylo
prunned_trees = as.multiPhylo( pruned_trees_list)

### distance
dist = RF.dist(
  tree1 = prunned_trees,
  tree2 = NULL, 
  normalize = T
)

### PCOA
pcoa = pcoa(dist, correction="none", rn=NULL)
pcoa_df = as.data.frame(cbind(pcoa$vectors))
### get % var 
pc_rel_var = pcoa$values$Relative_eig
### axis names
pc_axis_1 = paste0("PCoA (", round(pc_rel_var[1]*100, 2), "%)" )
pc_axis_2 = paste0("PCoA (", round(pc_rel_var[2]*100, 2), "%)" )

### plot pcoa
pcoa_plot = ggplot(data = pcoa_df,
       aes(x=as.numeric(Axis.1),
           y=as.numeric(Axis.2)
           )) +
  
  geom_point(size = 2, alpha = 0.5) +
  labs(x= pc_axis_1, y= pc_axis_2)+
  
  guides(color = guide_legend(title="",
                             ncol = 5,
                             byrow = TRUE) ) + 
  
  theme(panel.background=element_rect(fill="white"),
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"),
        axis.title=element_text(size=12,face="bold"),
        legend.position = "bottom"
        )

### export plot
tiff(paste0(dtype, "pcoa_ml_trees.tiff"), 
     units="cm", width=10, height=9, res=600)
 pcoa_plot
dev.off()

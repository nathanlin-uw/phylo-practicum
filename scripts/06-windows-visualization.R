# Load everything in
library(ape)
library(phangorn)
library(phytools)
library(ggplot2)

getwd() #Check the working directory. we want to be in the results/RAxML/10Mb-concatenation folder
# setwd("phylo-practicum/results/RAxML/10Mb-concatenation/") 


# For author-provided concatenation files
# List all .bestTree files. $ ensures the end of the name
tree_files <-list.files(pattern="\\.raxml.bestTree$") 

# Set up list for all trees
trees <- list() # list with all the trees
class(trees)<- "multiPhylo" #make it a multiphylo object for ease of use with other 

# Read each tree in
i<-1
for(tree_file in tree_files){ ##go thru each file and read the tree
  trees[[i]]<- read.tree(tree_file)
  i<-i+1
}

# Re-reroot all our gene trees by the respective outgroup
for(i in 1:length(trees)){
  trees[[i]]<- root(trees[[i]],
                         outgroup = "H_vulgare_HVens23",
                         resolve.root=TRUE)
  trees[[i]]<-chronos(trees[[i]]) ## make ultrametric for nicer densitree
}

# Create new parsimony consensus supertree
st<-superTree(trees)
st<-root(st,"H_vulgare_HVens23",resolve.root = T)
plot(st)

# Plot densitree 
densiTree(trees,consensus=st,scaleX=T,type='cladogram', alpha=0.1)

### For our own sliding windows concatenation
trees2 <- read.nexus("../../../data/Wheat_Relative_History_Data_Glemin_et_al/Densitree_OneCopyGenes-modified.nex")

#re-reroot all our gene trees by the respective outgroup
for(i in 1:length(trees2)){
  trees2[[i]]<- root(trees2[[i]],
                         outgroup = "H_vulgare_HVens23",
                         resolve.root=TRUE)
  trees2[[i]]<-chronos(trees2[[i]]) ## make ultrametric for nicer densitree
}

st2<-superTree(trees2)
st2<-root(st2,"H_vulgare_HVens23",resolve.root = T)
plot(st2)

densiTree(trees2,consensus=st2,scaleX=T,type='cladogram', alpha=0.1)
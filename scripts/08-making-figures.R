# Loading everything in
library(ape)
library(phangorn)

# BiocManager::install("ggtree")
library(ggtree)
library(dplyr)

# Set the working directory to the results/ folder
# setwd("//wsl.localhost/Ubuntu/home/ntlin/phylo-practicum/results")

###### FIGURE 1A ######
# Read the full concatenation tree and supertree in
full_concat_tree = read.tree(file="RAxML/full_concatenation/triticeae_allindividuals_OneCopyGenes.fasta.raxml.bestTree")
supertree = read.tree(file="07-supertree.tre")

# Root the trees
full_concat_tree = root(full_concat_tree, outgroup="H_vulgare_HVens23", resolve.root=TRUE)
supertree = root(supertree, outgroup="H_vulgare_HVens23", resolve.root=TRUE)

# Suppose your trees are called full_concat_tree and supertree
# First, ladderize both trees
full_concat_tree <- ladderize(full_concat_tree)
supertree <- ladderize(supertree)

# Make sure they have the same tip labels (they should both have 47 tips)
common_tips <- intersect(full_concat_tree$tip.label, supertree$tip.label)

length(common_tips)
length(full_concat_tree$tip.label)
length(supertree$tip.label)

# Reorder the second tree to match the tip order of the first
supertree <- reorder.phylo(supertree, order = "cladewise")
supertree <- rotateConstr(supertree, full_concat_tree$tip.label) # rotate to match full_concat_tree tip order

# Plot side by side
par(mfrow = c(1, 2))
plot(full_concat_tree, main = "Full concatenation", cex = 0.8)
plot(supertree, main = "Supertree", cex = 0.8)

# Calculate the RF distance to verify that they really are different trees
RF.dist(full_concat_tree, supertree)

# Reproduce edge colors on the full concatenation tree
species_colors <- c(
  "Ae_umbellulata" = "yellow",
  "Ae_caudata" = "orange",
  "Ae_comosa" = "darkorange3",
  "Ae_uniaristata" = "sienna4",
  "Ae_bicornis" = "purple",
  "Ae_longissima" = "pink",
  "Ae_sharonensis" = "mediumpurple1",
  "Ae_searsii" = "maroon2",
  "Ae_tauschii"    = "red",
  "T_boeoticum" = "darkgreen",
  "T_urartu" = "green",
  "Ae_speltoides" = "blue4",
  "Ae_mutica" = "steelblue1"
)

# Need to grab the species from the tip labels
  # Our species names are Genusabbreviation_specificepithet
  # Our tip names are Genusabbreviation_specificepithet_everythingelse
species_from_tip_name <- function(tip_name) {
  # Split by underscores, it gives us a list so we need to grab the first and only element
  split_tip_name <- strsplit(tip_name, "_")[[1]]
  # Put the first two parts back together
  paste(split_tip_name[1:2], collapse="_")
} 

# Apply our previous function to all the tip labels
tip_species <- as.factor(sapply(full_concat_tree$tip.label, species_from_tip_name))

p <- ggtree(full_concat_tree)

# Loop over species to color clades
for(species in names(species_colors)) {
  # Get tips belonging to this species
  tips <- full_concat_tree$tip.label[tip_species == species]
  # Get MRCA node
  node <- getMRCA(full_concat_tree, tips)
  if(!is.null(node)) {
    p <- p + geom_hilight(node = node, fill = species_colors[species], alpha = 0.3)
  }
}

# Add tip labels
p <- p + geom_tiplab()
p

# Manually rotate clades to mimic figure 1A
p + geom_text2(aes(subset = !isTip, label = node), hjust = -0.3)

# Rotate the clades
p <- ggtree(full_concat_tree)

p2 <- rotate(p, node = 55)
p2 <- rotate(p2, node = 62)
p2 <- rotate(p2, node = 50)
p2 <- rotate(p2, node = 49)
p2 <- rotate(p2, node = 69)
p2 <- rotate(p2, node = 83)
p2 <- rotate(p2, node = 88)
p2 <- rotate(p2, node = 71)
p2 <- rotate(p2, node = 74)


for(sp in names(species_colors)) {
  tips <- full_concat_tree$tip.label[tip_species == sp]
  node <- getMRCA(full_concat_tree, tips)  # this works on phylo object
  if(!is.null(node)) {
    p2 <- p2 + geom_hilight(node = node, fill = species_colors[sp], alpha = 0.3)
  }
}

p2 <- p2 + geom_tiplab()
p2


###### FIGURE 1B ######
# Load things in 
library(ape)
library(phangorn)
library(phytools)
library(ggplot2)

# Get in the folder results/RAxML/10Mb-concatenation/ 
# setwd("//wsl.localhost/Ubuntu/home/ntlin/phylo-practicum/results/RAxML/10Mb-concatenation/")

# List all .bestTree files. $ ensures the end of the name
tree_files <- list.files(pattern="\\.raxml.bestTree$")

# This will be a list with all the trees inside
trees <- list() 
# Make it a multiphylo object
class(trees)<- "multiPhylo"  

# Iterate through each file and read in the tree
for (i in 1:length(tree_files)) {
  trees[[i]] <- read.tree(tree_files[i])
}

# Root all trees by our outgroup
for(i in 1:length(trees)){
  trees[[i]] <- root(trees[[i]], outgroup="H_vulgare_HVens23", resolve.root=TRUE)
  # Make the tree ultrametric (all tips equidistant from the root) for nicer densitree 
  trees[[i]] <- chronos(trees[[i]])
}

# Make a consensus parsimony supertree
windows_consensus_st <- superTree(trees)
windows_consensus_st <- root(windows_consensus_st, "H_vulgare_HVens23", resolve.root=TRUE)

# Plot the density tree
ultrametric_full_concat_tree = chronos(full_concat_tree)

png(filename="../../../figures/figure1b.png", width=1800, height=900, units="px")
par(mfrow=c(1,2), mar=c(0.1, 0.1, 0.1, 0.1))
plot(ultrametric_full_concat_tree, show.tip.label=FALSE)
densiTree(trees, consensus=full_concat_tree, direction='leftwards', scaleX=TRUE, type='cladogram', alpha=0.1)
dev.off()


###### New ASTRAL Coalescent-Based Tree ######
# Set the working directory to the results/ folder
# setwd("//wsl.localhost/Ubuntu/home/ntlin/phylo-practicum/results")

indiv_species_astral = read.tree(file="07-individual-species-tree-astral4.tre")
indiv_species_astral = root(indiv_species_astral, outgroup="H_vulgare_HVens23", resolve.root=TRUE)

par(mfrow=c(1,2), mar = c(0.1, 0.1, 0.1, 0.1))
plot(full_concat_tree)
plot(indiv_species_astral)

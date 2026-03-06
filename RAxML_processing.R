### Load in everything 
library(ape)
library(phangorn)
library(phytools)
library(ggplot2)


### Check that working directory is the RAxML/ folder inside the results/ folder
getwd()
# setwd("//wsl.localhost/Ubuntu/home/ntlin/phylo-practicum/results/RAxML")

# Retrieve the names of tree files
tree_files <- list.files(pattern="\\.raxml.bestTree$")


### Read in the tree files
# Make a list for all the trees
gene_trees <- list()
# Make this a "multiPhylo" object for using it in other things
class(gene_trees) <- "multiPhylo"

# Read each file in tree_files into its corresponding position in gene_trees
for (i in 1:length(tree_files)) {
  gene_trees[[i]] <- read.tree(tree_files[[i]])
}


### Basic correctness checks 
# Glemin et al. filtered the trees to only include those that could be rooted at these species in order:
  # H_vulgare, Er_bonaepartis, S_vavilovii, and Ta_caputMedusae
# Not all gene trees have every taxon, so we're going to check that each tree has at least one of these taxa
# We'll also try to reroot it

# Our appropriate root taxa
root_taxa <- c("H_vulgare_HVens23", "Er_bonaepartis_TB1", "S_vavilovii_Tr279", "Ta_caputMedusae_TB2")

# Extract first matching species for each tree and put it in gene_tree_outgroup
# Initialize the vector to have one position for each gene tree
gene_tree_outgroup <- rep(NA,length(gene_trees))
# Iterate over each gene tree
for (i in 1:length(gene_trees)) {
  # Pull out the specific gene tree
  gene_tree <- gene_trees[[i]]
  # What are all the outgroup taxa we could find 
  found_taxa <- root_taxa[root_taxa %in% gene_tree$tip.label]
  
  # Return the first one if it exists
  if (length(found_taxa) > 0) {
    gene_tree_outgroup[i] <- found_taxa[1]
  }
}

# Check if any gene trees don't have an outgroup taxon (this should return FALSE)
any(is.na(gene_tree_outgroup))

# Re-root all gene trees by respective outgroup
for (i in 1:length(gene_trees)) {
  gene_trees[[i]] <- root(gene_trees[[i]],
                         outgroup = gene_tree_outgroup[i],
                         resolve.root=TRUE)
  # Make this ultrameric for nicer densitree? whatever that means
  gene_trees[[i]] <- chronos(gene_trees[[i]])
}

# Plot one gene tree
plot(gene_trees[[1]])



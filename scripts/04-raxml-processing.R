### Load in everything 
library(ape)
library(phangorn)
library(phytools)
library(ggplot2)


### Check that working directory is the RAxML/ folder inside the results/ folder
getwd()
# setwd("//wsl.localhost/Ubuntu/home/ntlin/phylo-practicum/results/RAxML/gene_trees")

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
  print(paste("Finished", tree_files[[i]]))
}

# Make this file in the results folder
write.tree(gene_trees, file="../../04-all-gene-trees.tre")


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


### Quantifying how often each taxon is present across our genes
# Count how often each individual appears
all_labels <- unlist(lapply(gene_trees, function(x) x$tip.label)) 
# Organize as data frame for plotting
df <- as.data.frame(table(all_labels) / length(gene_trees)) 

# Plot it
ggplot(df, aes(x=all_labels, y=Freq)) +
  geom_col() +
  labs(x="Individual", y="Proportion")+
  theme(axis.text.x = element_text(angle = 90))

### Making an initial supertree to evaluate phylogenetic signal
# Supertrees can handle missing data so we can use the entire dataset as is
st <- superTree(gene_trees)
st <- root(st, "H_vulgare_HVens23", resolve.root = T)
plot(st)

### Making densitrees for evaluating gene conflict/concordance
## Initial densitree shows a lot of conflict
densiTree(gene_trees, consensus=st, scaleX=T, type='cladogram', alpha=0.05)

## Making a reduced-taxa densitree for better interpretability
# Reducing to only taxa present in all genes
common_tips <- Reduce(intersect, lapply(gene_trees, function(tr) tr$tip.label))
# Only 6 common tips
length(common_tips)

# Prune the trees by dropping tips that are not in common_tips
trees_pruned <- lapply(gene_trees[1:10], function(tr) drop.tip(tr, setdiff(tr$tip.label, common_tips)))
# We can do this to remove NULLs
# trees_pruned <- trees_pruned[!sapply(trees_pruned, is.null)]

# Plot the densitree
densityTree(trees_pruned,use.edge.length=FALSE,type="cladogram",nodes="centered")

## Making a densitree for one alignment to see how our 20 RAxML runs differ
# Read in a single alignment file
trees = read.tree(file="Ae_bicornis_Tr406_BIS2_Contig10132_simExt_macseNT_noFS_clean.aln.raxml.mlTrees")

# Root the trees
rtrees <- lapply(trees, function(tr) {
  root(tr, outgroup = "H_vulgare_HVens23", resolve.root = TRUE)
})

# Plot densitrees before and after rooting
densityTree(trees,type="cladogram",nodes="intermediate")
densityTree(rtrees,type="cladogram",nodes="intermediate")

# Plot topology-only densitrees with no branch lengths
densityTree(trees,use.edge.length=FALSE,type="cladogram",nodes="centered")
densityTree(rtrees,use.edge.length=FALSE,type="cladogram",nodes="centered")
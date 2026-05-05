library(ape)
library(phangorn)

# Make sure to start in the results folder
# setwd("//wsl.localhost/Ubuntu/home/ntlin/phylo-practicum/results")

### SUPERTRIPLETS SUPERTREE ###
# Got this supertree file from supertriplets
our_supertree <- read.tree(file="07-supertree.tre")
# Can go straight to plotting
plot(our_supertree)

# Root the tree
rooted_supertree = root(our_supertree, outgroup="H_vulgare_HVens23", resolve.root=TRUE)
plot(rooted_supertree)

# We could compare it to a parsimony tree
  # But it takes too long so I am not gonna do it
# gene_trees <- read.tree(file="04-all-gene-trees.tre")
# supertree_parsimony <- superTree(gene_trees)
# rooted_supertree_parsimony <- root(supertree_parsimony, outgroup="H_vulgare_HVens23", resolve.root=TRUE)
# plot(rooted_supertree_parsimony)


### COALESCENT SPECIES TREES ###
## At the individual level
# Got this species tree file from ASTRAL4
astral4_indiv_species_tree <- read.tree(file="07-individual-species-tree-astral4.tre")
plot(astral4_indiv_species_tree)
rooted_astral4_indiv_species_tree <- root(astral4_indiv_species_tree, outgroup="H_vulgare_HVens23", resolve.root=TRUE)
plot(rooted_astral4_indiv_species_tree)

# At the species level
# Got this species tree file from ASTRAL4
astral4_sp_species_tree <- read.tree(file="07-species-level-species-tree-astral4.tre")
plot(astral4_sp_species_tree)
rooted_astral4_sp_species_tree <- root(astral4_sp_species_tree, outgroup="H_vulgare", resolve.root=TRUE)
plot(rooted_astral4_sp_species_tree)

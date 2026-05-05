library(ape)

# We can start in the "scripts" directory
setwd("//wsl.localhost/Ubuntu/home/ntlin/phylo-practicum/scripts")

# First get individual names - remember that OneCopyGenes is our full genes list
genes_directory <- "../data/Wheat_Relative_History_Data_Glemin_et_al/OneCopyGenes/"
# This lists the files that match the regex pattern "ends with .aln"
  # Double backslash because we need to escape it in R to use it as an escape character...
gene_files <- paste(genes_directory, list.files(genes_directory, pattern="\\.aln$"), sep="")

# Make a placeholder for individual names
all_individual_names <- character()
# Populate it
for (i in 1:length(gene_files)) {
  # Remember that the .aln file is an aligned fasta 
  names_from_fasta_headers <- rownames(read.dna(gene_files[i], format="fasta"))
  # Adds all the current file's names to the "all names" vector, then keeps the unique ones
  all_individual_names <- unique(c(all_individual_names, names_from_fasta_headers))
}

# Sort alphabetically for consistency
all_individual_names <- sort(all_individual_names)
print(all_individual_names)

# Remove the tag at the end of the individual name to get the species name
# This regex matches the last underscore and everything after it for deletion (replace w/ "")
  # "underscore then one or more of anything that is not an underscore till end of string"
no_tag_names <- sub("_[^_]+$", "", all_individual_names)

# Map individuals to species (just "individual_name species_name")
mapping <- paste(all_individual_names, no_tag_names)
# Write this mapping to a file (one line per entry)
writeLines(mapping, "../results/07-species-mapping.txt")

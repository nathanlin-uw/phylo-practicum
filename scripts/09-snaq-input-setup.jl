## This is a Julia script that should be run from the scripts/ folder

## Loading things in
using PhyloNetworks
using SNaQ
using CSV, DataFrames

## Reading in our previous species mapping file
mappingfile = CSV.read("../results/07-species-mapping.txt", DataFrame; header=false, delim=' ')

## We need to swap the species and individual columns
rename!(mappingfile, :Column1 => :individual)
rename!(mappingfile, :Column2 => :species)

select!(mappingfile,[:species, :individual])

## Remove 3 out of the 4 outgroups to simplify the analysis
## We keep `H_vulgare_HVens23`, remove `Ta_caputMedusae_TB2`, `Er_bonaepartis_TB1`, `S_vavilovii_Tr279`
filter(row -> row.species in ["Ta_caputMedusae", "Er_bonaepartis", "S_vavilovii"], mappingfile)
size(mappingfile) ## (47,2)

## Write the mapping to a new file
CSV.write("../results/09-species-mapping.csv", mappingfile)

## mappingfile = CSV.read("../results/09-species_mapping.csv", DataFrame)
taxonmap = Dict(r[:individual] => r[:species] for r in eachrow(mappingfile)) # as dictionary

## Now we read the gene trees and compute the CF table
trees = readmultinewick("../results/04-all-gene-trees.tre")
length(trees) ## 8708

for gt in trees
  for badtip in ["Ta_caputMedusae_TB2", "Er_bonaepartis_TB1", "S_vavilovii_Tr279"]
    if badtip in tiplabels(gt)
      deleteleaf!(gt, badtip)
    end
  end
end

# Write out new gene trees with the outgroups removed
writemultinewick(trees, "../results/09-all-gene-trees-snaq.tre")

## trees = readmultinewick("../results/09-all-gene-trees-snaq.tre")

## Creating and writing CF table:
df_sp = tablequartetCF(countquartetsintrees(trees, taxonmap; showprogressbar=false)...);
keys(df_sp)  # columns names
CSV.write("../results/09-tableCF-species.csv", df_sp); 

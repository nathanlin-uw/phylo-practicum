Recreating the analysis of Glemin et al. (2019) Pervasive hybridizations in the history of wheat relatives.

All scripts are stored in the `scripts/` folder.

### Getting the data ###
I navigated to https://www.agap-ge2pop.org/wheat-relative-history/ in the acknowledgements of the paper, downloaded the full archive, and moved it to and unzipped everything in the data/ folder.

### Setting up required software ###
Setting up a software folder
- This will contain all my executables so I will only need to add that folder to my PATH
- Opened up my \~/.bashrc file (\~ is shorthand for /home/ntlin/ folder) -- this configures shell setup things on startup
- Added this line to .bashrc: export PATH="$PATH:~/phylo-practicum/software"
- Now all executables in ~/phylo-practicum/software should be accessible

RAxML (Next Generation) v1.2.2
- Downloaded the Linux binary at https://github.com/amkozlov/raxml-ng/releases/download/1.2.2/raxml-ng_v1.2.2_linux_x86_64.zip
- Moved it to software/, unzipped it, and moved raxml-ng from the zip folder to software/
- Did `sudo chmod +x raxml-ng` to make it an executable
- Now I can use the command `raxml-ng`

SuperTriplets v1.1
- Downloaded SuperTriplets_v1.1.jar from https://www.agap-ge2pop.org/supertriplets/download/ 
- Placed it in a folder supertriplets/, which I moved to the software/ folder
- Now I can run this with `java -jar ~/phylo-practicum/software/supertriplets/SuperTriplets_v1.1.jar`

ASTRAL-IV v1.24.4.8
- Downloaded the Linux binary at https://github.com/chaoszhang/ASTER 
- Copied everything into the software/ folder
- Did `sudo chmod +x astral4` to make it an executable
- Now I can use the command `astral4`

Julia v1.20.1
- Did `curl -fsSL https://install.julialang.org | sh` in the terminal, then restarted the terminal
- Now I can type in `julia` to start the julia environment, and `exit()` to exit
- I used the Julia packages "PhyloNetworks", "SNaQ", "PhyloPlots", "CSV", and "DataFrames"
- In Julia, I can type `]` to get into package mode, then type `add [PackageName]`; backspace to go back to normal julia mode

HyDe v1.0.2
- Downloaded the files with `git clone https://github.com/pblischak/HyDe.git`
- Copied the HyDe folder into the software/ folder
- In the HyDe/ folder, I ran `python3 -m pip install -r requirements.txt` and `python3 -m pip install .`
- Now I can run HyDe with `run_hyde.py` anywhere in the terminal

R and Rstudio
I had R v4.5.1 and RStudio 2025.05.1+513 "Mariposa Orchid" Release for Windows already downloaded on my computer. In R, I downloaded the packages ape, phangorn, BiocManager, BiocManager:remotes, BiocManager:YuLab-SMU/treedataverse, MSCquartets, phytools, and ggplot2.

### Gene tree estimation with RAxML ###
I changed the script from the class repository, `scripts/04-raxml-gene-trees.sh` to make it so RAxML would not be repeated on the same files if it crashes and I need to run it again. For the "all_genes" run with the OneCopyGenes folder, it would save the names of the "finished" files in all_genes_log_completed_files.txt in the RAxML results directory. The core of this script is taking a specified data directory, iterating through its alignment files, and running raxml-ng on each, putting the output in `results/RAxML/gene_trees`. 

We are using the OneCopyGenes data, which consists of individual alignments for each locus.
For one gene, we'd do `raxml-ng --msa {filename.aln} --model GTR+G4`. We're using the General Time Reversible model so that sites can evolve at different rates with 4 gamma distributions for site variation. In the newest version of RAxML, we can do automatic model selection. 

Our main, full gene set is `OneCopyGenes`. For an abbreviated set of 10 genes, I made a subfolder of `data/Wheat_Relative_History_Data_Glemin_et_al/` called `OneCopyGenes-trimmed` and copied the first 10 genes into that folder using the following set of commands: `mkdir OneCopyGenes-trimmed; cd OneCopyGenes; ls | head -n10 | xargs -I {} cp "{}" ../OneCopyGenes-trimmed/`.

After running `04-raxml-gene-trees.sh` and generating gene trees in `results/RAxML/gene_trees/`, I visualized them using R code in the script `04-raxml-processing.R`. First, the script saves the list of all gene trees for future analyses in a file named `04-all-gene-trees.tre` in the `results/` folder. After that, this script reads the raxml.bestTree files in the gene trees directory into a list object, roots the trees based on either of the outgroups ("H_vulgare_HVens23", "Er_bonaepartis_TB1", "S_vavilovii_Tr279", or "Ta_caputMedusae_TB2"), then plots them in base R. The script also generates a summary of the proportion of gene trees different individuals appear in, creates and plots a supertree for overall signal and densitree for gene conflict, creates a new 10-gene densitree including only the taxa present in all genes for easier interpretability, and creates a densitree for one alignment to compare trees across the 20 RAxML runs. 

To summarize, this step uses the scripts `04-raxml-gene-trees.sh` and `04-raxml-processing.R` to produce, evaluate, and visualize gene trees, storing the RAxML output in `results/RAxML/gene_trees/` and creating a file containing the full set of gene trees, `results/04-all-gene-trees.tre`.

### Species tree with the supermatrix using the full set of genes ###
This step involved running RAxML on the supermatrix (concatenated gene set) provided by the authors, and comparing the result to our own concatenation.

The provided concatenation is located in `data/Wheat_Relative_History_Data_Glemin_et_al/`, and I ran `raxml-ng --msa triticeae_allindividuals_OneCopyGenes.fasta --model GTR+G4` in that folder, then, still in the same folder, moved the results to the `results/RAxML/full_concatenation/` folder with `mv *.raxml* ../../results/RAxML/full_concatenation/`.

To make my own concatenation, I used the R script `05-make-supermatrix.R`, which uses the ape package and creates a temporary Phylip before outputting a supermatrix called `manual-supermatrix.fasta` in the data directory. I then ran RAxML with the same GTR+G4 model.

To visualize the trees, I used similar code to that of `04-raxml-processing.R`, utilizing the `read.tree` and `root` function of the ape package and the base R plot function. To calculate RF distance, I used the phangorn package's function `RF.dist`.

### Species tree with the supermatrix using the 10Mb windows ###
In this step, I used concatenation of 10Mb windows rather than the set of full genes.

The provided concatenation alignments are located in `data/Wheat_Relative_History_Data_Glemin_et_al/Concatenation10Mb_OneCopyGenes/`, and I used the script `06-raxml-10mb-concatenation.sh` to run RAxML on each with a similar no-duplicate-runs approach to `04-raxml-gene-trees.sh`. 

I also created the 10Mb windows myself with the script `06-make-10Mb-windows.R`, which maps genes to chromosome positions and extracts and concatenates genes within each 10Mb sliding window of the Hordeum reference genome. I ran RAxML on these as well with `06-raxml-10mb-concatenation.sh`. 

For visualization of these trees, I used the script `06-windows-visualization.R`, which uses the packages ape, phangorn, phytools, and ggplot2 to read in, root, and make ultrametric each tree, create and root a supertree, and plot a densitree. The script also reads in, roots, makes ultrametric, and plots the supertree and densitree for the trees provided by the authors in `Densitree_OneCopyGenes.nex`. I modified the nexus file before using these trees, adding a semicolon after "begin trees", and changing taxon 0 to 47 in the translate block and trees.

### Species tree supertree/coalescent ###
I first ran SuperTriplets to produce a supertree using a triplets-based approach with the command `java -jar ~/software/SuperTriplets_v1.1.jar 04-all_gene_trees.tre 07-supertree.tre`, and compared it to the parsimony supertree with the R script `07-supertree-visualization.R`.

To make an individual level species tree using the coalescent model, I used ASTRAL4 with the command `astral4 -i 04-all_gene_trees.tre -o 07-individual-species-tree-astral4.tre`. This produced the tree file `07-individual-species-tree-astral4.tre` in the `results/` folder. To make a species level species tree using the coalescent model, I did some pre-processing (mapping individual names to species) using the R script `07-make-species-mapping.R`, then ran ASTRAL4 with `astral4 -i 04-all_gene_trees.tre -a 07-species_mapping.txt -o 07-species-tree-astral4.tre`. I visualized both of these species trees as well in `07-supertree-visualization.R`.

### Additional species tree visualizations ###
I attempted to reproduce figures 1A and 1B of Glemin et al. (the full concatenation tree/supertree and sliding window concatenation densitree), as well as a coalescent-based species tree, all in the script `08-making-figures.R`. 

For Figure 1A, the script uses the ape, ggtree, and dplyr packages to do some rotations of the ggtree plot and color clades by species. Notably, the authors had claimed that the supertree and concatenation tree were equivalent, but reordering and rotating the trees, as well as calculating RF distance, demonstrated that they were actually not equivalent.

For Figure 1B, the script uses the ape, phangorn, phytools, and ggplot2 packages to read in the trees, root them, and make them ultrametric. It then creates a supertree, roots it, and makes a densitree as out final product.

I plotted the ASTRAL tree in base R and compared it to the concatenation tree, finding that they are also not equivalent according to RF distance.

### Species network ###
Phylogenetic networks provide information on gene flow, allowing for better understanding of evolutionary histories beyond concatenation and coalescent methods alone. I used SNaQ to infer the phylogenetic network based on the starting ASTRAL species level species tree and a table of quartet concordance factors, and visualized it in Julia using the script `09-run-snaq.jl`. To set up the concordance factors table, I used the script `09-snaq-input-setup.jl`, which output that table of quartet concordance factors, `09-tableCF-species.csv`, from the individuals to species mapping file and the complete set of gene trees. 

The main tree here matched the tree from concatenation and the coalescent model, so that was a good sign.

### Hybrid detection with HyDe ###
I used HyDe to look for signals of hybridization among parent-parent-offspring triplet combinations. 

I first converted the fasta data to the Phylip format that HyDe uses using the `10a-to-phylip.py` script. This script uses BioPython's AlignIO's SequentialPhylipWriter, and wrote the Phylip output to `results/10-triticeae_allindividuals_OneCopyGenes.phylip`. I had to first `sudo apt install python3-pip` and install BioPython with `pip install BioPython` in the terminal -- that gave me an error that the environment is externally managed. What worked after that was just using conda with `conda activate base` (apparently BioPython was already installed in my base environment).

With the input file set up, I ran HyDe with `run_hyde.py -i ../results/10-triticeae_allindividuals_OneCopyGenes.phylip -m ../results/07-species_mapping.txt -o H_vulgare -n 47 -t 17 -s 11354214 --prefix 10-hyde`. This command tells HyDe that we are using H_vulgare as the outgroup, with 47 individuals sampled, 17 taxa, and 11354214 sites. It generated the result files `10-hyde-out-filtered.txt` and `10-hyde-out.txt`, which I moved to the `results/` folder. 

We were going to repeat this for the 10Mb windows, but the code for the Phylip conversion did not work, so we did not end up continuing with that. If we had, we would have used a script, `10-hyde-10mb.sh`, to iterate through each window and run HyDe on the respective gene set.

### Hybrid detection with MSCQuartets and Visualization ###
Our final week was an alternative method of detecting hybridization signal -- this one looks at whether quartets display tree-like evolution, and evaluates whether major/minor topology frequencies cannot be explained by Incomplete Lineage Sorting alone.

I used the script `11-run-mscquartets.R`, which uses the gene trees as input and creates a table of p values that can be visualized with a violin plot. This script used the MSCquartets package for the main analysis, then dplyr for cleaning up the table and ggplot2 for creating the violin plot. 

Additionally, I had visualized the results from HyDe to recreate figures 3A and 3B using the script `11-visualize-hyde.R`. This script uses the dplyr, purrr, and ggplot2 packages to clean up the HyDe output, go through the positions on chromosome 3, filter for hybrids, identify hybridization signals for each topology, filter for focal species, and then make the violin boxplots for hybridization index for various taxon combinations and scatterplot for hybridization index along the third chromosome.

# Thanks for a great semester!!

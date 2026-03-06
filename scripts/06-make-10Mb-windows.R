library(ape)

library(glue)

getwd() #should be in the code folder
setwd("//wsl.localhost/Ubuntu/home/ntlin/phylo-practicum/scripts") # if not in the correct folder

genes_dir <- "../data/Wheat_Relative_History_Data_Glemin_et_al/OneCopyGenes/"
concat_10_dir <-"../data/Wheat_Relative_History_Data_Glemin_et_al/OneCopyGenes-10mb/"

dir.create(concat_10_dir)

gene_map_loc <- "../data/Wheat_Relative_History_Data_Glemin_et_al/MappedOnHordeumChromosomes_OneCopyGenes.txt"
gene_map <- read.table(gene_map_loc,header = T)
mb10 <- 10000000 #10mb

window_num<-1
for(chromosome in unique(gene_map$Chrom)){

  chrom_genes<- gene_map[gene_map$Chrom==chromosome,] #get all genes on the chromosome
  chrom_len <- max(chrom_genes$Pos)
  
  start_pos<-1
  end_pos <-mb10
  
  while(end_pos < chrom_len){
    #find genes in the window
    window_genes <- chrom_genes[chrom_genes$Pos>=start_pos & chrom_genes$Pos< end_pos,]$Name

    # I want to see how many genes are in each chromosome/position combo (each window)
    print(glue("We're at chromosome {chromosome} position {start_pos} and we have {length(window_genes)} genes in this window"))
    
    
    if(length(window_genes>2)){ #only look at windows with 3+ genes
      
      window_genes <- paste(genes_dir,window_genes,sep='')
      alignments <- lapply(window_genes, read.dna, format = "fasta")
      
      
      
      ##THIS GIVES AN w ERROR - same name in alignment file 
      #concat_window <- do.call(cbind, c(alignments, fill.with.gaps = TRUE))
      
      #first remove genes with duplicate records
      # Clean the list by keeping only unique names for each gene
      alignments <- lapply(alignments, function(x) {
        x[!duplicated(rownames(x)), ]
      })
      concat_window <- do.call(cbind, c(alignments, fill.with.gaps = TRUE))
      
      
      write.dna(concat_window, file = paste(concat_10_dir,'window_',window_num,'.aln'), format = "fasta")
      window_num<- window_num+1
    }
    start_pos<-start_pos+mb10
    end_pos <-end_pos+mb10
  }
  
}


##################################################
### Auxilliary script to check duplicate names ###
##################################################

genes <- list.files(genes_dir,pattern = '\\.aln$')


dup_inds <-rep(NA,length(genes))

i<-1
for(gene in genes){
  aln <- read.dna(file=paste(genes_dir,gene,sep=''),format='fasta')
  dup_inds[i]<-any(duplicated(rownames(aln)))
  i<-i+1
  if(i %% 500 ==0){
    print(i)
  }
}



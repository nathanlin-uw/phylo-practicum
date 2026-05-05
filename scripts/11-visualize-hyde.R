library(dplyr)
library(purrr)
library(ggplot2)

parents1 = c("Ae_mutica", "Ae_speltoides")## B Clade
parents2 = c("T_boeoticum", "T_urartu") #A Clade
hybrids = c("Ae_bicornis", "Ae_longissima", "Ae_sharonensis", "Ae_searsii", "Ae_tauschii", "Ae_caudata", "Ae_umbellulata", "Ae_comosa", "Ae_uniaristata") # D clade



hyde_chrom_dir <-"../results/HyDe/10Mb-concatenation-ch3/"
chrom_files<-list.files(hyde_chrom_dir,pattern='out.txt')

all_hyde <-list()
for(chrom_file in chrom_files){ # go thru all positions on Chromosome 3
  hyde_file<-paste(hyde_chrom_dir,chrom_file,sep='')
  hyde_result <-read.table(hyde_file,header = T)

  chrom_pos <- sub(".*_(\\d+)-out\\.txt", "\\1", chrom_file)
  
  filtered_results <- hyde_result %>%
    filter(Hybrid %in% hybrids) %>% #get the rows with hybrid species as the hybrid
    filter( #specify the parents
      (P1 %in% parents1 & P2 %in% parents2) | 
        (P1 %in% parents2 & P2 %in% parents1)
    )  %>%
    mutate(P1_signal = ABBA,
           P2_signal = AABB,
           P12_signal= ABAB)
    if(nrow(filtered_results)==0){
      next
    }
  
  all_hyde[[chrom_pos]]<-filtered_results
}
all_hyde <- bind_rows(all_hyde,.id='pos') #combine results into a single data frame

mutica <- all_hyde %>% 
  filter((P1 %in% "Ae_mutica") | (P2 %in% "Ae_mutica"))  %>%
  mutate(B_signal = ifelse(P1 %in% parents1, P1_signal , P2_signal), #Signal for the topology (B,D),A
         A_signal = ifelse(P1 %in% parents2, P1_signal , P2_signal), #Signal for the topology (A,D),B
         AB_signal= P12_signal, #Signal for the topology (A,B),D
         h_index = (A_signal-AB_signal)/(A_signal+B_signal-(2*AB_signal)) # produces values outside [0,1]
         ) %>%
  mutate(Gamma = ifelse((Gamma<0.0) | (Gamma>1.0),NA,Gamma)) %>%
  mutate(h_index = ifelse((h_index<0.0) | (h_index>1.0),NA,h_index)) 



ggplot(mutica, aes(x = as.numeric(pos), y = Gamma, color)) +
  stat_summary(fun.data = "mean_cl_normal", #confidence interval
               geom = "errorbar", 
               width = 0.2, 
               color = "royalblue") +
  
  stat_summary(fun = "mean",  #mean dots
               geom = "point", 
               size = 2, 
               color = "black") +
  
  # Add the horizontal dashed line at 0.5
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  
  # Set Y-axis limits to represent probability
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 10)) +
  theme_bw() +
  labs(
    x = "Genomic Position (Window Index)", 
    y = expression(gamma),
    title = "Inheritance Probability across Chromosome 3"
  ) 

#filter entries with the focal species.
A_B_species<- c("Ae_mutica","T_boeoticum")
mutica_boeoticum <- mutica %>% 
  filter((P1 %in% A_B_species) & (P2 %in% A_B_species) )
  
##Make the order of species match Fig3A
target_order <- c(
  "Ae_bicornis", "Ae_longissima", "Ae_sharonensis", "Ae_searsii", 
  "Ae_tauschii", "Ae_caudata", "Ae_umbellulata", "Ae_comosa", "Ae_uniaristata"
)
mutica_boeoticum <- mutica_boeoticum %>%
  mutate(Hybrid = factor(Hybrid, levels = target_order))



ggplot(mutica_boeoticum, aes(x = Hybrid, y = Gamma, fill = Hybrid)) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    color = NA
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  scale_fill_manual(values = hybrid_colors, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dotted",
    color = "black"
  ) +
  labs(
    x = "D focal species",
    y = expression(gamma)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

# If we wanted to filter by p value
significant_mutica_boeoticum <- mutica_boeoticum %>%
	filter(Pvalue<1e-6)

# If we wanted to specify focal species
rm(list = ls())
setwd("~/sepsis/workspace/karius/datasets/")
require(tidyr)
df <- read.csv("karius_genus_pathogens.csv", sep = ",", stringsAsFactors = F)
df <- df[df$pathogen %in% c("Escherichia coli", "Enterococcus faecalis", "Proteus mirabilis"), ]
X <- df[, c('Lymphocryptovirus', 'Campylobacter', 'Shigella', 
            'Bacillus', 'Veillonella', 'Betapolyomavirus', 
            'Tannerella', 'Salmonella', 'Aeromonas', 
            'Bifidobacterium', 'Alphatorquevirus', 'Shewanella', 
            'Cytomegalovirus', 'Escherichia', 'Prevotella', 
            'Lachnoclostridium', 'Proteus', 'Enterobacter', 
            'Burkholderia', 'Pandoraea', 'Blautia', 
            'Enterococcus', 'Stenotrophomonas')]

rowsums <- apply(X, 1, sum)
X <- X / rowsums

# Mean RA of septic and healthy
means <- aggregate(X, by=list(df$y), mean)
sdevs <- aggregate(X, by=list(df$y), sd)

require(ggplot2)
require(reshape2)
plot_df <- data.frame(Escherichia = X$Escherichia, 
                      Enterococcus = X$Enterococcus, 
                      Proteus = X$Proteus, 
                      y = df$pathogen, 
                      check.names = F, stringsAsFactors = F)
plot_df <- melt(plot_df)

ggplot(plot_df, aes(x = y, y = value, fill = variable)) +
  geom_boxplot() +
  labs(x = "Causative Pathogen", y = "Relative Abundance", fill = "Feature") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave("../results/boxplot_genera.png", dpi = 600, width = 7)

rm(list = ls())
setwd("~/sepsis/workspace/karius/datasets/")
df <- read.csv("karius_genus_pathogens.csv", sep = ",", stringsAsFactors = F)
df$pathogen[df$pathogen == "Human herpesvirus 4"] <- "Lymphocryptovirus"
df$pathogen[df$pathogen == "Human herpesvirus 1"] <- "Simplexvirus"
df$pathogen[df$pathogen == "Human herpesvirus 5"] <- "Cytomegalovirus"

# Partition by Genus
require(tidyr)
split_cols <- separate(df, col = pathogen, sep = " ", into = c("Genus", "Species"))
y <- split_cols$Genus
y <- y[y != "none"]

require(reshape2)
plot_df <- melt(table(y))

require(ggplot2)
require(pals)
ggplot(plot_df, aes(x = y, y = value, fill = y)) +
  geom_bar(stat = "identity") +
  labs(y = "No. of Samples") +
  scale_fill_manual(values = as.vector(polychrome(length(unique(y))))) +
  geom_text(aes(label = value), position = position_dodge(width=0.9), vjust=-0.25, size = 3) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank())

ggsave("../results/confirmed_pathogens.png", dpi = 600, width = 6, height = 4)
  


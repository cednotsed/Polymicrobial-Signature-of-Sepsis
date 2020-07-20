rm(list = ls())
setwd("~/git_repos/Polymicrobial-Signature-of-Sepsis/")
df <- read.csv("results/SHAP_values_per_infection.csv", sep = ",", stringsAsFactors = F)
head(df)
df$color <- "Positive"
df$color[df$shap < 0] <- "Negative"

require(ggplot2)

plt1  <- ggplot(df, aes(x = Genus, y = shap, color = color)) +
         geom_point() +
         geom_hline(yintercept = 0, linetype = "dotted") +
         labs(y = "SHAP Values", x = "Genus") +
         theme_bw() + 
         theme(axis.ticks.x = element_blank(),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
               legend.position = "none")

## Boxplot of genera abundance
df2 <- read.csv("datasets/karius_genus_raw_maxi.csv", sep = ",", stringsAsFactors = F)
genera <- unique(df$Genus)
df2 <- df2[]
# ggsave("results/shap_by_infection.png", dpi = 600, width = 6, height = 4)
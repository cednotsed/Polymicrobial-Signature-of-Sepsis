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
         theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               legend.position = "none")
         

## Boxplot of genera abundance
df$abundance <- log(df$abundance)

plt2  <- ggplot(df, aes(x = Genus, y = abundance)) +
         geom_boxplot() +
         labs(y = "ln(read count)", x = "Genus") +
         geom_hline(yintercept = 0, linetype = "dotted") +
         theme_bw() +
         theme(axis.ticks.x = element_blank(),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
               legend.position = "none")

require(ggpubr)         
ggarrange(plt1, plt2, 
          ncol = 1, nrow = 2, 
          common.legend = T, 
          align = "v",
          labels = "auto",
          label.y = 1.1)
ggsave("results/shap_by_infection.png", dpi = 600, width = 6, height = 5)

# Plot most abundant genera in NTCs
rm(list=ls())
setwd('~/git_repos/Polymicrobial-Signature-of-Sepsis')

df <- read.csv('datasets/kapusta_genus_raw_silva.csv')
df <- df[df$y == "ntc", ]
df <- df[, c(colnames(df)[!sapply(df, function(x){all(x == 0)})])]
means <- apply(df[, 2:ncol(df)], 2, mean)
means <- means[order(means, decreasing = T)][1:20]

require(ggplot2)
require(reshape2)
melted <- melt(means)
melted$variable <- rownames(melted)
melted$variable <- factor(melted$variable, levels = melted$variable[order(melted$value, decreasing = T)])
ggplot(melted, aes(x = variable, y = value, fill = value)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,2), "cm")) + 
  labs(x = "Genus", y = "Average Read Count (n = 3)") +
  scale_fill_gradient(high = "red", low = "blue")

ggsave("results/external_validation/kapusta_ntc_most_abundant.png", dpi=600, width = 6, height = 5)

rm(list = ls())
setwd("~/git_repos/Polymicrobial-Signature-of-Sepsis/datasets")
maxi <- read.csv("kapusta_order_raw_maxi.csv", stringsAsFactors = F)
silva <- read.csv("kapusta_order_raw_silva.csv", stringsAsFactors = F)
common_cols <- intersect(colnames(maxi), colnames(silva))
common_cols <- common_cols[common_cols != "y"]

y <- as.character(maxi[, "y"])
y[y == "healthy"] <- "Healthy"
y[y == "septic"] <- "Septic"
y[y == "ntc"] <- "NTC"

maxi <- maxi[, 2:ncol(maxi)]
maxi <- maxi[, common_cols]
maxi <- maxi / apply(maxi, 1, sum) * 100
maxi <- aggregate(maxi, by=list(y), mean)

silva <- silva[, 2:ncol(silva)]
silva <- silva[, common_cols]
silva <- silva / apply(silva, 1, sum) * 100
silva <- aggregate(silva, by=list(y), mean)

require(ggplot2)
require(reshape2)
require(RColorBrewer)
pal <- brewer.pal(8, "Set2")

# Maxikraken db
melted_maxi <- melt(maxi)

# Get most abundant orders
means <- aggregate(melted_maxi$value, by=list(melted_maxi$variable), mean)
orders <- as.character(means$Group.1)[order(means$x, decreasing = T)][1:8]

melted_maxi <- melt(maxi[, colnames(maxi) %in% c("Group.1", orders)])
melted_silva <- melt(silva[, colnames(silva) %in% c("Group.1", orders)])

# Set factor sequence
melted_maxi$variable <- factor(melted_maxi$variable, levels = orders)
melted_silva$variable <- factor(melted_silva$variable, levels = orders)

plt1 <- ggplot(melted_maxi, aes(y = value, x=Group.1, fill=variable)) +
  geom_bar(stat="identity") +
  labs(fill = "Order", x = "MaxiKraken2", y = "Mean Relative Abundance (%)") +
  scale_fill_manual(values = pal) +
  theme(plot.margin = unit(c(0,0,0,1), "cm"))

plt2 <- ggplot(melted_silva, aes(y = value, x=Group.1, fill=variable)) +
  geom_bar(stat="identity") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,0,0,1), "cm")) +
  labs(fill = "Order", x = "Silva") +
  scale_fill_manual(values = pal)

# Plot no. of minimisers in database
silva <- read.csv("silva_290620_info.txt", sep = "\t", stringsAsFactors = F, header = F)
maxi <- read.csv("maxikraken2_info.txt", sep = "\t", stringsAsFactors = F, header = F)

# Trim whitespace
silva$V6 <- trimws(silva$V6)
maxi$V6 <- trimws(maxi$V6)

# Get orders
silva_df <- silva[silva$V6 %in% orders, c("V1", "V6")]
silva_df$Database <- "Silva"
maxi_df <- maxi[maxi$V6 %in% orders, c("V1", "V6")]
maxi_df$Database <- "MaxiKraken2"

plot_df <- rbind(silva_df, maxi_df)
plot_df$V6 <- factor(plot_df$V6, levels = orders)

plt3 <- ggplot(plot_df, aes(y = V1, x = Database, fill = V6)) +
        geom_bar(stat="identity") +
        labs(fill = "Order", x = "Database", y = "Proportion of minimisers (%)") +
        scale_fill_manual(values = pal) +
        theme(plot.margin = unit(c(0,0,0,1), "cm"))

require(ggpubr)
ggarrange(plt1, plt2, plt3, 
          ncol = 2, nrow = 2, 
          common.legend = T, 
          align = "hv",
          labels = "auto")
ggsave("../results/silva_versus_maxi_barplot.png", dpi = 600, width = 8, height = 5)

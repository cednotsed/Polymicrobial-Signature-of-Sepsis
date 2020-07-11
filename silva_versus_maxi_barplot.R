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
melted_maxi <- melted_maxi[melted_maxi$variable %in% orders,]

# Get most abundant orders
means <- aggregate(melted_maxi$value, by=list(melted_maxi$variable), mean)
orders <- as.character(means$Group.1)[order(means$x, decreasing = T)][1:8]

# Silva db
melted_silva <- melt(silva)
melted_silva <- melted_silva[melted_silva$variable %in% orders,]

plt1 <- ggplot(melted_maxi, aes(y = value, x=Group.1, fill=variable)) +
  geom_bar(stat="identity") +
  labs(fill = "Order", x = "MaxiKraken", y = "Mean Relative Abundance (%)") +
  scale_fill_manual(values = pal)

plt2 <- ggplot(melted_silva, aes(y = value, x=Group.1, fill=variable)) +
  geom_bar(stat="identity") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(fill = "Order", x = "Silva") +
  scale_fill_manual(values = pal)

require(ggpubr)
ggarrange(plt1, plt2, ncol = 2, nrow = 1, common.legend = T, align = "hv")
ggsave("../results/silva_versus_maxi_barplot.png", dpi = 600, width = 8, height = 5)

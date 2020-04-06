rm(list = ls())
setwd("~/sepsis/workspace/karius/datasets/")
require(tidyr)
df <- read.csv("karius_genus_pathogens.csv", sep = ",", stringsAsFactors = F)
X <- df[, c('Lymphocryptovirus', 'Campylobacter', 'Shigella', 
            'Bacillus', 'Veillonella', 'Betapolyomavirus', 
            'Tannerella', 'Salmonella', 'Aeromonas', 
            'Bifidobacterium', 'Alphatorquevirus', 'Shewanella', 
            'Cytomegalovirus', 'Escherichia', 'Prevotella', 
            'Lachnoclostridium', 'Proteus', 'Enterobacter', 
            'Burkholderia', 'Pandoraea', 'Blautia', 
            'Enterococcus', 'Stenotrophomonas')]


# Partition by Genus
split_cols <- separate(df, col = pathogen, sep = " ", into = c("Genus", "Species"))
y <- split_cols$Genus
y[y == "Human"] <- "Herpesviridae"
y[y == "Cytomegalovirus"] <- "Herpesviridae"
y[y == "none"] <- "Healthy Control"

# Count number of samples with non-zero k-mer counts
get_number <- function(df) {
  test <- data.frame(df)
  test[test > 0] <- 1
  test[test == 0] <- 0
  test <- apply(test, 2, sum)
  return(test)
}

numbers <- get_number(X)

# Get RA
rowsums <- apply(X, 1, sum)
X <- X / rowsums

# Remove Signals
X_new <- X[, colnames(X) != "Escherichia"]
X_new <- X_new[, colnames(X_new) != "Proteus"]
X_new <- X_new[, colnames(X_new) != "Enterococcus"]
X_new <- X_new[, colnames(X_new) != "Cytomegalovirus"]
X_new <- X_new[, colnames(X_new) != "Lymphocryptovirus"]


# t-SNE
set.seed(66)
require(Rtsne)

perp <- 15

tsne <- Rtsne(as.matrix(X),
              verbose = T,
              perplexity = perp,
              max_iter = 10000,
              pca = F,
              theta = 0)

tsne_new <- Rtsne(as.matrix(X_new),
              verbose = T,
              perplexity = perp,
              max_iter = 10000,
              pca = F,
              theta = 0)

plot_ori <- data.frame(tsne$Y, Pathogen=y)
plot_new <- data.frame(tsne_new$Y, Pathogen=y)

# Params
size <- 2
pal <- rep("grey70", 22)
pal[8] <- "red"
pal[11] <- "orange"
pal[16] <- "purple"
pal[7] <- "blue"
pal[12] <- "black"

require(ggpubr)
plt1 <- ggplot(plot_ori, aes(x = X1, y = X2, color = Pathogen)) +
        geom_point(size = size, alpha = 0.8) +
        scale_color_manual(values = pal) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.ticks = element_blank(),
              legend.title = element_blank())

plt2 <- ggplot(plot_new, aes(x = X1, y = X2, color = Pathogen)) +
  geom_point(size = size, alpha = 0.8) +
  scale_color_manual(values = pal) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank())

require(ggpubr)

ggarrange(plt1, plt2, 
          ncol = 2, nrow = 1, 
          common.legend = T,
          legend = "bottom",
          labels = "auto",
          align = "hv",
          hjust = 0)

ggsave("../results/pca_tsne_septic.png", dpi = 600, width = 8, height = 5)

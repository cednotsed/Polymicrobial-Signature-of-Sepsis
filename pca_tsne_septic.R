rm(list = ls())
setwd("~/git_repos/Polymicrobial-Signature-of-Sepsis/")
require(tidyr)
df <- read.csv("datasets/karius_genus_raw.csv", sep = ",", stringsAsFactors = F)
X <- df[, c('Campylobacter', 'Shigella', 'Escherichia',
            'Prevotella', 'Enterobacter', 'Burkholderia',
            'Streptococcus', 'Moraxella', 'Bacteroides',
            'Stenotrophomonas')]
y <- df$y

# Partition by Genus
split_cols <- separate(df, col = pathogen, sep = " ", into = c("Genus", "Species"))
pathogens <- split_cols$Genus
pathogens[pathogens == "Human"] <- "Herpesviridae"
pathogens[pathogens == "Cytomegalovirus"] <- "Herpesviridae"
pathogens[pathogens == "none"] <- "Healthy Control"

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
to_remove <- c('Escherichia', 'Streptococcus', 'Mycobacterium', 'Cytomegalovirus',
              'Staphylococcus', 'Proteus', 'Klebsiella', 'Pseudomonas',
              'Moraxella', 'Enterococcus', 'Enterobacter', 'Citrobacter',
              'Haemophilus', 'Fusobacterium', 'Salmonella', 'Serratia',
              'Aerococcus', 'Campylobacter', 'Lymphocryptovirus', 'Simplexvirus')

X_new <- X[, !(colnames(X) %in% to_remove)]

# t-SNE
set.seed(66)
require(Rtsne)
require(vegan)

# Get Bray-Curtis pairwise matrix
bc <- vegdist(X, method = "bray")
bc_new <- vegdist(X_new, method = "bray")

perp <- 40

tsne <- Rtsne(bc,
              verbose = T,
              perplexity = perp,
              max_iter = 10000,
              is_distance = T,
              pca = T,
              theta = 0)

tsne_new <- Rtsne(bc_new,
              verbose = T,
              perplexity = perp,
              max_iter = 10000,
              is_distance = T, 
              pca = T,
              theta = 0)

plot_ori <- data.frame(tsne$Y, Pathogen = pathogens)
plot_new <- data.frame(tsne_new$Y, Pathogen = pathogens)

# Params
size <- 2

pal <- rep("grey70", 22)
pal[which(levels(plot_ori$Pathogen) == "Escherichia")] <- "red"
# pal[11] <- "orange"
# pal[16] <- "purple"
# pal[7] <- "blue"
# pal[12] <- "black"

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

ggsave("results/decontam/pca_tsne_septic.png", dpi = 600, width = 8, height = 5)

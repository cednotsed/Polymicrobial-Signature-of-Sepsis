rm(list=ls())
#library(devtools)
#install_github("hallucigenia-sparsa/seqgroup")
require('seqgroup')
setwd('~/git_repos/Polymicrobial-Signature-of-Sepsis')
df <- read.csv('datasets/karius_genus_raw.csv')
meta <- read.csv("datasets/pathogen_list.csv")
meta <- meta$Genus

pathogens <- c('Escherichia', 'Prevotella', 'Bacteroides', 'Lactobacillus', 
               'Stenotrophomonas', 'Lymphocryptovirus', 'Shigella', 
               'Burkholderia', 'Cellulomonas', 'Paracoccus')

X <- df[, 3:ncol(df)]
X <- X[, colnames(X) %in% meta]
septic <- X[df$y == 'septic' & df$pathogen == "Escherichia coli", ]
healthy <- X[df$y == 'healthy', ]
healthy <- sample(healthy, size = length(septic))

rho_vec <- c()
p_vec <- c()
name_vec <- c()
test_df <- rbind(septic, healthy)

for (pathogen in colnames(X)) {
  test <- cor.test(test_df$Escherichia, test_df[, pathogen], method = "spearman", exact = F)
  rho <- test$estimate
  p <- test$p.value
  if (!is.na(rho)) {
    rho_vec <- c(rho_vec, rho)
    p_vec <- c(p_vec, p)
    name_vec <- c(name_vec, pathogen)
  }
}

results <- data.frame(rho = rho_vec, p = p_vec)
rownames(results) <- name_vec
results$p <- results$p * nrow(results)
results <- results[results$p < 0.05, ]
results <- results[order(results$rho, decreasing = T), ]

edges <- c()
weights <- c()
color <- c()

for (pathogen in rownames(results)[2:nrow(results)]) {
  edges <- c(edges, "Escherichia", pathogen)
  rho <- results[pathogen, "rho"]
  weights <- c(weights, rho)
  color <- c(color, ifelse(rho < 0, "green", "black"))
}

g <- graph(edges=edges, directed = F)

E(g)$weight <- weights
E(g)$color <- ifelse(E(g)$weight < 0, "red", "green")
l <- layout_with_fr(g, weights=weights)

plot(g,
     vertex.frame.color = "white",
     vertex.size = 8,
     vertex.label.color = "black",
     vertex.label.cex = 1.1,
     vertex.label.font = 1.8,
     vertex.label.family = "Arial",
     vertex.label.dist = 0.2,
     edge.width = abs(E(g)$weight) * 10,
     layout=layout.star)

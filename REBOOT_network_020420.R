rm(list=ls())
#library(devtools)
#install_github("hallucigenia-sparsa/seqgroup")
require('seqgroup')
setwd('sepsis/workspace/karius/datasets')
df <- read.csv('karius_genus_pathogens.csv')
X <- df[, c('Lymphocryptovirus', 'Campylobacter', 'Shigella', 
            'Bacillus', 'Veillonella', 'Betapolyomavirus', 
            'Tannerella', 'Salmonella', 'Aeromonas', 
            'Bifidobacterium', 'Alphatorquevirus', 'Shewanella', 
            'Cytomegalovirus', 'Escherichia', 'Prevotella', 
            'Lachnoclostridium', 'Proteus', 'Enterobacter', 
            'Burkholderia', 'Pandoraea', 'Blautia', 
            'Enterococcus', 'Stenotrophomonas')]

septic <- X[df$y == 'septic', ]
healthy <- X[df$y == 'healthy', ]
healthy <- sample(healthy, size = length(septic))

# Tranpose to get taxa as rows samples as columns
septic <- t(septic) 
healthy <- t(healthy)

# Count number of samples with non-zero k-mer counts
get_number <- function(df) {
        test <- data.frame(df)
        test[test > 0] <- 1
        test[test == 0] <- 0
        test <- apply(test, 2, sum)
        return(test)
}

numbers <- get_number(X)

# Plot Septic network
set.seed(69)
edges <- 40
pval <- 0.05
width <- 1

# Plot septic graph
septic_graph <- barebonesCoNet(septic,
                               norm = T,
                               methods = c("bray", "spearman", "kld"),
                               init.edge.num = edges,
                               pval.cor = F,
                               pval.T = pval,
                               renorm = T,
                               permutandboot = T,
                               iters = 1000,
                               bh = T, 
                               plot = F,
                               keep.filtered = F,
                               # min.occ = ncol(septic) / 10,
                               verbose = T)

par(mar = c(0,0,0,0))
png(file = "../results/sepsis_colored_network.png",width = 10, height = 8, units = 'in', res=300)

layout <- layout_(septic_graph, with_dh(weight.edge.lengths = edge_density(septic_graph)/1000))
plot(septic_graph,
     vertex.frame.color = "white",
     vertex.size = 8,
     vertex.label.color = "black",
     vertex.label.cex = 0.9,
     vertex.label.font = 1.8,
     vertex.label.family ="Arial",
     vertex.label.dist = 0.2,
     edge.width = E(septic_graph)$weight * width,
     edge.curved = 0.2,
     layout = layout,
     main = 'Microbial Network for Septic Patients')

dev.off()

# Plot healthy network
healthy_graph <- barebonesCoNet(healthy,
                                norm = T,
                                methods = c("bray", "spearman", "kld"),
                                init.edge.num = edges,
                                pval.cor = F,
                                pval.T = pval,
                                renorm = T,
                                permutandboot = T,
                                iters = 1000,
                                bh = T, 
                                plot = F,
                                keep.filtered = F,
                                # min.occ = ncol(healthy) / 10,
                                verbose = T)

layout <- layout_(healthy_graph, with_dh(weight.edge.lengths = edge_density(healthy_graph)/100))
par(mar = c(0,0,0,0))
png(file = "../results/healthy_colored_network.png", width = 10, height = 8, units = 'in', res = 300)

plot(healthy_graph,
     vertex.frame.color = "white",
     vertex.size = 8,
     vertex.label.color = "black",
     vertex.label.cex = 0.9,
     vertex.label.font = 1.8,
     vertex.label.family = "Arial",
     vertex.label.dist = 0.2,
     edge.width = E(healthy_graph)$weight * width,
     edge.curved = 0.2,
     layout = layout,
     main='Microbial Network for Healthy Patients')

dev.off()

# Substract healthy correlations
edge_septic <- as_ids(E(septic_graph))
edge_healthy <- as_ids(E(healthy_graph))
idx <- intersect(edge_septic, edge_healthy) # Get boolean index from septic that are in healthy
septic_filtered <- delete.edges(septic_graph, E(septic_graph)[idx])

blue <- "#33CCFF"
V(septic_filtered)$color <- "#999999"
lung <- V(septic_filtered)$name %in% c("Burkholderia", "Pandoraea", "Prevotella", "Aeromonas", "Enterobacter", "Stenotrophomonas", "Bacillus")
gut <- V(septic_filtered)$name %in% c("Proteus", "Enterococcus", "Lachnoclostridium", "Salmonella", "Shigella", "Blautia")
oral <- V(septic_filtered)$name %in% c("Bifidobacterium", "Veillonella", "Tannerella")
V(septic_filtered)[lung]$color <- "#33CCFF"
V(septic_filtered)[gut]$color <- "#CC9966"
V(septic_filtered)[oral]$color <- "#9933CC"
layout <- layout_(septic_filtered, with_dh(weight.edge.lengths = edge_density(septic_filtered)/500))

par(mar=c(0,0,0,0))
png(file="../results/septic_corrected_colored_network.png",width = 10, height = 8, units = 'in', res=1200)
plot(septic_filtered,
     vertex.frame.color = "white",
     vertex.size = 8,
     vertex.label.color = "black",
     vertex.label.cex = 1.1,
     vertex.label.font = 1.8,
     vertex.label.family = "Arial",
     vertex.label.dist = 0.2,
     edge.width = E(septic_filtered)$weight * width,
     edge.curved = 0.2,
     layout = layout)

legend("topright", legend=c("Oral", "Lung", "Gut"), 
       col=c("#9933CC", "#33CCFF", "#CC9966"),
       border = "black",
       pch = 19, cex = 1.2)

dev.off()

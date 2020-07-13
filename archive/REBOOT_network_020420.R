rm(list=ls())
#library(devtools)
#install_github("hallucigenia-sparsa/seqgroup")
require('seqgroup')
setwd('~/git_repos/Polymicrobial-Signature-of-Sepsis')
df <- read.csv('datasets/karius_genus_raw.csv')
X <- df[, c('Enterobacter', 'Moraxella', 'Stenotrophomonas', 'Prevotella',
             'Escherichia', 'Burkholderia', 'Bacteroides', 'Campylobacter',
             'Enterococcus', 'Streptococcus')]

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
edges <- 10
pval <- 0.05
width <- 8

# Plot septic graph
septic_graph <- barebonesCoNet(septic,
                               norm = T,
                               methods = c("spearman"),
                               init.edge.num = edges,
                               pval.cor = F,
                               pval.T = pval,
                               renorm = T,
                               permutandboot = T,
                               iters = 1000,
                               bh = T, 
                               plot = F,
                               keep.filtered = F,
                               verbose = T)

par(mar = c(0,0,0,0))
png(file = "results/decontam/sepsis_colored_network.png",width = 10, height = 8, units = 'in', res=300)

# Layout
color_weight <- E(septic_graph)$weight
color_weight[color_weight == "red"] <- color_weight[color_weight == "red"] * -3
color_weight[color_weight == "green"] <- color_weight[color_weight == "green"] * 3 
l <- layout_with_fr(septic_graph, weights=color_weight)

plot(septic_graph,
     vertex.frame.color = "white",
     vertex.size = 8,
     vertex.label.color = "black",
     vertex.label.cex = 0.9,
     vertex.label.font = 1.8,
     vertex.label.family ="Arial",
     vertex.label.dist = 0.2,
     edge.width = E(septic_graph)$weight * width,
     edge.curved = 0,
     layout = l,
     main = 'Microbial Network for Septic Patients')

dev.off()

# Plot healthy network
healthy_graph <- barebonesCoNet(healthy,
                                norm = T,
                                methods = c("spearman"),
                                init.edge.num = edges,
                                pval.cor = F,
                                pval.T = pval,
                                renorm = T,
                                permutandboot = T,
                                iters = 1000,
                                bh = T, 
                                plot = F,
                                keep.filtered = F,
                                verbose = T)

# Layout
color_weight <- E(healthy_graph)$weight
color_weight[color_weight == "red"] <- color_weight[color_weight == "red"] * -3
color_weight[color_weight == "green"] <- color_weight[color_weight == "green"] * 3 
l2 <- layout_with_fr(healthy_graph, weights=color_weight)

par(mar = c(0,0,0,0))
png(file = "results/decontam/healthy_colored_network.png", width = 10, height = 8, units = 'in', res = 300)

plot(healthy_graph,
     vertex.frame.color = "white",
     vertex.size = 8,
     vertex.label.color = "black",
     vertex.label.cex = 0.9,
     vertex.label.font = 1.8,
     vertex.label.family = "Arial",
     vertex.label.dist = 0.2,
     edge.width = E(healthy_graph)$weight * width,
     edge.curved = 0,
     layout = l2,
     main='Microbial Network for Healthy Patients')

dev.off()

# Substract healthy correlations
edge_septic <- E(septic_graph)
edge_healthy <- E(healthy_graph)
to_remove <- graph.intersection(septic_graph, healthy_graph) # Get boolean index from septic that are in healthy
idx <- which(as_ids(edge_septic) %in% as_ids(E(to_remove)))
septic_filtered <- delete.edges(septic_graph, idx)
# blue <- "#33CCFF"
# V(septic_filtered)$color <- "#999999"
# lung <- V(septic_filtered)$name %in% c("Burkholderia", "Pandoraea", "Prevotella", "Aeromonas", "Enterobacter", "Stenotrophomonas", "Bacillus")
# gut <- V(septic_filtered)$name %in% c("Proteus", "Enterococcus", "Lachnoclostridium", "Salmonella", "Shigella", "Blautia")
# oral <- V(septic_filtered)$name %in% c("Bifidobacterium", "Veillonella", "Tannerella")
# V(septic_filtered)[lung]$color <- "#33CCFF"
# V(septic_filtered)[gut]$color <- "#CC9966"
# V(septic_filtered)[oral]$color <- "#9933CC"

# Layout
color_weight <- E(septic_filtered)$weight
color_weight[color_weight == "red"] <- color_weight[color_weight == "red"] * -3
color_weight[color_weight == "green"] <- color_weight[color_weight == "green"] * 3 
l3 <- layout_with_fr(septic_filtered, weights=color_weight)

par(mar=c(0,0,0,0))
png(file="results/decontam/septic_corrected_colored_network.png",width = 10, height = 8, units = 'in', res=1200)
plot(septic_filtered,
     vertex.frame.color = "white",
     vertex.size = 8,
     vertex.label.color = "black",
     vertex.label.cex = 1.1,
     vertex.label.font = 1.8,
     vertex.label.family = "Arial",
     vertex.label.dist = 0.2,
     edge.width = E(septic_filtered)$weight * width,
     edge.curved = 0,
     layout = l3)

# legend("topright", legend=c("Oral", "Lung", "Gut"), 
#        col=c("#9933CC", "#33CCFF", "#CC9966"),
#        border = "black",
#        pch = 19, cex = 1.2)

dev.off()

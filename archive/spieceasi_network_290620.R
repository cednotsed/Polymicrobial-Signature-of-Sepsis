# Simple decontamination networks
rm(list=ls())
#library(devtools)
#install_github("hallucigenia-sparsa/seqgroup")
setwd('~/git_repos/Polymicrobial-Signature-of-Sepsis')

df <- read.csv('datasets/kapusta_grumaz_karius_genus_raw.csv')
# meta <- read.csv('datasets/simple_decontam_pathogens.csv', header = F, stringsAsFactors = F)
# meta <- meta[, 1]

# df <- df[df$pathogen %in% c("none", "Escherichia coli"), ]

# X <- df[, meta]
X <- df[, 2:ncol(df)]

# Count number of samples with non-zero k-mer counts
get_number <- function(df) {
        test <- data.frame(df)
        test[test > 0] <- 1
        test[test == 0] <- 0
        test <- apply(test, 2, sum)
        return(test)
}

numbers <- get_number(X)
to_retain <- names(numbers)[numbers > nrow(df) / 2]
X <- X[, to_retain]

septic <- X[df$y == 'septic', ]
healthy <- X[df$y == 'healthy', ]
healthy <- sample(healthy, size = length(septic))

# Tranpose to get taxa as rows samples as columns
septic <- t(septic) 
healthy <- t(healthy)

# Plot Septic network
set.seed(69)
edges <- 50
pval <- 0.05
width <- 8

# Plot septic graph
require(seqgroup)
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
png(file = "results/simple_decontam/sepsis_simple_decontam_network.png",width = 10, height = 8, units = 'in', res=300)

l <- layout_with_fr(septic_graph, weights=E(septic_graph)$weight * 1)

plot(septic_graph,
     vertex.frame.color = "white",
     vertex.size = 8,
     vertex.label.color = "black",
     vertex.label.cex = 0.9,
     vertex.label.font = 1.3,
     vertex.label.family ="Arial",
     vertex.label.dist = 0.2,
     edge.width = E(septic_graph)$weight * width,
     edge.curved = 0.01,
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

l2 <- layout_with_fr(healthy_graph, weights=E(healthy_graph)$weight * 1)
par(mar = c(0,0,0,0))
png(file = "results/simple_decontam/healthy_simple_decontam_network.png", width = 10, height = 8, units = 'in', res = 300)

plot(healthy_graph,
     vertex.frame.color = "white",
     vertex.size = 8,
     vertex.label.color = "black",
     vertex.label.cex = 0.5,
     vertex.label.font = 1.8,
     vertex.label.family = "Arial",
     vertex.label.dist = 0.2,
     edge.width = E(healthy_graph)$weight * width,
     edge.curved = 0,
     layout = l2,
     main='Microbial Network for Healthy Patients')

dev.off()

# Subtract healthy correlations
edge_septic <- E(septic_graph)
edge_healthy <- E(healthy_graph)
to_remove <- graph.intersection(septic_graph, healthy_graph) # Get boolean index from septic that are in healthy
idx <- which(as_ids(edge_septic) %in% as_ids(E(to_remove)))
septic_filtered <- delete.edges(septic_graph, idx)

# Remove vertices without edges
bad.vs <- V(septic_filtered)[degree(septic_filtered) == 0] 
septic_filtered <- delete.vertices(septic_filtered, bad.vs)
# blue <- "#33CCFF"
# V(septic_filtered)$color <- "#999999"
# lung <- V(septic_filtered)$name %in% c("Burkholderia", "Pandoraea", "Prevotella", "Aeromonas", "Enterobacter", "Stenotrophomonas", "Bacillus")
# gut <- V(septic_filtered)$name %in% c("Proteus", "Enterococcus", "Lachnoclostridium", "Salmonella", "Shigella", "Blautia")
# oral <- V(septic_filtered)$name %in% c("Bifidobacterium", "Veillonella", "Tannerella")
# V(septic_filtered)[lung]$color <- "#33CCFF"
# V(septic_filtered)[gut]$color <- "#CC9966"
# V(septic_filtered)[oral]$color <- "#9933CC"
color_weight <- E(septic_filtered)$weight
color_weight[color_weight == "red"] <- color_weight[color_weight == "red"] * -3
color_weight[color_weight == "green"] <- color_weight[color_weight == "green"] * 3 
l3 <- layout_with_fr(septic_filtered, weights=color_weight)

par(mar=c(0,0,0,0))
png(file="results/simple_decontam/septic_corrected_simple_decontam_network.png",width = 10, height = 8, units = 'in', res=600)
plot(septic_filtered,
     vertex.frame.color = "white",
     vertex.size = 8,
     vertex.label.color = "black",
     vertex.label.cex = 0.7,
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



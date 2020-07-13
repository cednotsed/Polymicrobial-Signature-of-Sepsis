# Simple decontamination networks
rm(list=ls())
#library(devtools)
#install_github("hallucigenia-sparsa/seqgroup")
setwd('~/git_repos/Polymicrobial-Signature-of-Sepsis')

df <- read.csv('datasets/karius_genus_raw_maxi.csv')


# Get simple decontam feature space
meta <- read.csv('datasets/simple_decontam_pathogens.csv', header = F, stringsAsFactors = F)
meta <- meta$V1

X <- df[, 2:ncol(df)]
X <- X[, colnames(X) %in% meta]

require(SpiecEasi)
require(igraph)
require(Matrix)

get_graph <- function(X) {
  set.seed(10010)
  
  # Run SparCC
  sparcc.amgut <- sparcc(X)

  ## Define arbitrary threshold for SparCC correlation matrix for the graph
  sparcc.graph <- sparcc.amgut$Cor
  sparcc.graph[sparcc.graph < 0.2] <- 0
  diag(sparcc.graph) <- 0
  sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
  
  ## Create igraph object
  vertex.names <- setNames(seq(ncol(X)), colnames(X))
  ig.sparcc <- adj2igraph(sparcc.graph, vertex.attr = vertex.names)
  V(ig.sparcc)$name <- colnames(X)
  
  bad.vs <- V(ig.sparcc)[degree(ig.sparcc) == 0] 
  sparcc.filt <- delete.vertices(ig.sparcc, bad.vs)
  
  return(sparcc.filt)
}

# Get septic/healthy graphs
X_healthy <- X[df$y == "healthy", ]
X_septic <- X[df$y == "septic", ]

healthy_g <- get_graph(X_healthy)
septic_g <- get_graph(X_septic)

# Subtract healthy co-occurences
edge_septic <- E(septic_g)
to_remove <- graph.intersection(septic_g, healthy_g) # Get boolean index from septic that are in healthy
idx <- which(as_ids(edge_septic) %in% as_ids(E(to_remove)))
corrected_g <- delete.edges(septic_g, idx)

# Remove edgeless vertices
bad.vs <- V(corrected_g)[degree(corrected_g) == 0] 
corrected_g <- delete.vertices(corrected_g, bad.vs)

confirmed <- c('Escherichia', 'Streptococcus', 'Mycobacterium', 'Cytomegalovirus',
               'Staphylococcus', 'Proteus', 'Klebsiella', 'Pseudomonas',
               'Moraxella', 'Enterococcus', 'Enterobacter', 'Citrobacter',
               'Haemophilus', 'Fusobacterium', 'Salmonella', 'Serratia',
               'Aerococcus', 'Campylobacter', 'Lymphocryptovirus', 'Simplexvirus')

CR <- c('Leptotrichia', 'Clostridium', 'Veillonella', 'Bacteroides',
        'Klebsiella', 'Shewanella', 'Prevotella', 'Capnocytophaga',
        'Stenotrophomonas', 'Escherichia', 'Burkholderia', 'Enterococcus',
        'Psychrobacter', 'Enterobacter', 'Alphatorquevirus', 'Bacillus',
        'Agrobacterium')

# Annotate vertex colors
V(corrected_g)$color <- "#999999"
confirmed_v <- V(corrected_g)$name %in% confirmed
CR_v <- V(corrected_g)$name %in% CR
CR_confirmed_v <- V(corrected_g)$name %in% intersect(confirmed, CR)
V(corrected_g)[confirmed_v]$color <- "red"
V(corrected_g)[CR_v]$color <- "blue"
V(corrected_g)[CR_confirmed_v]$color <- "purple"


# Plot
png(file = "results/sparCC_networks.png", 
    width = 10, 
    height = 8,
    units = 'in',
    res = 300)

clust <- cluster_edge_betweenness(corrected_g)
clust$membership[clust$membership != 1] <- NA
clust$membership[clust$names %in% c("Helicobacter", "Haemophilus")] <- 1
plot(corrected_g, 
     margin = c(0, -0.5, 0.5, 0),
     layout = layout.fruchterman.reingold(corrected_g), 
     vertex.size = 5,
     vertex.label.color = "black",
     vertex.frame.color = NA,
     vertex.label.dist = 1,
     edge.width = E(corrected_g)$weight * 10,
     mark.groups = communities(clust), 
     mark.col = c("#66FFFF"),
     mark.border = c("#66FFFF"))

legend("topright", legend=c("Confirmed pathogen", 
                            "In CR feature space", 
                            "Confirmed pathogen and in CR feature space"),
       col=c("red", "blue", "purple"),
       border = "black",
       pch = 19, cex = 1.2)

text(-0.4, 0.2, "Oral commensals", font = 2)

dev.off()


# Simple decontamination networks
rm(list=ls())
setwd('~/git_repos/Polymicrobial-Signature-of-Sepsis')

df <- read.csv('datasets/kapusta_grumaz_karius_genus_raw.csv')


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

CR <- c('Morganella', 'Blautia', 'Raoultella', 'Cellulosimicrobium',
        'Campylobacter', 'Alloprevotella', 'Megasphaera', 'Bacteroides',
        'Shewanella', 'Salmonella', 'Citrobacter', 'Cellulomonas',
        'Stenotrophomonas', 'Oerskovia', 'Enterobacter', 'Cupriavidus',
        'Rhodococcus', 'Clostridioides', 'Klebsiella', 'Pandoraea',
        'Cronobacter')

# Annotate vertex colors
V(corrected_g)$color <- "#999999"

# Plot
png(file = "results/sparCC_networks_pooled.png", 
    width = 10, 
    height = 8,
    units = 'in',
    res = 300)

clust <- cluster_louvain(corrected_g)
clust$membership[clust$membership != 4] <- NA
clust$membership[clust$names %in% c("Parabacteroides", "Alistipes")] <- NA
plot(corrected_g, 
     margin = c(0, -0.5, 0.5, 0),
     layout = layout.fruchterman.reingold(corrected_g, 
                                          weight = E(corrected_g)$weight * 1 / 2000, 
                                          niter = 1000), 
     vertex.size = 2,
     vertex.label.color = "black",
     vertex.frame.color = NA,
     vertex.label.dist = 0,
     vertex.label.cex = 0.6,
     edge.width = E(corrected_g)$weight * 10,
     mark.groups = communities(clust), 
     mark.col = c("#66FFFF"),
     mark.border = c("#66FFFF"),
     asp = 0,
     xlim = c(-1.5, 1.1),
     ylim = c(-1.1, 1)
)

text(0.9, 1, "Oral commensals", font = 2)

dev.off()


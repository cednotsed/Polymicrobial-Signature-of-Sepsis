# Simple decontamination networks
rm(list=ls())
#library(devtools)
#install_github("hallucigenia-sparsa/seqgroup")
setwd('~/git_repos/Polymicrobial-Signature-of-Sepsis')

df <- read.csv('datasets/kapusta_grumaz_karius_genus_raw.csv')
meta <- read.csv('datasets/simple_decontam_pathogens.csv', header = F, stringsAsFactors = F)
meta <- meta$V1

X <- df[, 2:ncol(df)]
X <- X[, colnames(X) %in% meta]
# Count number of samples with non-zero read counts
get_number <- function(df) {
  test <- data.frame(df)
  test[test > 0] <- 1
  test[test == 0] <- 0
  test <- apply(test, 2, sum)
  return(test)
}

numbers <- get_number(X)

names(numbers)[numbers < nrow(X) / 3]

library(SpiecEasi)
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

library(igraph)
bad.vs <- V(ig.sparcc)[degree(ig.sparcc) == 0] 
sparcc.filt <- delete.vertices(ig.sparcc, bad.vs)

l1 <- layout.fruchterman.reingold(sparcc.filt)

plot(sparcc.filt, 
     layout=l1, 
     vertex.size=2, 
     edge.width = E(ig.sparcc)$weight * 10)


###This file is used to process the data
###10 variables from topic modelling
sina_content <- read.delim("./Real\ network\ analysis/content.txt", header = FALSE, sep = " ")[,1:10]

###Use edge information to create adjacency matrix
sina_edge <- read.delim("./Real\ network\ analysis/edge.txt", head = FALSE, sep = "")

###create a directed graph
n <- nrow(sina_content)
T <- ncol(sina_content)
sina_A <- matrix(0, nrow = n, ncol = n)
for (i in 1:nrow(sina_edge)){
  sina_A[sina_edge[i,1],sina_edge[i,2]] = 1
}





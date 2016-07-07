require(igraph)
setwd("~/Documents/Network Science and Applications/project")
# load("Data/enron.RData")
# save.image(file = "Data/enron.RData")

#### load the graph ####
setwd("~/Documents/Network Science and Applications/project")
g <- read.graph("Data/email-Enron.txt", format = "edgelist", directed = FALSE)
g <- simplify(g) # delete multiple edges
V(g) # number of vertices 36692
E(g) # number of edges 183831

#### extract the subgraph for the giant component ####
decomposed_g <- decompose.graph(g) # 1065 components
giant_g <- decomposed_g[[which.max(sapply(decomposed_g, vcount))]]

V(giant_g) # 33696, 92%
E(giant_g) # 180811, 98%
average.path.length(giant_g) # 4.025164, which is even shorter than 6
diameter(giant_g) # 13

#### plot network graphs ####
# whole graph
png("Images/whole_graph.png")
plot.igraph(g, vertex.size = 5, vertex.label = NA, vertex.color="SkyBlue2")
dev.off()

# giant component
png("Images/giant_graph.png")
plot.igraph(giant_g, vertex.size = 5, vertex.label = NA, vertex.color="SkyBlue2")
dev.off()


#### degree ####
# fit power-law distribution
d <- degree(giant_g)
mean(d) # 10.7319
pl_fit <- power.law.fit(d)
pl_fit$alpha # 1.98
pl_fit$xmin # 5

# log-log plot
pdf("Images/degree_distribution.pdf")
par(mai=c(1, 1, .5, .5))
plot(degree.distribution(giant_g), log = "xy", xlab = "log(degree)", ylab = "log(frequency)", cex.lab = 1.5)
dev.off()

#### centrality ####
# closeness centrality
clo_centr <- closeness(giant_g)

# betweeness centrality
bet_centr <- betweenness(giant_g, directed = FALSE)
edge_bet_centr <- edge_betweenness(giant_g, directed = FALSE)

# eigenvector centrality
eigen_centr <- evcent(giant_g)$vector

#### clustering coefficient ####
clu_coef <- transitivity(giant_g, type = "average") # 0.71

#### detect communities ####
# louvain method
enron_communities <- cluster_louvain(giant_g)
community_member <- membership(enron_communities)
community_size <- sizes(enron_communities) # 239 communities
modularity(enron_communities) # 0.5969575
# sort(community_size, decreasing = TRUE)

# histogram of community size
pdf("Images/hist_community_louvain.pdf")
par(mai=c(1, 1, .5, .5))
hist(log(community_size), xlab = "log(size)", cex.lab = 1.5, main = "")
dev.off()

# label propagation
enron_communities2 <- cluster_label_prop(giant_g)
community_member2 <- membership(enron_communities2)
community_size2 <- sizes(enron_communities2) # 903 communities
modularity(enron_communities2) # 0.3136913
# sort(community_size2, decreasing = TRUE)

# histogram of community size
pdf("Images/hist_community_label.pdf")
par(mai=c(1, 1, .5, .5))
hist(log(community_size2), xlab = "log(size)", cex.lab = 1.5, main = "")
dev.off()

#### simulation ####
# find the optimum path to the target
find.optimum.path <- function(g, member, from, to, k = 10) {
  source("R/k_shortest_paths.R")
  suppressWarnings(spaths <- k.shortest.paths(g, from = from, to = to, k = k))
  shortest_dist <- spaths[[1]]$dist
  vertex_pass <- list()
  same_community <- rep(0, k)
  score <- rep(0, k)
  
  list_sp <- list()
  for (i in 1:10) {
    list_sp[[i]] <- as.vector(spaths[[i]]$vert[[1]])
  }
  rep <- length(unique(list_sp))
  
  if (shortest_dist > 4) {
    reachable = 0
    for (i in 1:k) {
      if (spaths[[i]]$dist == shortest_dist) {
        vertex_pass[[i]] <- as.vector(spaths[[i]]$vert[[1]])[2:5]
        same_community[i] <- sum(member[vertex_pass[[i]]] == member[to])  
      } else {
        vertex_pass[[i]] <- vector(mode = "integer")
      }
    }
  } else {
    reachable = 1
    for (i in 1:k) {
      if (spaths[[i]]$dist == shortest_dist) {
        vertex_pass[[i]] <- as.vector(spaths[[i]]$vert[[1]])[- c(1, shortest_dist + 1)]
        same_community[i] <- sum(member[vertex_pass[[i]]] == member[to])
      } else {
        vertex_pass[[i]] <- vector(mode = "integer")
      }
    }
  }
  for (i in which(same_community == max(same_community))) {
    score[i] <- sum(eigen_centr[vertex_pass[[i]]])
  }
  list(reachable = reachable, rep,
       vert = as.vector(spaths[[which.max(score)]]$vert[[1]]))
}

find.optimum.path2 <- function(g, member, from, to, k = 10) {
  source("R/k_shortest_paths.R")
  suppressWarnings(spaths <- k.shortest.paths(g, from = from, to = to, k = k))
  shortest_dist <- spaths[[1]]$dist
  vertex_pass <- list()
  same_community <- rep(0, k)
  score <- rep(0, k)
  
  list_sp <- list()
  for (i in 1:10) {
    list_sp[[i]] <- as.vector(spaths[[i]]$vert[[1]])
  }
  rep <- length(unique(list_sp))
  
  if (shortest_dist > 4) {
    reachable = 0
    for (i in 1:k) {
      if (spaths[[i]]$dist == shortest_dist) {
        vertex_pass[[i]] <- as.vector(spaths[[i]]$vert[[1]])[1:5]
        for (j in 1:4) {
          same_community[i] <- same_community[i] + 
            (member[vertex_pass[[i]][j]] == member[vertex_pass[[i]][j + 1]])
        }
      } else {
        vertex_pass[[i]] <- vector(mode = "integer")
      }
    }
  } else {
    reachable = 1
    for (i in 1:k) {
      if (spaths[[i]]$dist == shortest_dist) {
        vertex_pass[[i]] <- as.vector(spaths[[i]]$vert[[1]])
        for (j in 1:shortest_dist) {
          same_community[i] <- same_community[i] + 
            (member[vertex_pass[[i]][j]] == member[vertex_pass[[i]][j + 1]])
        }
      } else {
        vertex_pass[[i]] <- vector(mode = "integer")
      }
    }
  }
  for (i in which(same_community == max(same_community))) {
    score[i] <- sum(eigen_centr[vertex_pass[[i]]])
  }
  list(reachable = reachable, rep,
       vert = as.vector(spaths[[which.max(score)]]$vert[[1]]))
}

vertex_seq <- 1:33696

# Louvain method
set.seed(1128)
for (i in 1:10) {
  vertex_begin <- sample(vertex_seq, 1)
  community_index <- which(community_member == community_member[vertex_begin])
  vertex_end <- sample(vertex_seq[- community_index], 1)
  optimum_path <- find.optimum.path(giant_g, community_member, vertex_begin, vertex_end)
  optimum_path2 <- find.optimum.path2(giant_g, community_member, vertex_begin, vertex_end)
}

# label propagation
set.seed(1128)
for (i in 1:10) {
  vertex_begin <- sample(vertex_seq, 1)
  community_index <- which(community_member2 == community_member2[vertex_begin])
  vertex_end <- sample(vertex_seq[- community_index], 1)
  optimum_path <- find.optimum.path(giant_g, community_member2, vertex_begin, vertex_end)
  optimum_path2 <- find.optimum.path2(giant_g, community_member2, vertex_begin, vertex_end)
}


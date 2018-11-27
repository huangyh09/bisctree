# Chech if a clonal configuration satisfies a perfect phylogeny
# Perfect phylogeney means a mutation only happens once and never recurrents,
# which is equavlent as finit site assumption. 

# Algorithm from Dan Gusfield, 1991. Efficient algorithm for inferring 
# evolutionary trees.
# Implemented by Yuanhua Huang, 10 June 2018, on airplane.

# Lamma: a configuration matrix has a perfect phylogeny if and only if for every
# pair of clone i an j, either they are disjoint or one cotains another.

#' Check phylogeny for clonal configuration matrix
#' 
#' If a perfect phylogeny exists, build it. Otherwise, return NA or its near
#' perfect phylogeny.
#' 
#' @param Config A binary matrix, N SNVs x K clones.
#' @export
#' 
phylogeny_maker <- function(Config){
  #TODO: remove variants that do not occur in any clone.
  if (any(is.na(Config)) || any(is.null(Config))) {
    stop("[phylogeny_maker] Error: NA or NULL in Config, not supported.")
  }
  # Convert bool value to binary
  Config[Config == TRUE] <- 1
  Config[Config == FALSE] <- 0
  if (any(Config != 0 && Config != 1)) {
    print(Config)
    stop("[phylogeny_maker] Error: Config only supports binary value.")
  }
  
  N <- nrow(Config)
  K <- ncol(Config)
  
  # Step 1: order SNVs by its degree 
  SNV_degree <- rowSums(Config) + Config %*% 2**(seq(0, K - 1)) / (2**K)
  SNV_idx <- order(SNV_degree,  decreasing = TRUE)
  SNV_idx_use <- c(SNV_idx[1])
  Mut_list <- list(c(SNV_idx[1]))
  for (i in seq_len(N)) {
    if (i == 1) {next}
    if (SNV_degree[SNV_idx[i]] != SNV_degree[SNV_idx[i - 1]]) {
      SNV_idx_use <- c(SNV_idx_use, SNV_idx[i])
      Mut_list[[length(SNV_idx_use)]] <- c(SNV_idx[i])
    } else {
      Mut_list[[length(SNV_idx_use)]] <- c(Mut_list[[length(SNV_idx_use)]], 
                                           SNV_idx[i])
    }
  }
  Config_mut <- Config[SNV_idx_use, ]
  N_mut <- nrow(Config_mut)
  
  # Step 2: find the nearest ancestral SNV
  # variants not in any clone: ancestor with 0
  is_tree <- TRUE
  Ancestor_mat <- matrix(0, nrow = N_mut, ncol = K)
  Ancestor_SNV <- rep(0, N_mut)
  Ancestor_tmp <- rep(0, K)
  for (i in seq_len(N_mut)) {
    for (j in seq_len(K)) {
      if (Config_mut[i,j] == 1) {
        Ancestor_mat[i,j] <- Ancestor_tmp[j]
        Ancestor_tmp[j] <- i
      }
    }
    
    Ancestor_SNV[i] <- max(Ancestor_mat[i,])
    if (sum(Config_mut[i,]) == 0) {
      next
    }
    if (Ancestor_SNV[i] != min(Ancestor_mat[i, which(Config_mut[i,] == 1)])) {
      is_tree <- FALSE
    }
  }
  
  if (is_tree) {
    tree <- make_phylo(Config_mut, Ancestor_SNV, Mut_list)
    tree[["P"]] <- NULL
    tree[["Z"]] <- Config
  } else {
    tree <- NULL
  }

  return_list <- list("tree" = tree, 
                      "is_tree" = is_tree,
                      "Mut_list" = Mut_list,
                      "Config_mut" = Config_mut,
                      "Ancestor_mat" = Ancestor_mat,
                      "Ancestor_SNV" = Ancestor_SNV)
  return_list
}



library(ape)
#' Create ape::phylo object from the mutation configuation in phylogeny function
#' 
#' All inputs are from phylogeny function outputs
#' 
#' @param Config_mut A N_mut x K binary matrix. Here a mutation means a set of 
#' mutations have the same clonal occurrence.
#' @param Ancestor_SNV A vector of N_mut variable, indicating the parent 
#' mutation
#' @param Mut_list A list of SNVs for each mutation set
#' @import ape
#' 
make_phylo <- function(Config_mut, Ancestor_SNV, Mut_list) {
  # Assume always rooted tree, i.e., min(Ancestor_SNV) == 0.
  #TODO: remove variants that do not occur in any clone.
  
  K <- ncol(Config_mut)
  N_mut <- nrow(Config_mut)
  edge.label <- c()
  
  mut_degree <- rowSums(Config_mut)
  edge <- matrix(0, nrow = 0, ncol = 2)
  #tip edge without mutation: inner_node (old), tip_node(new)
  tip_edge <- matrix(0, nrow = 0, ncol = 2)
  node.ids <- rep(NA, N_mut)
  for (k in seq_len(ncol(Config_mut))) {
    # base clone
    if (sum(Config_mut[,k]) == 0) {
      edge <- rbind(edge, c(K + 1, k))
      edge.label <- c(edge.label, "")
      next
    }
    # add tips; if no unique mutation, add an edge
    idx_use <- which(Config_mut[,k] == 1)
    idx_min <- idx_use[which.min(mut_degree[idx_use])]
    if (mut_degree[idx_min] == 1) {
      node.ids[idx_min] <- k
    } else{
      tip_edge <- rbind(tip_edge, c(idx_min, k))
    }
  }
  
  # update node ids
  idx_Nnode <- which(is.na(node.ids))
  Nnode <- length(idx_Nnode) + 1      #always having root node
  node.ids[idx_Nnode] <- order(Ancestor_SNV[idx_Nnode]) + K + 1

  # add edges for tips without unique mutations
  for (i in seq_len(nrow(tip_edge))) {
    edge.label <- c(edge.label, "")
    edge <- rbind(edge, c(node.ids[tip_edge[i,1]], tip_edge[i,2]))
  }
  # add edges for non-tip mutations
  SNV_info <- matrix(0, nrow = 0, ncol = 3)
  for (i in seq_len(length(Ancestor_SNV))) {
    if (Ancestor_SNV[i] == 0) {
      edge <- rbind(edge, c(K + 1, node.ids[i]))
    } else{
      edge <- rbind(edge, c(node.ids[Ancestor_SNV[i]], node.ids[i]))
    }
    edge.label <- c(edge.label, paste( length(Mut_list[[i]]), "SNVs"))
    # keep SNV info
    for (j in seq_len(length(Mut_list[[i]]))) {
      SNV_info <- rbind(SNV_info, c(Mut_list[[i]][j], edge[nrow(edge), ]))
    }
  }
  colnames(SNV_info) <- c("sna", "sna.st.node", "sna.ed.node")
  
  # order siblings
  sort_idx <- order(edge %*% c(1, 1/(max(edge) + 1)))
  edge <- edge[sort_idx, ]
  edge.label <- edge.label[sort_idx]
  # branch.ids order is wrong; please use plot_tree
  # branch.ids <- c(edge.label[1:K], "Root", edge.label[(K + 1):nrow(edge)])
  tree <- list("sna" = SNV_info,
               "edge" = edge,
               "Nnode" = Nnode, 
               "edge.label" = edge.label,
               "tip.label" = paste0("C", seq(ncol(Config_mut))))
  class(tree) <- "phylo"
  tree
}


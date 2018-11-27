# Infer evolution tree from bulk DNA-seq with importance sampling


#' Generate samples of Psi from a beta proposal distribution
#'
#' @param proposal A list with two matrices: alpha and beta, both with shape of
#' n_cluster x n_sample. They are the paramters of beta distribution.
#' @param n_iter An integer, the number of samples to generate.
#'
generate_Psi <- function(proposal, n_iter=1){
  Psi_all <- list()
  log_Lik <- rep(0, n_iter)
  for(s in seq_len(n_iter)){
    Psi_mat <- matrix(1, nrow=nrow(proposal$alpha), ncol=ncol(proposal$alpha))
    Psi_new <- matrix(stats::rbeta(n=Psi_mat, proposal$alpha, proposal$beta),
                      nrow=nrow(proposal$alpha), ncol=ncol(proposal$alpha))
    Psi_all[[s]] <- Psi_new
    log_Lik[s] <- sum(stats::dbeta(Psi_new, proposal$alpha, proposal$beta, log=T))
  }

  return_list <- list("log_Lik"=log_Lik, "Psi_all"=Psi_all)
  return_list
}


#' log likelihood based on binomial distribution
#'
#' @param A a NxM matrix of alternative read counts. N: number of SNV, M: number
#' of samples.
#' @param D a NxM matrix of total read counts, sum of alternative and reference.
#' @param Psi a KxM matrix of binomial parameters of K clones across M samples.
#' @param om a vector with length of K. The prior the prevalence of K clones.
#'
Binom_logLik <- function(A, D, Psi, om){
  # TODO: see if it can be perform via snv so, for loop on ncol(A)
  prob_sum <- rep(0.0, nrow(A))
  for(i in seq_len(nrow(A))){
    #bin_prob <- rowSums(dbinom(A[i,], D[i,], Psi) * om)
    for(k in seq_len(nrow(Psi))){
      prob_sum[i] <- prob_sum[i] + prod(dbinom(A[i,], D[i,], Psi[k,])) * om[k]
    }
  }
  log_Lik <- sum(log(prob_sum))
  log_Lik
}


#' Log prior of Psi based on Dirichlet distribution and Sum Condition
#'
#' @param Psi a matrix of allele prevalence of each clone in each sample
#' @param B a perfect phylonogy tree, square and binary, matrix.
#' @param mu a matrix of hyper-parameters of dirichlet distribution
#'
Dirichlet_logLik <- function(Psi, B_inv_list, mu){
  #Note, no need to Psi_tmp[1,] <- 1, since base clone is not included.

  idx_perm <- gtools::permutations(n=nrow(Psi), r=nrow(Psi))
  log_Lik <- matrix(-Inf, nrow=nrow(idx_perm), ncol=length(B_inv_list))

  # check Sum Condition
  for(n in seq_len(length(B_inv_list))){
    B_inv <- as.matrix(B_inv_list[[n]])
    for(i in seq_len(nrow(idx_perm))){
      Psi_tmp <- Psi[idx_perm[i,], ]
      Phi <- B_inv %*% Psi_tmp

      if (min(Phi) >= 0 & max(colSums(Phi)) <= 0.5){
        Phi <- t(Phi) / colSums(t(Phi))
        log_Lik[i, n] <- sum(log(gtools::ddirichlet(Phi/rowSums(Phi), mu)))
        #log_Lik[i, n] <- 0
      }
    }
  }
  log_Lik
}


#' Infer evolution tree from bulk sequencing
#'
#' @export
#'
bulk_tree_IS <- function(A, D, B_inv_list, proposal, n_iter=1000,
                            om=NULL, mu=NULL, ga=NULL){
  # TODO: 1) adaptive importance sampling, 2) check the sum condition: using
  # multiple Psi in permutations or just once. Note, for 5-clone tree 5, the
  # leave switch only count once.

  # check input
  A <- as.matrix(A)
  D <- as.matrix(D)
  proposal$alpha <- as.matrix(proposal$alpha)
  proposal$beta <- as.matrix(proposal$beta)

  # check the size of A, D, proposal, B_list
  if (nrow(A) != nrow(D) || ncol(A) != ncol(D) ||
      nrow(proposal$alpha) != nrow(proposal$alpha) ||
      ncol(proposal$alpha) != ncol(proposal$alpha) ||
      ncol(proposal$alpha) != ncol(A) ||
      nrow(proposal$alpha) != nrow(B_inv_list[[1]])){
    stop(paste0("A and D must have the same size;\n",
                "proposal$alpha and proposal$beta must have the same size;\n",
                "ncol(A) and ncol(proposal$alpha) must have the same length;\n",
                "nrow(B) and nrow(proposal$alpha) must have the same length."))
  }

  # check missing arguments (default: uniform prior)
  K <- nrow(B_inv_list[[1]])           # number of clones
  S <- length(B_inv_list)              # number of valide trees
  if(is.null(om)){om <- rep(1/K, K)}   # uniform prior for clustering C
  if(is.null(mu)){mu <- rep(1,   K)}   # uniform prior for dirichlet
  if(is.null(ga)){ga <- rep(1/S, S)}   # uniform prior for trees

  dat_samp <- generate_Psi(proposal, n_iter)
  lik_samp <- dat_samp$log_Lik
  Psi_samp <- dat_samp$Psi_all

  log_Lik_all <- matrix(-Inf, nrow=n_iter, ncol=S)
  for(i in seq_len(n_iter)){
    Psi <- Psi_samp[[i]]
    log_prior <- Dirichlet_logLik(Psi, B_inv_list, mu)
    log_prior <- apply(log_prior, MARGIN=c(2), max) #check max or sum
    log_Lik_all[i, ] <- log_prior + Binom_logLik(A, D, Psi, om) - lik_samp[i]
  }
  log_Lik_all
  #matrixStats::colLogSumExps(log_Lik_all) - log(nrow(log_Lik_all))
}


sciClone_proposal <-function(A, D, tree_dir=NULL){
  N <- nrow(A)
  M <- ncol(A)
  vaf_list <- {}
  for(j in seq_len(M)){
    vaf_tmp <- data.frame(chrom=rep("chr1", N), pos=seq_len(N),
                          ref=D[,j]-A[,j], alt=A[,j], vaf=A[,j]/D[,j]*100)
    vaf_list[[j]] <- vaf_tmp
  }
  samp_names <- paste0("sample", seq_len(M))
  sc = sciClone(vafs=vaf_list, sampleNames=samp_names)

  K <- ncol(sc@clust$mu/sc@clust$alpha) + 1
  pp_alpha <- t(sc@clust$mu/sc@clust$alpha)
  pp_beta <- t(sc@clust$nu/sc@clust$beta)

  # pp_alpha <- t(sc@clust$mu*sc@clust$beta)
  # pp_beta <- t(sc@clust$nu*sc@clust$alpha)
  print(sc@clust$cluster.means)
  print(t(pp_alpha / (pp_alpha+pp_beta)))

  proposal = list()
  proposal[["alpha"]] <- pp_alpha
  proposal[["beta"]] <- pp_beta

  if(is.null(tree_dir)){
    return_list <- list("sc_obj"=sc, "proposal"=proposal)
  }else{
    dat_tmp <- read.table(paste0(tree_dir, "/Matrix", K, ".txt"),
                          skip=1)
    B_list <- list()
    B_inv_list <- list()
    for(tt in seq_len(nrow(dat_tmp)/ncol(dat_tmp))){
      first_idx <- ncol(dat_tmp) * (tt-1) + 1
      last_idx <- ncol(dat_tmp) * tt
      B_list[[tt]] <- as.matrix(dat_tmp[first_idx:last_idx, ])
      B_list[[tt]] <- B_list[[tt]][2:K, 2:K]
      B_inv_list[[tt]] <- solve(B_list[[tt]])
    }
    return_list <- list("sc_obj"=sc, "proposal"=proposal,
                        "B_list"=B_list, "B_inv_list"=B_inv_list)
  }

  return_list
}

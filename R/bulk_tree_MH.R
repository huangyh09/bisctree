# infer phylogeny from bulk sequencing data

#' Log posterior of a pair of Psi and Tree
#'
#' based on Dirichlet distribution and Sum Condition
#' log likelihood based on binomial distribution
#'
#' @param A a NxM matrix of alternative read counts. N: number of SNV, M: number
#' of samples.
#' @param D a NxM matrix of total read counts, sum of alternative and reference.
#' @param Psi a KxM matrix of clonal fractions in M samples.
#' @param B a perfect phylonogy tree, square and binary, matrix.
#' @param mu a matrix of hyper-parameters of dirichlet distribution
#' @param om a vector with length of K. The prior of the prevalence of K clones.
#' @param ga a vector with length of S. The prior of the S valide tree.
#'
log_post <- function(A, D, Psi, B, C=NULL, om=NULL, mu=NULL, ga=NULL){
  # check input and the size of A, D, Psi, B
  A <- as.matrix(A)
  D <- as.matrix(D)
  if (nrow(A) != nrow(D) || ncol(A) != ncol(D) ||
      ncol(Psi) != ncol(A) || nrow(Psi) != nrow(B)) {
    stop(paste0("A and D must have the same size;\n",
                "A and Psi must have the same ncol;\n",
                "B and Psi must have the same nrow."))
  }

  # calculate prior and its uniform default
  K <- nrow(B)                         # number of clones
  if (is.null(om)) {om <- rep(1/K, K)}   # uniform prior for clustering C
  if (is.null(mu)) {mu <- rep(1,   K)}   # uniform prior for dirichlet
  # if(is.null(ga)){ga <- rep(1,   S)} # uniform prior for trees

  tree_logPrior <- 0
  Diri_logPrior <- 0
  # Diri_logPrior <- sum(log(gtools::ddirichlet(Psi, mu)))

  # calculate VAF and binomial loglikelihood

  ### TODO: prob_mat for variant clustering

  VAF = B %*% Psi / 2.0
  prob_sum <- rep(0.0, nrow(A))
  for (i in seq_len(nrow(A))) {
    if (!is.null(C)) {
      prob_sum[i] <- prod(dbinom(A[i,], D[i,], VAF[C[i],]))  #* om[k]
      next
    }

    #bin_prob <- rowSums(dbinom(A[i,], D[i,], Psi) * om)
    for (k in seq_len(K)) {  #no base clone
      prob_sum[i] <- prob_sum[i] + prod(dbinom(A[i,], D[i,], VAF[k,])) * om[k]
    }
  }
  bin_logLik <- sum(log(prob_sum))

  # return log posterior
  log_post <- tree_logPrior + Diri_logPrior + bin_logLik
  log_post
}


#' Probability desensity function of Logistic-normal distribution
#'
#' @param x a fractional vector with sum to 1: length of K
#' @param mu the mean of the normal distribution: length of K-1
#' @param cov the covariance of the normal distribution: K-1 x K-1
#' @param log the bool variable. If yes, return log pdf
#'
#' @example
#' logistic_normal_pdf(x=c(0.22,0.4), mu=log(c(0.2,0.4)/0.4),
#'                     cov = matrix(c(0.1,0,0,0.1), nrow=2, ncol=2))
logistic_normal_pdf <- function(x, mu, cov, log = FALSE){
  x <- c(x, 1 - sum(x))
  if (sum(x) < 0.9999 || sum(x) > 1.0001) {
    stop("It requires sum(x) is 1.")
  }
  K <- length(x)
  x_logit <- matrix(log(x[1:(K - 1)]/x[K]) - mu, nrow = K - 1)
  cov_det <- det(cov)
  cov_inv <- solve(cov)

  if (cov_det < 0) {
    stop("Error: det(cov) is negative!")}
  rt_pdf <- (-0.5*log(2*pi*cov_det) - sum(log(x)) -
               0.5*(t(x_logit) %*% cov_inv %*% x_logit))

  if (log == FALSE) {rt_pdf = exp(rt_pdf)}
  rt_pdf
}


#' Metroplis-Hasting algorithm to sample tree and clonal fractions.
#'
#' @export
bulk_tree_MH <- function(A, D, B_list, n_iter=1000, trans_mat=NULL,
                         sigma=0.1, collaps_C=TRUE){
  A <- as.matrix(A)
  D <- as.matrix(D)

  N <- nrow(A)
  M <- ncol(A)
  S <- length(B_list)
  K <- nrow(B_list[[1]])
  if (is.null(trans_mat)) {
    trans_mat <- matrix(1.0/S, nrow = S, ncol = S)}
  cov <- diag(K) * sigma

  if (collaps_C) { Cluster_var = NULL}

  P_all <- rep(0, n_iter)
  Psi_all <- list()
  tree_all <- rep(0, n_iter)
  accept_all <- rep(FALSE, n_iter)

  # random initialization (TODO: EM initalization)
  Phi_cur <- matrix(0, nrow = K, ncol = M)    #logit Psi
  Psi_cur <- t(t(exp(Phi_cur)) / (colSums(exp(Phi_cur)) + 1))
  tree_cur <- sample(S, size = 1)
  P_cur <- log_post(A, D, Psi_cur, B_list[[tree_cur]])

  # sampling
  samp_info <- matrix(0, nrow = n_iter, ncol = 8)
  for (i in seq(1, n_iter)) {
    #proposal a new state
    for (m in seq(0, M)) {
      Phi_new <- Phi_cur
      Psi_new <- Psi_cur
      tree_new <- tree_cur
      if (m == 0) { #& last_accept==T
        # if (i < 0) {
        #   tree_prob <- rep(0, S)
        #   for (ii in seq_len(length(B_list))) {
        #     tree_prob[ii] <- log_post(A, D, Psi_cur, B_list[[ii]])
        #   }
        #   tree_prob <- exp(tree_prob - max(tree_prob))
        #   tree_prob <- tree_prob / sum(tree_prob)
        #   trans_mat[,] <- t(t(trans_mat[,]) + tree_prob)
        #   tree_cur <- sample(S, 1, replace = FALSE, prob = tree_prob)
        #   tree_all[i] <- tree_cur
        #   next
        # }
        trans_mat <- trans_mat / rowMeans(trans_mat)
        tree_new <- sample(S, size = 1, prob = trans_mat[tree_cur,])
        Q_cur <- 0 #log(trans_mat[tree_new, tree_cur])
        Q_new <- 0 #log(trans_mat[tree_cur, tree_new])
      } else{
        Phi_new[,m] <- mvtnorm::rmvnorm(n = 1, Phi_cur[,m], cov)
        Psi_new[,m] <- exp(Phi_new[,m]) / (sum(exp(Phi_new[,m])) + 1)
        Q_cur <- logistic_normal_pdf(Psi_cur[,m], Phi_new[,m], cov, log = T)
        Q_new <- logistic_normal_pdf(Psi_new[,m], Phi_cur[,m], cov, log = T)
      }
      P_new <- log_post(A, D, Psi_new, B_list[[tree_new]])
      if (P_new > -0.05) {print(c(i, m, P_new))}

      #accept ratio
      alpha <- min(1, exp(P_new - P_cur + Q_cur - Q_new))
      threshold <- runif(1, min = 0, max = 1)
      if (alpha >= threshold) {
        P_all[i] <- P_new
        tree_all[i] <- tree_new
        Psi_all[[i]] <- Psi_new
        accept_all[i] <- accept_all[i] + 1

        P_cur <- P_new
        Psi_cur <- Psi_new
        Phi_cur <- Phi_new
        tree_cur <- tree_new
      }else{
        P_all[i] <- P_cur
        tree_all[i] <- tree_cur
        Psi_all[[i]] <- Psi_cur
      }
    }
    samp_info[i,] <- c(i, accept_all[i], tree_new, alpha, P_new, P_cur, Q_cur, Q_new)
  }
  #print(samp_info)
  return_list <- list(Psi_all = Psi_all, tree_all = tree_all, P_all = P_all,
                      accept_all = accept_all)
  return_list
}


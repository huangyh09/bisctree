## This file is for binomial distributions

#' EM algorithm for estimating binomial mixture model
#' 
#' 
#' @param k A vector of integers. number of success
#' @param n A vector of integers. number of trials
#' @param n_components A number. number of components
#' @param tol numeric(1), tolerance value for convergence between iterations
#' 
#' @return a list containing \code{p}, a vector of floats between 0 and 1 giving
#' the estimated success probability for each component, \code{psi}, estimated
#' fraction of each component in the mixture, and \code{prob}, the matrix of 
#' fitted probabilities of each observation belonging to each component.
#' 
#' @import stats
#' 
#' @export
#' 
#' @examples
#' n1 <- array(sample(1:30, 50, replace = TRUE))
#' n2 <- array(sample(1:30, 200, replace = TRUE))
#' k1 <- apply(n1, 1, rbinom, n = 1, p = 0.5)
#' k2 <- apply(n2, 1, rbinom, n = 1, p = 0.01)
#' RV <- mixBinom(c(k1, k2), c(n1, n2))
mixBinom <- function(k, n, n_components = 2, tol = 1e-03, min_iter=200,
                     max_iter=10000, uniform_prior=FALSE, n_initial=20){
  #TODO: it's very easy to reach local optima in mutliple samples
  k <- as.matrix(k)
  n <- as.matrix(n)
  if (nrow(k) != nrow(n) || ncol(k) != ncol(n)) {
    stop("n and k must be of the same size.\n")
  }
  S <- nrow(n)
  M <- ncol(n)
  
  ## Random initialzation on parameters
  #p_share <- 
  
  p <- matrix(stats::runif(n_components * M, 0.0, 1.0), nrow=n_components)
  
  if(uniform_prior){ 
    n_param <- n_components * M
    psi <- rep(1/n_components, n_components)
  }else{ 
    n_param <- n_components * (M+1)
    psi <- prop.table(stats::runif(n_components, 0.0, 1.0))
  }
  prob_mat <- t(stats::rmultinom(S, size=1, prob=psi))
  
  log_like <- -Inf
  log_like_new <- -1e-100
  
  best_p <- p
  best_psi <- psi
  best_prob_mat <- prob_mat
  best_log_like <- log_like_new
  run_gap <- as.integer(min_iter/n_initial)
  print(run_gap)
  
  ## Iterations
  for(run_time in seq_len(max_iter)){
    if(run_time > min_iter && (log_like_new - log_like) < tol){break}
    
    log_like <- log_like_new
    
    ## M-step
    for (j in seq_len(n_components)) {
      p[j,] <- rowSums(t(k)*prob_mat[,j]) / rowSums(t(n)*prob_mat[,j])
      if(uniform_prior==FALSE){ psi[j] <- sum(prob_mat[,j]) / S}
    }
    ### correct empty clusters
    # if(run_time < 50){
    #   p[which(is.na(p))] <- stats::runif(sum(is.na(p)), 0.0, 1.0)
    # }
    if(run_time < min_iter){
      p[which(is.na(p))] <- stats::runif(sum(is.na(p)), 0.0, 1.0)
    }
    
    
    ## E-step
    like_mat <- matrix(0, nrow=S, ncol=n_components)
    for (j in seq_len(n_components)) {
      for(m in seq_len(M)){
        lik_tmp <- stats::dbinom(k[,m], size = n[,m], prob = p[j,m], log = T)
        like_mat[,j] <- like_mat[,j] + lik_tmp
      }
    }
    prob_mat <- exp(like_mat - apply(like_mat,1,max))
    prob_mat <- t(t(prob_mat) * psi)
    #prob_mat[which(prob_mat < 1e-20)] <- 1e-20
    prob_mat <- prob_mat / rowSums(prob_mat)
    #prob_mat[which(is.na(prob_mat))] <- -Inf
    
    # logLikelihood
    log_like_new <- 0
    like_mat[which(is.na(like_mat))] <- -Inf
    for(i in seq_len(S)){
      log_like_tmp <- matrixStats::logSumExp(log(psi) + like_mat[i,])
      log_like_new <- log_like_new + log_like_tmp
    }
    
    ### TODO: multiple short initialization
    # if((run_time < min_iter) && ((run_time %% run_gap) == 0)){
    #   print(c(run_time, run_time %% run_gap))
    #   if(log_like_new > best_log_like){
    #     best_log_like <- log_like_new
    #     best_prob_mat <- prob_mat
    #     best_p <- p
    #     best_psi <- psi
    #   }
    #   log_like_new <- -1e-100
    #   prob_mat <- t(stats::rmultinom(S, size=1, prob=psi))
    #   p <- matrix(stats::runif(n_components * M, 0.0, 1.0), nrow=n_components)
    #   if(!uniform_prior){
    #     psi <- prop.table(stats::runif(n_components, 0.0, 1.0))
    #   }
    # }else if(run_time == min_iter){
    #   log_like_new <- best_log_like
    #   prob_mat <- best_prob_mat
    #   p <- best_p
    #   psi <- best_psi
    # }
  }
  BIC <- -2*log_like + n_param * (log(S) - log(2*pi))
  
  ## return values
  return_list <- list("p" = p, "psi" = psi, "prob" = prob_mat, "BIC" = BIC,
                      "log_like" = log_like)
  return_list
}


#' Predicted probability from learned binomial mixture model
#' 
#' @param k A vector of integers. number of success
#' @param n A vector of integers. number of trials
#' @param p a vector of binomial success probabilities
#' @param psi A float between 0 and 1. fraction of each component
#' 
#' @export
#' 
#' @examples
#' n1 <- array(sample(1:30, 50, replace = TRUE))
#' n2 <- array(sample(1:30, 200, replace = TRUE))
#' k1 <- apply(n1, 1, rbinom, n = 1, p = 0.5)
#' k2 <- apply(n2, 1, rbinom, n = 1, p = 0.01)
#' RV <- mixBinom(c(k1, k2), c(n1, n2))
#' prob <- predMixBinom(3, 10, RV$p, RV$psi)
predMixBinom <- function(k, n, p, psi) {
    if (length(p) != length(psi)) {
        stop("p and psi must be of the same length (number of mixture components)\n")
    }
    if (length(k) != length(n)) {
        stop("n and k must be of the same length (number of observations)\n")
    }
    prob_test <- matrix(nrow = length(k), ncol = length(p))
    for (j in seq_along(p))
        prob_test[, j] <- psi[j] * stats::dbinom(k, size = n, prob = p[j], 
                                                 log = FALSE) 
    prob_test <- prob_test / rowSums(prob_test)
    prob_test
}


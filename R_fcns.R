library(tidyverse)
library(dplyr)
library(magrittr)
library(reshape)
library(MASS)


###########################################
#####  Create a logger for the loops  #####
###########################################

logger <- function (msg, ...) {
  cat(as.character(Sys.time()), "-", sprintf(msg, ...))
}


###########################################
#####  Compute the trace of a matrix  #####
###########################################

trace <- function(A) {
  p <- dim(A)[1] 
  tr <- 0 
  
  for (j in 1:p) {
    diagonal <- A[j,j]
    tr <- tr + diagonal
  }
  return(tr[[1]])
}


##########################################################
#####  Specify the interval of ridge penalty values  #####
##########################################################

lambda_values <- function(Xs, q = 0.5, l_negative, l_positive, case) {
  X <- as.matrix(Xs)
  n <- nrow(X)
  p <- ncol(X)
  
  X_svd <- svd(X)
  d <- X_svd$d
  threshold <- quantile(d, probs = c(q))
  
  # lower bound for lambda
  d_tilde <- (1 - as.numeric(d < threshold[[1]]))*d
  min_dt <- min(d_tilde[which(d_tilde != 0)], na.rm = TRUE)
  lower <- -(sqrt(p/n) - 1)*(sqrt(p/n) - 1)*(min_dt)
  
  # set-up grid points for lambda
  if (lower < 0) {
    l_bdd_n <- min(lower,-sqrt(.Machine$double.eps))
    u_bdd_n <- max(lower,-sqrt(.Machine$double.eps))
    l_bdd_p <- min(upper,sqrt(.Machine$double.eps))
    u_bdd_p <- max(upper,sqrt(.Machine$double.eps))
    
    neg_lambdas <- seq(l_bdd_n, u_bdd_n, length = l_negative - 1)
    
    if (case == 1) {
      pos_lambdas <- seq(l_bdd_p, u_bdd_p, length = l_positive)
      lambdas <- sort(c(neg_lambdas, 0, pos_lambdas))  
    } else if (case == 2) {
      pos_step <- (u_bdd_n - l_bdd_n)/(l_negative - 1)
      pos_lambdas <- seq(l_bdd_p, u_bdd_p, by = 20*pos_step)
      lambdas <- sort(c(neg_lambdas, 0, pos_lambdas))
    } else if (case == 3) {
      pos_lambdas <- c(10^(seq(-9, 5, length = l_positive)))
      lambdas <- sort(c(neg_lambdas, 0, pos_lambdas))
    }
  } else {
    pos_lambdas <- seq(lower, upper, length = l_positive)
    lambdas <- sort(c(0, pos_lambdas))
  }
  
  n_lambda <- length(lambdas)
  output <- list(lambdas = lambdas, n_lambda = n_lambda, X = X,
                 d = d, d_tilde = d_tilde)
  return(output)
}


##########################################################
##### Compute the inverse of Tikhonov regularization #####
##########################################################

inv_Tikhonov <- function(Xs, lambdas_scalar) {
  X <- as.matrix(Xs)
  n <- nrow(X)
  p <- ncol(X)
  
  if (n > p) {
    TR <- as.matrix(t(X) %*% X + diag(n*lambdas_scalar, p))
    inv_TR <- ginv(TR, tol = sqrt(.Machine$double.eps))
    W <- as.matrix(inv_TR %*% t(X))
    L_lambda <- as.matrix(X %*% W)
  } else {
    TR_alt <- as.matrix(X %*% t(X) + diag(n*lambdas_scalar, n))
    inv_TR <- ginv(TR_alt, tol = sqrt(.Machine$double.eps))
    W <- as.matrix(t(X) %*% inv_TR)
    if (lambdas_scalar != 0) {
      L_lambda <- as.matrix(X %*% W)
    } else {
      L_lambda <- diag(1, n) - inv_TR
    }
  }

  output <- list(inv_TR = inv_TR, W = W, L_lambda = L_lambda)
  return(output)
}


##########################################################
#####  Create the covariance matrix for Simulations  #####
##########################################################

make_cov <- function(cov_type, p, rho = 0.1) {
  if (cov_type == "ID") {
    Sigma <- diag(1, nrow = p)
  } else if (cov_type == "CSYM") {
    Sigma <- matrix(rho, nrow = p, ncol = p)
    diag(Sigma) <- 1
  } else if (cov_type == "NCDIAG") {
    Sigma <- diag(log(2:(p+ 1)), nrow = p, ncol = p)
  } else if (cov_type == "AR") {
    Sigma <- diag(1, nrow = p, ncol = p)
    Sigma[1, ] <- rho^(0:(p - 1))
    for (j in 2:p) {
      Sigma[j, ] <- c(Sigma[1, j:p], Sigma[1, 1:(j - 1)])
    }
  }
  return(Sigma)
}


##########################################################
#####    Generate the data for simulation studies    #####
##########################################################

generate_data <- function(n, p, cov_X, rho_X, data_X, sim_X, sim_beta, 
                          mu, sigma, sparsity, seed) {
  if (sim_X == TRUE) {
    set.seed(seed)
    Sigma <- make_cov(cov_X, p, rho_X)
    X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigma)
    Xs <- scale(X, center = TRUE, scale = TRUE)
    
    if (sim_beta == TRUE) {
      beta <- as.matrix(rnorm(p, mu, sigma))
      if (sparsity > 0) beta[1:round(p * sparsity)] <- 0   
    } else {
      beta <- as.matrix(rep(mu, p))
    }
    
    stopifnot(ncol(Xs) == nrow(beta))
  } else {
    Xs <- as.matrix(data_X)
    if (sim_beta == TRUE) {
      set.seed(seed)
      beta <- as.matrix(rnorm(ncol(Xs), mu, sigma))
      if (sparsity > 0) beta[1:round(ncol(Xs) * sparsity)] <- 0  
    } else {
      beta <- as.matrix(rep(mu, ncol(Xs)))
    }
    stopifnot(ncol(Xs) == nrow(beta))
  }
  
  output <- list(Xs = Xs, beta = beta)
  return(output)
}


##########################################################
##### Compute the optimal lambda using GCV and LOOCV #####
##########################################################

lambda_optimal <- function(Xs, Y_perm, lambdas, n_perms) {
  start <- Sys.time()
  n_perms_p1 <- n_perms + 1
  n_lambda   <- length(lambdas)
  
  X <- as.matrix(Xs)
  Y <- as.matrix(Y_perm[[1]])
  n <- nrow(X)
  p <- ncol(X)
  
  gcv <- matrix(NA, nrow = n_perms_p1, ncol = n_lambda)
  loo <- matrix(NA, nrow = n_perms_p1, ncol = n_lambda)
  pct_rank <- array(NA, c(2, n_perms_p1, n_lambda))
  lambda_opt <- matrix(NA, nrow = n_perms_p1, ncol = 4)
  colnames(lambda_opt) <- c("GCV_idx", "lhat_GCV", "LOO_idx", "lhat_LOO")
  
  for (l in 1:n_lambda) {
    inv_TReg <- inv_Tikhonov(X, lambdas[l])
    L_lambda <- inv_TReg$L_lambda
    inv_TR <- inv_TReg$inv_TR
    
    for(m in 1:n_perms_p1) {
      beta_hat <- as.numeric(inv_TReg$W %*% Y_perm[[m]])
      
      gcv_num <- as.numeric(Y - X %*% beta_hat)
      gcv_den <- 1 - trace(L_lambda)/n
      loo_den <- diag(diag(1,n) - diag(L_lambda))
      
      gcv[m,l] <- sum((gcv_num/gcv_den)^2)/n 
      loo[m,l] <- sum((gcv_num/loo_den)^2)/n
      
      if (is.na(gcv[m,l] || is.na(loo[m,l])) == TRUE) {
        gcv_num <- as.matrix(inv_TR %*% Y)
        gcv_den <- trace(inv_TR)/n
        loo_den <- diag(inv_TR)
        
        gcv[m,l] <- sum((gcv_num/gcv_den)^2)/n
        loo[m,l] <- sum((gcv_num/loo_den)^2)/n
      } 
      
      if (is.nan(gcv[m,l] || is.nan(loo[m,l])) == TRUE) {
        gcv_num <- as.matrix(inv_TR %*% Y)
        gcv_den <- trace(inv_TR)/n
        loo_den <- diag(inv_TR)
        
        gcv[m,l] <- sum((gcv_num/gcv_den)^2)/n
        loo[m,l] <- sum((gcv_num/loo_den)^2)/n
      }
    }
  }
  
  for (m in 1:n_perms_p1) {
    l_idx <- 1:n_lambda
    idx_gcv <- which(gcv[m,] == min(gcv[m,], na.rm = TRUE))
    idx_loo <- which(loo[m,] == min(loo[m,], na.rm = TRUE))
    
    if(length(idx_gcv) == 1) {
      gcv_min <- l_idx[idx_gcv]
    } else if(length(idx_gcv) > 1) {
      gcv_min <- min(idx_gcv)
    } else {
      gcv_min <- NA
    }
    
    if(length(idx_loo) == 1) {
      loo_min <- l_idx[idx_loo]
    } else if(length(idx_loo) > 1) {
      loo_min <- min(idx_loo)
    } else {
      loo_min <- NA
    }
    
    if(is.na(gcv_min) == TRUE) {
      lambda_opt_gcv <- NA
    } else if(is.na(loo_min) == TRUE) {
      lambda_opt_loo <- NA
    } else {
      lambda_opt_gcv <- lambdas[gcv_min]
      lambda_opt_loo <- lambdas[loo_min]
    }
    
    lambda_opt[m,1] <- gcv_min
    lambda_opt[m,2] <- lambda_opt_gcv
    lambda_opt[m,3] <- loo_min
    lambda_opt[m,4] <- lambda_opt_loo
  }
  
  for(l in 1:n_lambda) {
    pct_rank_gcv <- percent_rank(gcv[,l])
    pct_rank_loo <- percent_rank(loo[,l])
    pct_rank[1,,l] <- pct_rank_gcv
    pct_rank[2,,l] <- pct_rank_loo
  }
  
  output <- list(lambdas = lambdas, gcv = gcv, loo = loo,
                 lambda_opt = lambda_opt, pct_rank = pct_rank)
  return(output)
}

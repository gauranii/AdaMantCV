source("R_fcns.R")


###########################################################
##### Adaptive Mantel with Cross Validation Algorithm #####
###########################################################

adamant_cv <- function(Xs, Ys, q, l_negative, l_positive, 
                       case, n_perms, seed) {
  start <- Sys.time()
  X <- as.matrix(Xs)
  Y <- as.matrix(Ys)
  n <- nrow(X)
  p <- ncol(X)
  
  stopifnot(nrow(X) == nrow(Y))
  
  #########################################################
  ##### Specify the range of ridge penalty parameters #####
  #########################################################
  
  lv <- lambda_values(X, q, l_negative, l_positive, case)
  lambdas <- lv$lambdas
  d_tilde <- lv$d_tilde
  
  n_perms_p1 <- n_perms + 1
  n_lambda   <- length(lambdas)
  
  ##############################################
  ##### Generate permutations of rows of Y #####
  ##############################################
  
  X_svd <- svd(X)
  U <- X_svd$u
  
  set.seed(seed)
  Y_perm <- prepend(map(1:n_perms, function(b) {Y[sample(1:n),]}), list(Y))
  Z_perm <- map(Y_perm, ~ crossprod(U, .x))
  
  #########################################################################
  ##### Compute the test statistic T and percent rank P within lambda #####
  #########################################################################
  
  weights <- map(lambdas, ~ pmax(0, d_tilde^2/(.x + d_tilde^2)))
  if (0 %in% lambdas) {
    weights[[which(lambdas == 0)]] <- rep(1, length(d_tilde))
  } else if (Inf %in% lambdas) {
    weights[[which(lambdas == Inf)]] <- d_tilde^2
  }
  
  out <- tibble(lambda = rep(lambdas, each = n_perms_p1),
                b = rep(1:(n_perms_p1), n_lambda), 
                MT = NA, gcv = NA, loo = NA, Pval_gamma = NA)
  
  for(l in 1:n_lambda) {
    out$MT[(((l - 1)*n_perms_p1) + 1):(l*(n_perms_p1))] <-
      unlist(map(Z_perm, ~crossprod(weights[[l]], .x^2)))
  }
  
  out <- out %>%
    group_by(lambda) %>%
    mutate(Pval_lambda = 1 - percent_rank(MT))
  
  test_stat <- matrix(out$MT, nrow = n_perms_p1, ncol = n_lambda)
  pct_rank_ts <- matrix(out$Pval_lambda, nrow = n_perms_p1, ncol = n_lambda)
  
  est_MOM <- matrix(NA, nrow = 2, ncol= n_lambda)
  P_gamma <- matrix(NA, nrow = n_perms_p1, ncol= n_lambda)
  
  for(l in 1:n_lambda) {
    est_MOM[1,l] <- mean(test_stat[2:n_perms_p1,l], na.rm = TRUE)
    est_MOM[2,l] <- var(test_stat[2:n_perms_p1,l], na.rm = TRUE)
  }
  
  mean_MOM <- est_MOM[1,]
  var_MOM  <- est_MOM[2,]
  beta_MOM  <- mean_MOM/var_MOM
  alpha_MOM <- mean_MOM^2/var_MOM
  
  for(l in 1:n_lambda) {
    for(m in 1:n_perms_p1) {
      P_gamma[m,l] <- pgamma(test_stat[m,l], alpha_MOM[l], rate = beta_MOM[l], 
                             lower.tail = FALSE)
    }
  }
  
  ##########################################################
  ##### Compute the optimal lambda using GCV and LOOCV #####
  ##########################################################
   
  l_opt <- lambda_optimal(Xs, Y_perm, lambdas, n_perms)
  out$gcv <- as.vector(l_opt$gcv)
  out$loo <- as.vector(l_opt$loo)
  out$Pval_gamma <- as.vector(P_gamma)
  out_gcv <- out_loo <- out_gam <- out
  
  compute_gcv <- l_opt$gcv
  compute_loo <- l_opt$loo
  summary_cv  <- l_opt$lambda_opt
  pct_rank_cv <- l_opt$pct_rank
  
  ####################################################################
  ##### Compute adaptive Mantel Test P value using GCV and LOOCV #####
  ####################################################################
  
  Pval_gcv <- out_gcv %>%
    group_by(b) %>%
    summarize(Pval_b_gcv = Pval_lambda[which.min(gcv)]) %>%
    transmute(Pval_adapt_gcv = percent_rank(Pval_b_gcv), b = b) %>%
    filter(b == 1) %$%
    Pval_adapt_gcv
  
  Pval_loo <- out_loo %>%
    group_by(b) %>%
    summarize(Pval_b_loo = Pval_lambda[which.min(loo)]) %>%
    transmute(Pval_adapt_loo = percent_rank(Pval_b_loo), b = b) %>%
    filter(b == 1) %$%
    Pval_adapt_loo
  
  Pval_gam_gcv <- out_gam %>%
    group_by(b) %>%
    summarize(Pval_g_gcv = Pval_gamma[which.min(gcv)]) %>%
    transmute(Pval_g_adapt_gcv = percent_rank(Pval_g_gcv), b = b) %>%
    filter(b == 1) %$%
    Pval_g_adapt_gcv
  
  Pval_gam_loo <- out_gam %>%
    group_by(b) %>%
    summarize(Pval_g_loo = Pval_gamma[which.min(loo)]) %>%
    transmute(Pval_g_adapt_loo = percent_rank(Pval_g_loo), b = b) %>%
    filter(b == 1) %$%
    Pval_g_adapt_loo
  
  lambda_gcv <- out_gcv %>% 
    group_by(b) %>% 
    slice(which.min(gcv)) %$%
    lambda
  
  lambda_loo <- out_loo %>%
    group_by(b) %>% 
    slice(which.min(loo)) %$%
    lambda
  
  Pval_b <- as.matrix(cbind(
    out_gcv %>%
      group_by(b) %>%
      summarize(Pval_b_gcv = Pval_lambda[which.min(gcv)]),
    out_loo %>%
      group_by(b) %>%
      summarize(Pval_b_loo = Pval_lambda[which.min(loo)]),
    out_gam %>%
      group_by(b) %>%
      summarize(Pval_g_gcv = Pval_gamma[which.min(gcv)]),
    out_gam %>%
      group_by(b) %>%
      summarize(Pval_g_loo = Pval_gamma[which.min(loo)])))
  Pval_b <- Pval_b[,c(-3,-5,-7)]
  Pval_pct_rank <- cbind(Pval_b, percent_rank(Pval_b[,2]), 
                         percent_rank(Pval_b[,3]), percent_rank(Pval_b[,4]), 
                         percent_rank(Pval_b[,5]), as.numeric(lambda_gcv),
                         as.numeric(lambda_loo))
  end <- Sys.time()
  print(end - start)
  
  output <- list(lambdas = lambdas, test_stat = test_stat, out = out,
                 pct_rank_ts = pct_rank_ts, gcv = compute_gcv, 
                 loo = compute_loo, summary_cv = summary_cv,
                 pct_rank_cv = pct_rank_cv, P_gamma = P_gamma,
                 Pval_gcv = Pval_gcv, Pval_loo = Pval_loo, 
                 Pval_gam_gcv = Pval_gam_gcv, Pval_gam_loo = Pval_gam_loo,
                 Pval_b = Pval_b, Pval_pct_rank = Pval_pct_rank)
  return(output)
}

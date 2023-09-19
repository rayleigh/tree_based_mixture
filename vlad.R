library(MCMCpack)

calc_extension_param <- function(K, alpha = rep(1, K), num_iter = 2000) {
  draws <- MCMCpack::rdirichlet(num_iter, alpha)
  v_list <- kmeans(draws, K)$centers
  gamma_val = 0
  for (k in 1:K) {
    gamma_val = gamma_val + norm(v_list[k,] - 1 / K, type = "2")
  }
  sqrt(K^2 - K) / gamma_val
}


vlad <- function(data, K) {
  data_center = colMeans(data)
  centered_data <- scale(data, scale = F)
  svd_info <- svd(centered_data)
  u_list <- svd_info$u[, 1:(K - 1)]
  v_list <- kmeans(u_list, K)$centers
  topics_m <- matrix(0, nrow = K, ncol = ncol(data))
  gamma_val = calc_extension_param(K)
  for (i in 1:K) {
    topics_m[i,] <- data_center + 
      gamma_val * svd_info$v[, 1:(K - 1)] %*% (svd_info$d[1:(K - 1)] * v_list[i,])
  }
  return(topics_m)
}


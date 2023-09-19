get_r_k <- function(data, centroid) {
  max(apply(data, 1, function(row) {norm((centroid - row), type = "2")}))
}

proj_on_s <- function(beta_k, doc, K, remain_ind, dist = F) {
  thetas <- rep(0, K)
  if (length(remain_ind) == 1) {
    thetas[remain_ind[1]] = 1
    if (dist) {
      return(norm((doc - beta_k[, remain_ind[1]]), type = "2"))
    } else {
      return(thetas)
    }
  }
  beta_0 <- beta_k[, remain_ind[1]]
  center_beta_k <- apply(beta_k[, remain_ind[-1], drop = F], 2, function(col) {col - beta_0})
  alphas <- solve(t(center_beta_k) %*% center_beta_k) %*% t(center_beta_k) %*% (doc - beta_0)
  if (sum(alphas) <= 1 && all(alphas >= 0)) {
    if (dist) {
      return(norm((doc - beta_0 - center_beta_k %*% alphas), type = "2"))
    } else {
      thetas[remain_ind] <- c(1 - sum(alphas), alphas)
      return(thetas)
    }
  } else if (any(alphas < 0)) {
    remain_ind <- c(remain_ind[1], remain_ind[which(alphas > 0) + 1])
    return(proj_on_s(beta_k, doc, K, remain_ind, dist))
  } else {
    remain_ind <- remain_ind[-1]
    return(proj_on_s(beta_k, doc, K, remain_ind, dist))
  }
}

adj_center <- function(c, center, gen_center) {
  new_center <- gen_center + c * (center - gen_center)
  new_center[new_center < 0] <- 0
  new_center <- new_center / sum(new_center)
}

get_proj_dist_for_centers <- function(r, k, k_center, gen_center, all_centers, data) {
  all_centers[k, ] <- adj_center(r, k_center, gen_center)
  K <- nrow(all_centers)
  symplex_dist <- mean(apply(data, 1, function(obs) {proj_on_s(t(all_centers), obs, K, 1:K, T)}))
}

tune_center_for_k <- function(k, k_center, gen_center, all_centers, data, labels, tune = T) {
  m_k <- get_r_k(data[labels == k,, drop = F], gen_center) / 
    norm((gen_center - k_center), type = "2")
  r = m_k
  if (tune) {
    r <- optimize(get_proj_dist_for_centers, lower = 0, upper = m_k, 
                  k = k, k_center = k_center, gen_center = gen_center, 
                  all_centers = all_centers, data = data)$minimum
  }
  cat(r, "\n")
  all_centers[k, ] <- adj_center(r, k_center, gen_center)
  all_centers
}

weighted_k_means <- function(X, mu, K, num_obs) {
  error = c(Inf)
  i=1
  repeat{
    labels <- apply(X, 1, function(x) {which.min(sapply(mu, function(mu_k) {norm(x - mu_k, type = "2")^2}))})
    mu = lapply(1:K, function(k) {
      colSums(X[labels == k,, drop = F] * num_obs[labels == k]) / sum(num_obs[labels == k])
    })
    error = c(error, sum(sapply(1:nrow(X), function(l) {
      norm(X[l,, drop = F] - mu[[labels[l]]], type = "2")^2 * num_obs[l]}))) 
    if ((error[i] - error[i + 1]) < 1e-6) {
        break 
    }
    i=i+1 
  }
  return(list("labels" = labels, "mu" = mu, "error" = error[-1])) 
}

run_tuned_gdm <- function(data, inc_cat, num_centers, weighted = F, tune = T, C = colMeans(data[, inc_cat])) {
  data$k <- NA
  all_centers <- matrix()
  if (weighted) {
    rand_mu0 <- lapply(sample(1:nrow(data), num_centers), function(i) {data[i,inc_cat]})
    kmeans_results <- weighted_k_means(data[, inc_cat], rand_mu0, num_centers, data[, "weights"])
    data$k <- kmeans_results$labels
    all_centers <- do.call(rbind, kmeans_results$mu)
    kmeans_results$centers <- all_centers
  } else {
    kmeans_results <- kmeans(data[, inc_cat], centers = num_centers)
    data$k <- kmeans_results$cluster
    all_centers <- kmeans_results$centers
  }
  for (k in 1:num_centers) {
    all_centers <- tune_center_for_k(k, kmeans_results$centers[k,], C, all_centers, 
                                     data[, inc_cat], data$k, tune)
  }
  return(t(all_centers))
}


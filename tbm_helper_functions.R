#Tune VLAD
tune_vlad <- function(data, K) {
  data_center = colMeans(data)
  centered_data <- scale(data, scale = F)
  svd_info <- svd(centered_data)
  u_list <- svd_info$u[, 1:(K - 1)]
  kmeans_results <- kmeans(u_list, K, nstart = 25)
  v_list <- kmeans_results$centers
  topics_m <- matrix(0, nrow = K, ncol = ncol(data))
  all_centers <- t(svd_info$v[, 1:(K - 1)] %*% 
    (diag(svd_info$d[1:(K - 1)]) %*% t(v_list)))
  all_centers <- sweep(all_centers, 2, data_center, "+")
  for (k in 1:K) {
    all_centers <- tune_center_for_k(
      k, all_centers[k, ], data_center, all_centers, data, kmeans_results$cluster)
  }
  topics_m <- all_centers
  return(topics_m)
}


#Topic path
generate_topic_permutations <- function(num_topics) {
  perms <- do.call(expand.grid, lapply(2:num_topics, function(i) {2:num_topics}))
  perms <- perms[!apply(perms, 1, function(row) {any(duplicated(row))}),]
  cbind(1, perms)
}

create_possible_topic_path_m <- function(J, K) {
  possible_paths_m <- create_possible_path_m(J, K)
  topic_path_m <- create_topic_path_m(J, K)
  translate_path_to_topic_m(possible_paths_m, topic_path_m, J)
}

translate_path_to_topic_m <- function(path_m, topic_path_m, J) {
  for (i in 1:nrow(path_m)) {
    path_var_ind <- c(1, rep(0, ncol(path_m) - 1))
    for (j in 2:J) {
      path_var_ind[j] = path_m[i, j]
      #Checks which rows equals each other
      topic_ind <- which(colSums(t(topic_path_m[, 1:J]) == path_var_ind) == J)
      path_m[i, j] = topic_path_m[topic_ind, J + 1]
    }
  }
  return(path_m)
}

create_possible_path_m <- function(J, K) {
  possible_paths_list <- list(c(1))
  for (j in 2:J) {
    tmp <- possible_paths_list
    new_possible_paths_list <- list()
    for (path in tmp) {
      new_possible_paths_list <- c(new_possible_paths_list,
                                   lapply(1:K, function(k) {c(path, k)}))
    }
    possible_paths_list <- new_possible_paths_list
  }
  do.call(rbind, possible_paths_list)
}

create_topic_path_m <- function(J, K) {
  num_topics = sum(K^(0:(J - 1)))
  topic_doc_path_m <- matrix(0, nrow = num_topics, ncol = J + 1)
  index_list <- list(c(1))
  for (i in 1:num_topics) {
    item <- index_list[[i]]
    index_list <- c(index_list, lapply(1:K, function(k) {c(item, k)}))
    topic_doc_path_m[i, 1:length(item)] <- item
    topic_doc_path_m[i, J+1] <- i
  }
  return(topic_doc_path_m)
}


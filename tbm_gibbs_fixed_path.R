library(SALTSampler)
library(LaplacesDemon)

source("tbm_gibbs.R")

recover_unique_topics <- function(splitting_node, splitting_level, possible_paths_m) {
  all_paths <- which(possible_paths_m[, splitting_level] == splitting_node)
  unique(as.vector(possible_paths_m[all_paths,]))
}

count_path_choices <- function(sn_info_m, possible_paths_m) {
  path_choice_count <- rep(0, nrow(possible_paths_m))
  tmp_count <- table(sn_info_m[, 1])
  min_node <- min(possible_paths_m[, ncol(possible_paths_m)])
  path_choice_count[as.integer(names(tmp_count)) - min_node + 1] = tmp_count
  return(path_choice_count)
}

calc_doc_word_prob_subtree_node <- function(theta_d_for_doc, beta_prob_m, sn_info_m_for_doc, 
                                            possible_paths_m) {
  unique_topics <- recover_unique_topics(sn_info_m_for_doc[1], 
                                         sn_info_m_for_doc[2],
                                         possible_paths_m)
  doc_prob = colSums(sweep(beta_prob_m[unique_topics,], 1, theta_d_for_doc, "*"))
  return(doc_prob)
}

calc_log_posterior_subtree_node_same_depth_path <- 
    function(word_count_m, theta_d_m, beta_prob_m, 
             sn_info_m, node_prob,
             a_beta, a_theta, a_sn, possible_paths_m, 
             prior = T) {
      
  log_posterior_prob = 0
  if (prior) {
    log_posterior_prob =
      sum(ddirichlet(beta_prob_m, a_beta, log = T)) + 
      sum(ddirichlet(node_prob, a_sn, log = T))
    log_posterior_prob = log_posterior_prob +
      sum(apply(theta_d_m, 1, function(level_prob) {
        ddirichlet(theta_d_m, a_theta,
                   log = T)
      }))
  }
  
  # node_count <- rep(0, max(possible_paths_m))
  # nc_tmp <- table(sn_info_m[,1])
  # node_count[as.integer(names(nc_tmp))] <- nc_tmp
  log_posterior_prob = log_posterior_prob + 
    sum(count_path_choices(sn_info_m, possible_paths_m) * log(node_prob))
  
  for (d in 1:nrow(word_count_m)) {
    doc_prob <- calc_doc_word_prob_subtree_node(
      theta_d_m[d, 1:sn_info_m[d, 3]], beta_prob_m, sn_info_m[d,], 
      possible_paths_m)
    log_posterior_prob = log_posterior_prob + sum(word_count_m[d,] * log(doc_prob))
  }
  return(log_posterior_prob)
}

sample_theta_d_node_subtree <- function(theta_d_m, beta_prob_m,
                                        sn_info_m, node_prob,
                                        word_count_m,
                                        possible_paths_m,
                                        a_beta, a_theta_list, a_sn) {
  
  for (d in 1:nrow(theta_d_m)) {
    level_sample_order <- sample(sn_info_m[d, 3])
    prev_theta_d <- Logit(theta_d_m[d, 1:sn_info_m[d, 3]])
    prev_log_post_for_doc <- 
      calc_log_posterior_subtree_node_same_depth_path(
        word_count_m[d,, drop = F], 
        theta_d_m[d,, drop = F], beta_prob_m,
        sn_info_m[d,, drop = F], node_prob,
        a_beta, a_theta_list, a_sn, possible_paths_m)
    next_theta_d <- prev_theta_d
    next_log_post_for_doc <- prev_log_post_for_doc
    for (level in level_sample_order) {
      prop_info <- PropStep(prev_theta_d, level, rep(0.1, length(prev_theta_d)))
      next_theta_d = prop_info[1:length(prev_theta_d)]
      next_log_post_for_doc <- calc_log_posterior_subtree_node_same_depth_path(
        word_count_m[d,, drop = F], 
        matrix(invlogit(next_theta_d), nrow = 1), beta_prob_m,
        sn_info_m[d,, drop = F], node_prob,
        a_beta, a_theta_list, a_sn, possible_paths_m)
      log_accept_prob = attr(prop_info, "dbt") + 
        next_log_post_for_doc - prev_log_post_for_doc
      if (log(runif(1)) < log_accept_prob) {
        prev_log_post_for_doc = next_log_post_for_doc
        prev_theta_d <- next_theta_d
      }
    }
    theta_d_m[d, 1:sn_info_m[d, 3]] <- invlogit(prev_theta_d)
  }
  return(theta_d_m)
}

sample_theta_d_with_node_node_subtree <- function(
  theta_d_m, beta_prob_m, sn_info_m, node_prob,
  word_count_m, possible_paths_m,
  a_beta, a_theta_list, a_sn) {
  
  for (d in 1:nrow(theta_d_m)) {
    new_sn <- sample(max(possible_paths_m), 1, prob = node_prob)
    new_level <- which(apply(possible_paths_m, 2, function(col) {
      any(col == new_sn)
    }))
    new_num_topics <- length(recover_unique_topics(new_sn, new_level, possible_paths_m))
    new_sn_info <- matrix(c(new_sn, new_level, new_num_topics), nrow = 1)
    new_theta_d <- matrix(rdirichlet(1, a_theta_list[[new_sn]]), nrow = 1)
    prev_log_post_for_doc <- 
      calc_log_posterior_subtree_node(
        word_count_m[d,, drop = F], 
        theta_d_m[d,, drop = F], beta_prob_m,
        sn_info_m[d,, drop = F], node_prob,
        a_beta, a_theta_list, a_sn, possible_paths_m, prior = F)
    next_log_post_for_doc <- 
      calc_log_posterior_subtree_node(
        word_count_m[d,, drop = F], 
        new_theta_d, beta_prob_m,
        new_sn_info, node_prob,
        a_beta, a_theta_list, a_sn, possible_paths_m, prior = F)
    if (log(runif(1)) < (next_log_post_for_doc - prev_log_post_for_doc)) {
      sn_info_m[d,] <- new_sn_info
      theta_d_m[d, 1:sn_info_m[d, 3]] <- new_theta_d
    }
  }
  return(list(theta_d_m, sn_info_m))
}

sample_beta_prob_m_node_subtree <- function(theta_d_m, beta_prob_m, 
                                             sn_info_m, node_prob,
                                             word_count_m,
                                             possible_paths_m,
                                             a_beta, a_theta_list, a_sn) {
  
  topic_sample_order <- sample(nrow(beta_prob_m))
  for (k in topic_sample_order) {
    affected_docs <- which(
      apply(sn_info_m, 1, function(row) {
        any(k %in% recover_unique_topics(row[1], row[2], possible_paths_m))
      }))
    if (length(affected_docs) == 0) { 
      beta_prob_m[k,] <- rdirichlet(1, a_beta)
      next
    } 

    prev_beta_k <- Logit(beta_prob_m[k,])
    # prev_log_post <- calc_log_posterior_subtree_node(
    #   word_count_m[affected_docs,, drop = F], 
    #   theta_d_m[affected_docs,, drop = F], beta_prob_m,
    #   sn_info_m[affected_docs,, drop = F], node_prob, 
    #   a_beta, a_theta_list, a_sn, possible_paths_m)
    prev_log_post <- calc_log_posterior_subtree_node_same_depth_path(
      word_count_m[affected_docs,, drop = F],
      theta_d_m[affected_docs,, drop = F], beta_prob_m,
      sn_info_m[affected_docs,, drop = F], node_prob,
      a_beta, a_theta_list, a_sn, possible_paths_m)
    next_beta_k <- prev_beta_k
    next_log_post <- prev_log_post
    
    word_sample_order <- sample(ncol(beta_prob_m))
    check_nums <- floor(1:10 * ncol(beta_prob_m) / 10)
    #word_sample_order <- 
    #  sample(ncol(beta_prob_m), size = floor(0.1 * ncol(beta_prob_m)))
    for (i in 1:length(word_sample_order)) {
    #for (word in word_sample_order) {
      word = word_sample_order[i]
      #prop_info <- PropStep(prev_beta_k, word, rep(0.1, length(prev_beta_k)))
      prop_info <- PropStep(next_beta_k, word, rep(0.1, length(prev_beta_k)))
      next_beta_k = prop_info[1:length(prev_beta_k)]
      if (!(i %in% check_nums)) {
        next
      }
      beta_prob_m_tmp <- beta_prob_m
      beta_prob_m_tmp[k,] <- invlogit(next_beta_k)
      next_log_post <- calc_log_posterior_subtree_node_same_depth_path(
        word_count_m[affected_docs,, drop = F], 
        theta_d_m[affected_docs,, drop = F], beta_prob_m_tmp,
        sn_info_m[affected_docs,, drop = F], 
        node_prob, a_beta, a_theta_list, a_sn, 
        possible_paths_m)
      log_accept_prob = attr(prop_info, "dbt") + 
        next_log_post - prev_log_post
      if (log(runif(1)) < log_accept_prob) {
        prev_log_post = next_log_post
        prev_beta_k <- next_beta_k
      } else {
        next_log_post = prev_log_post
        next_beta_k = prev_beta_k
      }
    }
    #beta_prob_m[k,] <- invlogit(prev_beta_k)
    beta_prob_m[k,] <- invlogit(next_beta_k)
  }
  return(beta_prob_m)
}

sample_c_d_same_depth_path <- function(theta_d_m, beta_prob_m, sn_info_m, node_prob,
                                       word_count_m, possible_paths_m, 
                                       a_beta, a_theta_list, a_sn) {
  
  J = ncol(possible_paths_m)
  for (d in 1:nrow(theta_d_m)) {
    log_prob <- sapply(possible_paths_m[, J], function(leaf_node_ind) {
      new_sn_info <- matrix(c(leaf_node_ind, J, J), nrow = 1)
      calc_log_posterior_subtree_node_same_depth_path(
        word_count_m[d,, drop = F],
        theta_d_m[d,, drop = F], beta_prob_m,
        new_sn_info, node_prob,
        a_beta, a_theta_list, a_sn, possible_paths_m, prior = F)
    })
    log_prob <- log_prob - max(log_prob)
    sn_info_m[d,] <- c(sample(possible_paths_m[, J], 1, prob = exp(log_prob)), J, J)
  }
  return(sn_info_m)
}

#theta_d_m: matrix of D x num_topics with each row a probability for each topic
#beta_prob_m: matrix of K x num_words cols with each row a topic
#sn_info_m: matrix of D x 3 (node, level, number of topics)
#word_count_m <- matrix D x num_words with each row counting the number of words
initialize_data_doc_prob_same_depth_path <- function(
  doc_word_list, num_words, possible_paths_m, a_beta, a_theta_list, a_sn) {
  
  #Create possible paths
  # possible_paths_m <- create_possible_topic_path_m(J, K)
  
  # results$theta_d_m, results$beta_prob_m, 
  # results$sn_info_m, results$pi_prob_m,
  # results$node_prob
  
  J = ncol(possible_paths_m)
  possible_paths_m <- possible_paths_m[order(possible_paths_m[, J]),]
  num_paths = nrow(possible_paths_m)
  num_topics = max(possible_paths_m)
  beta_prob_m <- matrix(1 / num_words, nrow = num_topics, ncol = num_words)
  node_prob <- rep(1 / num_paths, num_paths)
  
  sn_info_m <- matrix(J, nrow = length(doc_word_list), ncol = 3)
  sn_info_m[, 1] <- sample(possible_paths_m[, J], length(doc_word_list), replace = T)
  theta_d_m <- matrix(1 / J, nrow = length(doc_word_list), ncol = J)
  word_count_m <- matrix(0, nrow = length(doc_word_list), ncol = num_words)
  
  for (d in 1:length(doc_word_list)) {
    word_count_for_doc <- table(doc_word_list[[d]])
    word_count_m[d, as.integer(names(word_count_for_doc))] <- word_count_for_doc
  }
  return(list("theta_d_m" = theta_d_m,
              "beta_prob_m" = beta_prob_m,
              "node_prob" = node_prob,
              "sn_info_m" = sn_info_m,
              "word_count_m" = word_count_m,
              "possible_paths_m" = possible_paths_m))
}

sample_lda_model_same_depth_path <- 
  function(results, possible_paths_m, word_count_m,
           num_iter, a_beta, a_theta_list, a_sn,
           pooled = F, keep_all = T) {
    
    # theta_d_m, beta_prob_m, s_d_v, 
    # sn_info_m, pi_prob_m, node_prob,
    # word_count_m,
    # possible_paths_m,
    # a_beta, a_theta, a_sn, a_pi_list
    
    sampling_results <- list()
    for (i in 1:num_iter) {
      print(i)
 
      results$beta_prob_m <- sample_beta_prob_m_node_subtree(
        results$theta_d_m, results$beta_prob_m, 
        results$sn_info_m, results$node_prob, 
        word_count_m, possible_paths_m, 
        a_beta, a_theta_list, a_sn)
      print("Done beta_prob_m")
    
      results$theta_d_m <- sample_theta_d_node_subtree(
        results$theta_d_m, results$beta_prob_m, 
        results$sn_info_m, results$node_prob, 
        word_count_m, possible_paths_m, 
        a_beta, a_theta_list, a_sn)
      print("Done theta_d")
      
      results$sn_info_m <- sample_c_d_same_depth_path(
        results$theta_d_m, results$beta_prob_m, 
        results$sn_info_m, results$node_prob, 
        word_count_m, possible_paths_m, 
        a_beta, a_theta_list, a_sn)
      print("Done c_d")
      
      results$node_prob <- 
        rdirichlet(1, a_sn + count_path_choices(results$sn_info_m, possible_paths_m))
    
      if (keep_all || (i >= num_iter / 2)) {
        sampling_results <- c(sampling_results, list(results))
      }
    }
    return(sampling_results)
}

sample_hier_lda_same_level_depth <- function(doc_word_list, num_words, possible_paths_m, 
                                             a_beta, a_theta_list, a_sn, num_iter = 1000,
                                             num_chains = 1, keep_all = T) {
  
  results <- initialize_data_doc_prob_same_depth_path(
    doc_word_list, num_words, possible_paths_m, a_beta, a_theta_list, a_sn)
  word_count_m <- results$word_count_m
  results <- results[!(names(results) %in%
                         c("word_count_m", "possible_paths_m"))]
  if (num_chains == 1) {
    return(sample_lda_model_same_depth_path(
      results, possible_paths_m, word_count_m, num_iter, 
      a_beta, a_theta_list, a_sn, keep_all = keep_all))
  } else {
    return(mclapply(1:num_chains, function(i) {
      sample_lda_model_same_depth_path(
        results, possible_paths_m, word_count_m, num_iter, 
        a_beta, a_theta_list, a_sn, keep_all = keep_all)
    }, mc.cores = num_chains))
  }
}


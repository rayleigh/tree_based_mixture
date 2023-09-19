library(parallel)
library(MCMCpack)
library(tidyverse)

source("tbm_helper_functions.R")

build_alpha_expectation_m <- function(alpha_v_m) {
  alpha_expectation_m <- 
    matrix(0, nrow = nrow(alpha_v_m), ncol = ncol(alpha_v_m))
  for (i in 1:nrow(alpha_expectation_m)) {
    for (j in 1:ncol(alpha_expectation_m)) {
      alpha_expectation_m[i, j] = log_dir_expectation(j, alpha_v_m[i,])
    }
  }
  return(alpha_expectation_m)
}

create_topics_prob_for_doc <- function(
  c_d_probs_m_for_doc, doc_paths_m, J, K) {
  
  topic_prob_for_doc <- matrix(0, nrow = J, ncol = max(doc_paths_m))
  topic_prob_for_doc[1, 1] = 1
  topic_start <- 1
  topic_end <- 1
  for (j in 2:J) {
    level_topics <- doc_paths_m[, j]
    topic_prob_for_doc[j, sort(unique(level_topics))]  <-
      data.frame("prob" = c_d_probs_m_for_doc, 
                 "label" = level_topics) %>%
      group_by(label) %>% summarize(topic_prob = sum(prob)) %>%
      dplyr::select(topic_prob) %>% unlist()
  }
  return(topic_prob_for_doc)
}

log_dir_expectation <- function(i, alpha) {
  digamma(alpha[i]) - digamma(sum(alpha))
}

log_dir_constant <- function(alpha) {
  sum(lgamma(alpha)) - lgamma(sum(alpha))
}

vb_update_for_theta_d <- function(z_d_n_probs_for_doc, alpha_theta) {
  alpha_theta + rowSums(z_d_n_probs_for_doc)
}

vb_update_for_theta_d_summarized <- function(
    summarized_z_d_n_probs_for_doc, word_count_for_doc, alpha_theta) {
  alpha_theta + summarized_z_d_n_probs_for_doc %*% word_count_for_doc
}

vb_update_for_beta_k <- function(k, l_k, word_m, z_d_n_probs_m, c_d_probs_m,
                                 alpha_beta, doc_paths_m) {
  beta_k_update <- alpha_beta
  relevant_paths <- which(doc_paths_m[, l_k] == k)
  for (i in 1:length(alpha_beta)) {
    for (d in 1:nrow(word_m)) {
      relevant_words <- which(word_m[d,] == i)
      beta_k_update[i] = beta_k_update[i] + sum(z_d_n_probs_m[l_k, relevant_words, d] *
        sum(c_d_probs_m[d, relevant_paths]))
    }
  }
  return(beta_k_update)
}

vb_update_for_beta_k_summarized_all <- 
  function(vb_alpha_beta_m, word_count_m, summarized_z_d_n_probs_m, c_d_probs_m,
           alpha_beta, doc_paths_m) {
  
  next_vb_alpha_beta_m <- 
    matrix(alpha_beta, nrow = nrow(vb_alpha_beta_m),
           ncol = ncol(vb_alpha_beta_m), byrow = T)
  doc_topic_path_info <- 
    as.data.frame(doc_paths_m) %>% tidyr::gather() %>% 
    transmute(level = gsub("V", "", key), label = value)
  #Get level of topic
  topic_level <-
    doc_topic_path_info %>% unique() %>% 
    dplyr::select(level) %>% unlist() %>% as.numeric()
  for (d in 1:nrow(word_count_m)) {
    #Get path prob for topic
    doc_topic_path_prob <- 
      doc_topic_path_info %>% 
      mutate(prob = rep(c_d_probs_m[d,], ncol(doc_paths_m))) %>%
      group_by(label) %>%
      summarize(topic_prob = sum(prob)) %>%
      dplyr::select(topic_prob) %>% unlist()
    next_vb_alpha_beta_m <- next_vb_alpha_beta_m +
      (doc_topic_path_prob %o% word_count_m[d,]) * 
        summarized_z_d_n_probs_m[topic_level,, d]
  }
  return(unname(next_vb_alpha_beta_m))
}

vb_update_for_beta_k_summarized_doc_parallel <- 
  function(word_count_m_for_doc, c_d_probs_m_for_doc, 
           summarized_z_d_n_probs_m_for_doc, 
           doc_topic_path_info, topic_level, J) {
    
    doc_topic_path_prob <- doc_topic_path_info %>% 
      mutate(prob = rep(c_d_probs_m_for_doc, J)) %>%
      group_by(label) %>%
      summarize(topic_prob = sum(prob)) %>%
      dplyr::select(topic_prob) %>% unlist()
    return((doc_topic_path_prob %o% word_count_m_for_doc) * 
              summarized_z_d_n_probs_m_for_doc[topic_level,])
}

vb_update_for_beta_k_summarized_all_parallel <- 
  function(vb_alpha_beta_m, word_count_m, summarized_z_d_n_probs_m, c_d_probs_m,
           alpha_beta, doc_paths_m) {
    
  doc_topic_path_info <- 
    as.data.frame(doc_paths_m) %>% tidyr::gather() %>% 
    transmute(level = gsub("V", "", key), label = value)
  #Get level of topic
  topic_level <-
    doc_topic_path_info %>% unique() %>% 
    dplyr::select(level) %>% unlist() %>% as.numeric()
  vb_alpha_beta_m_list <- mclapply(1:nrow(word_count_m), function(d) {
    #Get path prob for topic
    vb_update_for_beta_k_summarized_doc_parallel(
      word_count_m[d,], c_d_probs_m[d,], 
      summarized_z_d_n_probs_m[,,d], doc_topic_path_info,
      topic_level, ncol(doc_paths_m))
  }, mc.cores = 2)
  next_vb_alpha_beta_m <- 
    sweep(rowSums(array(unlist(vb_alpha_beta_m_list), 
                  dim = c(nrow = nrow(vb_alpha_beta_m),
                          ncol = ncol(vb_alpha_beta_m),
                          nrow(word_count_m))),
            dims = 2), 2, alpha_beta, "+")
  return(unname(next_vb_alpha_beta_m))
}

vb_update_for_pi <- function(c_d_probs_m, alpha_pi) {
  alpha_pi + colSums(c_d_probs_m)
}

vb_update_for_z_d_n <- function(word, c_d_probs_m_for_doc,
                                vb_alpha_theta_for_doc, vb_alpha_beta_m, 
                                doc_paths_m, J, K) {
  log_z_update <- rep(0, J)
  topic_start <- 1
  topic_end <- 1
  for (j in 1:J) {
    log_z_update[j] = log_dir_expectation(j, vb_alpha_theta_for_doc)
    for (k in topic_start:topic_end) {
      relevant_paths <- which(doc_paths_m[, j] == k)
      log_z_update[j] = log_z_update[j] + 
        log_dir_expectation(word, vb_alpha_beta_m[k,]) * 
        sum(c_d_probs_m_for_doc[relevant_paths])
    }
    topic_start = topic_end + 1
    topic_end = K^j + topic_end
  }
  return(exp(log_z_update) / sum(exp(log_z_update)))
}

vb_update_for_z_d_n_all <- function(c_d_probs_m, 
                                    vb_alpha_theta_m, vb_alpha_beta_m, 
                                    num_words, doc_paths_m, J, K) {
  
  vb_alpha_beta_exp_m <- build_alpha_expectation_m(vb_alpha_beta_m)
  summarized_z_d_n_probs_m <- 
    array(0, c(J, num_words, nrow(c_d_probs_m)))
  for (d in 1:nrow(c_d_probs_m)) {
    vb_alpha_theta_exp_m_for_doc <- 
      build_alpha_expectation_m(vb_alpha_theta_m[d,, drop = F])
    doc_topic_prob <- 
      create_topics_prob_for_doc(c_d_probs_m[d,], doc_paths_m, J, K)
    level_word_prob_m <- doc_topic_prob %*% vb_alpha_beta_exp_m
    level_word_prob_m <- apply(level_word_prob_m, 2, function(col) {
      col <- col + vb_alpha_theta_exp_m_for_doc
      col <- col - max(col)
      readjust_prob(exp(col) / sum(exp(col)))
    })
    summarized_z_d_n_probs_m[,, d] <- level_word_prob_m
  }
  return(summarized_z_d_n_probs_m)
}

vb_update_for_z_d_n_all_parallel <- function(
  c_d_probs_m,vb_alpha_theta_m, vb_alpha_beta_m, 
  num_words, doc_paths_m, J, K) {
  
  vb_alpha_beta_exp_m <- build_alpha_expectation_m(vb_alpha_beta_m)
  summarized_z_d_n_probs_m_list <- mclapply(1:nrow(c_d_probs_m), function(d) {
    vb_alpha_theta_exp_m_for_doc <- 
      build_alpha_expectation_m(vb_alpha_theta_m[d,, drop = F])
    doc_topic_prob <- 
      create_topics_prob_for_doc(c_d_probs_m[d,], doc_paths_m, J, K)
    level_word_prob_m <- doc_topic_prob %*% vb_alpha_beta_exp_m
    level_word_prob_m <- apply(level_word_prob_m, 2, function(col) {
      col <- col + vb_alpha_theta_exp_m_for_doc
      col <- col - max(col)
      readjust_prob(exp(col) / sum(exp(col)))
    })
    level_word_prob_m
  }, mc.cores = 2)
  summarized_z_d_n_probs_m <- 
    array(unlist(summarized_z_d_n_probs_m_list), 
          dim = c(J, num_words, nrow(c_d_probs_m)))
  return(summarized_z_d_n_probs_m)
}

vb_update_for_c_d <- function(word_m_for_doc, z_d_n_probs_m_for_doc,
                              vb_alpha_beta_m, vb_alpha_pi, num_words, doc_paths_m, J, K) {
  log_topic_prob <- rep(0, sum(sapply(1:J - 1, function(j) {K^j})))
  topic_start <- 1
  topic_end <- 1
  for (j in 1:J) {
    for (k in topic_start:topic_end) {
      for (i in 1:num_words) {
        relevant_words <- which(word_m_for_doc == i)
        log_topic_prob[k] = log_topic_prob[k] + sum(z_d_n_probs_m_for_doc[j, relevant_words] * 
          log_dir_expectation(i, vb_alpha_beta_m[k,]))
      }
    }
    topic_start = topic_end + 1
    topic_end = K^j + topic_end
  }
  log_c_d_prob <- sapply(1:nrow(doc_paths_m), function(i) {
    path = doc_paths_m[i,]
    sum(log_topic_prob[path]) + log_dir_expectation(i, vb_alpha_pi)
  })
  log_c_d_prob = log_c_d_prob - max(log_c_d_prob)
  exp(log_c_d_prob) / sum(exp(log_c_d_prob))
}

vb_update_for_c_d_summarized <- 
  function(word_count_m_for_doc, 
           summarized_z_d_n_probs_for_doc,
           vb_alpha_beta_m, vb_alpha_pi, num_words, 
           doc_paths_m, J, K) {
  
  log_topic_prob <- rep(0, nrow(doc_paths_m))  
  vb_alpha_beta_exp_m <- build_alpha_expectation_m(vb_alpha_beta_m)  
  for (i in 1:length(log_topic_prob)) {
    log_topic_prob[i] =
      sum((summarized_z_d_n_probs_for_doc * vb_alpha_beta_exp_m[doc_paths_m[i,],]) %*%
      word_count_m_for_doc) + log_dir_expectation(i, vb_alpha_pi)
  }
  log_topic_prob <- log_topic_prob - max(log_topic_prob)
  return(exp(log_topic_prob) / sum(exp(log_topic_prob)))
}

vb_update_for_c_d_summarized_all <- 
  function(c_d_prob_m, word_count_m, 
           summarized_z_d_n_probs,
           vb_alpha_beta_m, vb_alpha_pi, num_words, 
           doc_paths_m, J, K) {
    
    vb_alpha_beta_exp_m <- build_alpha_expectation_m(vb_alpha_beta_m)
    vb_alpha_pi_exp_m <- build_alpha_expectation_m(matrix(vb_alpha_pi, nrow = 1))
    for (i in 1:ncol(c_d_prob_m)) {
      (summarized_z_d_n_probs * 
         array(vb_alpha_beta_exp_m[doc_paths_m[i,],], 
               c(J, num_words, nrow(word_count_m)))) %*% 
        array(t(word_count_m), c(num_words, 1, nrow(word_count_m)))
    }  
    
  }

vb_elbo_hier <- function(word_m, z_d_n_probs_m, c_d_probs_m,
                         vb_alpha_beta_m, vb_alpha_theta_m, vb_alpha_pi,
                         alpha_theta, alpha_beta, alpha_pi,
                         num_words, doc_paths_m, J, K) {
  log_vb_alpha_beta_e_m <- 
    t(apply(vb_alpha_beta_m, 1, function(row) {
      sapply(1:length(row), function(i) {
        log_dir_expectation(i, row)})  
    }))
  log_vb_alpha_theta_e_m <- 
    t(apply(vb_alpha_theta_m, 1, function(row) {
      sapply(1:length(row), function(i) {
        log_dir_expectation(i, row)})  
    }))
  log_vb_alpha_pi_e_m <- 
    sapply(1:length(vb_alpha_pi), function(i) 
      log_dir_expectation(i, vb_alpha_pi))
  topic_start <- 1
  topic_end <- 1
  #Expectation of log(P(pi | alpha_pi))
  elbo = sum((alpha_pi - 1) * log_vb_alpha_pi_e_m) -
            nrow(doc_paths_m) * log_dir_constant(alpha_pi)
  #Expectation of log(q(pi)) 
  elbo = elbo - (sum((vb_alpha_pi - 1) * log_vb_alpha_pi_e_m) -
                   log_dir_constant(vb_alpha_pi))
  for (j in 1:J) {
    for (k in topic_start:topic_end) {
      relevant_paths <- which(doc_paths_m[, j] == k)
      for (d in 1:nrow(word_m)) {
        for (i in 1:num_words) {
          relevant_words <- which(word_m[d,] == i)
          #Expectation of log(P(word | beta, theta, c_d))
          elbo = elbo + sum(z_d_n_probs_m[j, relevant_words, d] *
            sum(c_d_probs_m[d, relevant_paths]) * log_vb_alpha_beta_e_m[k, i])
        }
      }
      #Expectation of log(P(beta | alpha_beta))
      elbo = elbo + sum((alpha_beta - 1) * log_vb_alpha_beta_e_m[k,]) -
        log_dir_constant(alpha_beta)
      #Expectation of log(q(beta))
      elbo = elbo - (sum((vb_alpha_beta_m[k,] - 1) * log_vb_alpha_beta_e_m[k,]) -
                       log_dir_constant(vb_alpha_beta_m[k,]))
    }
    topic_start = topic_end + 1
    topic_end = K^j + topic_end
  }
  for (d in 1:nrow(word_m)) {
    #Expectation of log(P(z | theta))
    elbo = elbo + sum(rowSums(z_d_n_probs_m[,,d]) * 
                        log_vb_alpha_theta_e_m[d,])
    #Expectation of log(P(theta | alpha_theta))
    elbo = elbo + sum((alpha_theta - 1) * (log_vb_alpha_theta_e_m[d,])) -
      log_dir_constant(alpha_theta)
    #Expectation of log(P(c_d | pi))
    elbo = elbo + 
      sum(c_d_probs_m[d,] * log_vb_alpha_pi_e_m)
    #Expectation of log(q(z_d, n))
    relevant_inds <- which(word_m[d,] != -1)
    elbo = elbo - sum(z_d_n_probs_m[,relevant_inds,d] * 
                        log(z_d_n_probs_m[,relevant_inds,d]))
    #Expectation of log(q(c_d))
    elbo = elbo - sum(c_d_probs_m[d,] * log(c_d_probs_m[d,]))
    #Expectation of log(q(theta))
    elbo = elbo - (sum((vb_alpha_theta_m[d,] - 1) * log_vb_alpha_theta_e_m[d,]) -
                     log_dir_constant(vb_alpha_theta_m[d,]))
  }
  return(elbo)
}

vb_elbo_hier_summarized <- 
  function(word_count_m, summarized_z_d_n_probs_m,
           c_d_probs_m, vb_alpha_beta_m, vb_alpha_theta_m, vb_alpha_pi,
           alpha_theta, alpha_beta, alpha_pi) {
    
  log_vb_alpha_beta_e_m <- build_alpha_expectation_m(vb_alpha_beta_m)
  log_vb_alpha_theta_e_m <- build_alpha_expectation_m(vb_alpha_theta_m)
  log_vb_alpha_pi_e_m <- build_alpha_expectation_m(matrix(vb_alpha_pi, nrow = 1))
  
  #Expectation of log(P(word | beta, theta, c_d))
  elbo = sum(sweep(vb_alpha_beta_m, 2, alpha_beta, "-") * 
               log_vb_alpha_beta_e_m)
  
  #Expectation of log(P(z | theta))
  elbo = elbo + 
    sum(sweep(vb_alpha_theta_m, 2, alpha_theta, "-") * 
          log_vb_alpha_theta_e_m)

  #Expectation of log(P(c_d | pi))
  elbo = elbo + sum(c_d_probs_m %*% t(log_vb_alpha_pi_e_m))
  
  #Expectation of log(q(z_d, n))
  for (d in 1:nrow(word_count_m)) {
    elbo = elbo - 
      sum((summarized_z_d_n_probs_m[,,d] * 
            log(summarized_z_d_n_probs_m[,,d])) %*% word_count_m[d,])
  }

  #Expectation of log(q(c_d))
  elbo = elbo - sum(c_d_probs_m * log(c_d_probs_m))
  
  #Remaining cancel out
  
  return(elbo)
}

#Creates a data list to pass into vb
#word_m: D x max word in doc matrix with each row representing the words in a document in order
#word_count_m: D x num_words in doc matrix with each row representing the word_count
#z_d_n_probs_m: num levels x max word in doc x D tensor representing the probability
#of a level for each word in a document in the same order as word_m
#c_d_probs_m: D x num possible paths matrix representing the probability of a path
#for a document
#vb_alpha_beta_m: K x num_words representing the updated Dirichlet parameter for each topic
#vb_alpha_theta_m: D x num levels representing the updated Dirichlet parameter for the levels 
#in a document
#vb_alpha_pi: 1 x num_paths representing the updated Dirichlet parameter for the path probability 
#in a document
#possible_paths_m: D x num paths representing the possible paths
#doc_word_length: D vector listing the number of words in a document
initialize_data_for_vb <- function(
  doc_word_list, num_words, J, K, a_beta, a_theta, a_pi,
  pooled = F, num_topics = sum(K^(0:(J - 1))), parallel = F) {
  #Create possible paths
  if (pooled) {
    possible_paths_m <- as.matrix(generate_topic_permutations(num_topics))
  } else {
    possible_paths_m <- create_possible_topic_path_m(J, K)
  }
 
  print("Init data")  
  start_time <- proc.time()
  word_count_m <- matrix(0, nrow = length(doc_word_list), ncol = num_words)
  summarized_z_d_n_prob_m <- 
    array(0, c(ncol(possible_paths_m), num_words, length(doc_word_list)))
  vb_alpha_theta_m <- matrix(1, nrow = length(doc_word_list), 
                             ncol = ncol(possible_paths_m))
  c_d_prob_m <- rdirichlet(length(doc_word_list), rep(1, nrow(possible_paths_m)))
  print(proc.time() - start_time)
    
  start_time <- proc.time()
  for (d in 1:length(doc_word_list)) {
    word_count <- table(doc_word_list[[d]])
    word_count_m[d, as.integer(names(word_count))] <- word_count
    summarized_z_d_n_prob_m[,,d] <- t(rdirichlet(num_words, rep(1, ncol(possible_paths_m))))
    vb_alpha_theta_m[d,] <- 
      vb_update_for_theta_d_summarized(
        summarized_z_d_n_prob_m[,,d], word_count_m[d,], a_theta)
  }
  print(proc.time() - start_time)

  start_time <- proc.time()
  vb_alpha_pi <- vb_update_for_pi(c_d_prob_m, a_pi)
  vb_alpha_beta_m <- matrix(1, nrow = num_topics, ncol = num_words)
  if (parallel) {
    vb_alpha_beta_m <- vb_update_for_beta_k_summarized_all_parallel(
      vb_alpha_beta_m, word_count_m, summarized_z_d_n_prob_m,
      c_d_prob_m, a_beta, possible_paths_m)

  } else {
    vb_alpha_beta_m <- vb_update_for_beta_k_summarized_all(
      vb_alpha_beta_m, word_count_m, summarized_z_d_n_prob_m,
      c_d_prob_m, a_beta, possible_paths_m)
  }  
  print(proc.time() - start_time)

  return(list("word_count_m" = word_count_m,
              "summarized_z_d_n_probs_m" = summarized_z_d_n_prob_m, 
              "c_d_probs_m" = c_d_prob_m, 
              "vb_alpha_beta_m" = vb_alpha_beta_m, 
              "vb_alpha_theta_m" = vb_alpha_theta_m,
              "vb_alpha_pi" = vb_alpha_pi,
              "possible_paths_m" = possible_paths_m))
}

readjust_prob <- function(prob_v) {

    zero_ind <- which(prob_v < 1e-9)
    if (length(zero_ind) == 0) {
      return(prob_v)
    }
    prob_v[zero_ind] = 1e-9
    return(prob_v / sum(prob_v))
    # prob_v[-zero_ind] = prob_v[-zero_ind] - 
    #   1e-9 * length(zero_ind) / (length(prob_v) - length(zero_ind))
    # prob_v  
}

vb_update_hier <- function(word_count_m, summarized_z_d_n_probs_m, c_d_probs_m, 
                           vb_alpha_beta_m, vb_alpha_theta_m, vb_alpha_pi,
                           alpha_theta, alpha_beta, alpha_pi,
                           num_words, doc_paths_m, J, K,
                           tol = 1e-3, iter_max = 100) {
  prev_elbo = -Inf
  i = 1
  repeat {
    #Update z_d_n
     summarized_z_d_n_probs_m <- vb_update_for_z_d_n_all(
      c_d_probs_m, vb_alpha_theta_m, vb_alpha_beta_m, 
      num_words, doc_paths_m, J, K)
    
    #Update theta
    for (d in 1:nrow(vb_alpha_theta_m)) {
      vb_alpha_theta_m[d,] <- 
        vb_update_for_theta_d_summarized(
          summarized_z_d_n_probs_m[,,d], word_count_m[d,], alpha_theta)
    }
    
    #Update c_d
    for (d in 1:nrow(c_d_probs_m)) {
      c_d_probs_m[d,] <- 
        readjust_prob(vb_update_for_c_d_summarized(word_count_m[d,], 
                                     summarized_z_d_n_probs_m[,,d],
                                     vb_alpha_beta_m, vb_alpha_pi, num_words, 
                                     doc_paths_m, J, K))
    }
    
    #Update pi
    vb_alpha_pi <- vb_update_for_pi(c_d_probs_m, alpha_pi)
    
    #Update beta
    vb_alpha_beta_m <- vb_update_for_beta_k_summarized_all(
      vb_alpha_beta_m, word_count_m, summarized_z_d_n_probs_m, c_d_probs_m,
      alpha_beta, doc_paths_m)
    
    #Check ELBO
    next_elbo = vb_elbo_hier_summarized(
      word_count_m, summarized_z_d_n_probs_m,
      c_d_probs_m, vb_alpha_beta_m, vb_alpha_theta_m, vb_alpha_pi,
      alpha_theta, alpha_beta, alpha_pi)
    print(sprintf("%.9f",next_elbo))
    if (abs(next_elbo - prev_elbo) < tol || i >= iter_max) {
      break
    }
    i = i + 1
    prev_elbo = next_elbo
  }
  return(list("z_d_n_probs_m" = summarized_z_d_n_probs_m,
              "c_d_probs_m" = c_d_probs_m,
              "vb_alpha_beta_m" = vb_alpha_beta_m,
              "vb_alpha_theta_m" = vb_alpha_theta_m,
              "vb_alpha_pi" = vb_alpha_pi))
}

vb_update_hier_parallel <- 
    function(word_count_m, summarized_z_d_n_probs_m, c_d_probs_m, 
             vb_alpha_beta_m, vb_alpha_theta_m, vb_alpha_pi,
             alpha_theta, alpha_beta, alpha_pi,
             num_words, doc_paths_m, J, K,
             tol = 1e-3, iter_max = 100) {

  prev_elbo = -Inf
  i = 1
  repeat {
    #Update z_d_n
    print("zeta")
    start_time <- proc.time()
    summarized_z_d_n_probs_m <- vb_update_for_z_d_n_all_parallel(
      c_d_probs_m, vb_alpha_theta_m, vb_alpha_beta_m, 
      num_words, doc_paths_m, J, K)
    print(proc.time() - start_time)
    
    #Update theta
    print("beta")
    start_time <- proc.time()
    vb_alpha_theta_m <- 
      do.call(rbind, mclapply(1:nrow(vb_alpha_theta_m), function(d) {
        t(vb_update_for_theta_d_summarized(
          summarized_z_d_n_probs_m[,,d], word_count_m[d,], alpha_theta))
      }, mc.cores = 9))
    print(proc.time() - start_time)
    
    #Update c_d
    print("z")
    start_time <- proc.time()
    c_d_probs_m <- do.call(rbind, 
       mclapply(1:nrow(c_d_probs_m), function(d) {
         readjust_prob(vb_update_for_c_d_summarized(
           word_count_m[d,], summarized_z_d_n_probs_m[,,d],
           vb_alpha_beta_m, vb_alpha_pi, num_words, 
           doc_paths_m, J, K))
       }, mc.cores = 9))
    print(proc.time() - start_time)

    #Update pi
    print("pi")
    start_time <- proc.time()
    vb_alpha_pi <- vb_update_for_pi(c_d_probs_m, alpha_pi)
    print(proc.time() - start_time)
    
    #Update beta
    print("theta")
    start_time <- proc.time()
    vb_alpha_beta_m <- vb_update_for_beta_k_summarized_all_parallel(
      vb_alpha_beta_m, word_count_m, summarized_z_d_n_probs_m, c_d_probs_m,
      alpha_beta, doc_paths_m)
    print(proc.time() - start_time)

    #Check ELBO
    print("elbo")
    start_time <- proc.time()
    next_elbo = vb_elbo_hier_summarized(
      word_count_m, summarized_z_d_n_probs_m, c_d_probs_m, 
      vb_alpha_beta_m, vb_alpha_theta_m, vb_alpha_pi,
      alpha_theta, alpha_beta, alpha_pi)
    print(proc.time() - start_time)
    print(sprintf("%.9f",next_elbo))
    if (abs(next_elbo - prev_elbo) < tol || i >= iter_max) {
      break
    }
    i = i + 1
    prev_elbo = next_elbo
  }
  return(list("z_d_n_probs_m" = summarized_z_d_n_probs_m,
              "c_d_probs_m" = c_d_probs_m,
              "vb_alpha_beta_m" = vb_alpha_beta_m,
              "vb_alpha_theta_m" = vb_alpha_theta_m,
              "vb_alpha_pi" = vb_alpha_pi))
}

run_vb_hier <- function(doc_word_list, num_words, J, K, 
                        a_beta, a_theta, a_pi, tol = 1e-3, iter_max = 100) {
  param_list <- initialize_data_for_vb(
    doc_word_list, num_words, J, K, a_beta, a_theta, a_pi, pooled = F)
  vb_update_hier(param_list$word_count_m, 
                 param_list$summarized_z_d_n_probs_m, 
                 param_list$c_d_probs_m, param_list$vb_alpha_beta_m, 
                 param_list$vb_alpha_theta_m, param_list$vb_alpha_pi,
                 a_theta, a_beta, a_pi,
                 num_words, param_list$possible_paths_m, 
                 J, K, tol = tol, iter_max = iter_max)
}

run_vb_hier_parallel <- 
  function(doc_word_list, num_words, J, K, 
           a_beta, a_theta, a_pi, tol = 1e-3, iter_max = 100) {
  param_list <- initialize_data_for_vb(
    doc_word_list, num_words, J, K, a_beta, a_theta, a_pi, pooled = F, 
    parallel = T)
  vb_update_hier_parallel(
    param_list$word_count_m, param_list$summarized_z_d_n_probs_m, 
    param_list$c_d_probs_m, param_list$vb_alpha_beta_m, 
    param_list$vb_alpha_theta_m, param_list$vb_alpha_pi,
    a_theta, a_beta, a_pi, num_words, param_list$possible_paths_m, 
    J, K, tol = tol, iter_max = iter_max)
}

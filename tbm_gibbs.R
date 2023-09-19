library(parallel)

source("tbm_helper_functions.R")

#topic_word_count_m: topic x word matrix that keeps counts of number of words 
#assigned to topic for word
#topic_word_assignment_m: matrix with word, topic, and level assignment
#
sample_z_for_doc <- function(topic_word_count_m_for_doc,
                             topic_word_assignment_for_doc,
                             topic_choices_for_doc, a_beta, a_theta) {
  num_z = length(topic_choices_for_doc)
  z_counts <- rep(0, num_z)
  z_word_counts <- table(topic_word_assignment_for_doc[, 3])
  z_counts[which(1:length(topic_choices_for_doc) %in% 
                   names(z_word_counts))] <- z_word_counts
  for (i in 1:nrow(topic_word_assignment_for_doc)) {
    #Remove current word info
    curr_z_ind = topic_word_assignment_for_doc[i, 3]
    z_counts[curr_z_ind] = z_counts[curr_z_ind] - 1
    topic_word_count_m_for_doc[curr_z_ind, topic_word_assignment_for_doc[i, 1]] =
      topic_word_count_m_for_doc[curr_z_ind, topic_word_assignment_for_doc[i, 1]] - 1
    
    #Sample new topic assignment
    topic_word_denom <- (rowSums(topic_word_count_m_for_doc) + sum(a_beta))
    topic_probs <- (z_counts + a_theta) * 
      (topic_word_count_m_for_doc[, topic_word_assignment_for_doc[i, 1]] + 
         a_beta[topic_word_assignment_for_doc[i, 1]]) / topic_word_denom
    new_topic_ind <- sample(length(topic_choices_for_doc), 1, prob = topic_probs)

    #Update info
    z_counts[new_topic_ind] = z_counts[new_topic_ind] + 1
    topic_word_count_m_for_doc[new_topic_ind, topic_word_assignment_for_doc[i, 1]] =
      topic_word_count_m_for_doc[new_topic_ind, topic_word_assignment_for_doc[i, 1]] + 1
    topic_word_assignment_for_doc[i, 2] <- topic_choices_for_doc[new_topic_ind]
    topic_word_assignment_for_doc[i, 3] <- new_topic_ind
  }
  return(list("topic_word_count" = topic_word_count_m_for_doc,
              "topic_word_assignment" = topic_word_assignment_for_doc))
}

sample_z_for_a_level <- function(level, level_word_count_m, 
                               next_level_word_count_m_list,
                               total_level_word_count_m,
                               topic_word_count_m, total_topic_word_count,
                               doc_topic_assignment_m,
                               a_beta, a_theta) {
  for (i in 1:nrow(level_word_count_m)) {
    curr_topic = doc_topic_assignment_m[i, level]
    for (j in 1:length(level_word_count_m[i,])) {
      if (level_word_count_m[i, j] == 0) {
        next
      }
      for (k in 1:level_word_count_m[i, j]) {
        #Remove current word info
        next_level_word_count_m_list[[level]][i, j] = 
          next_level_word_count_m_list[[level]][i, j] - 1
        total_level_word_count_m[i, level] = total_level_word_count_m[i, level] - 1
        topic_word_count_m[curr_topic, j] = topic_word_count_m[curr_topic, j] - 1
        total_topic_word_count[curr_topic] = total_topic_word_count[curr_topic] - 1
        
        #Sample new topic assignment
        topic_word_denom <- total_topic_word_count[doc_topic_assignment_m[i,]] + sum(a_beta)
        topic_probs <- (total_level_word_count_m[i,] + a_theta) * 
          (topic_word_count_m[doc_topic_assignment_m[i,], j] + 
             a_beta[j]) / topic_word_denom
        new_topic_ind <- sample(ncol(total_level_word_count_m), 1, prob = topic_probs)
        
        #Update info
        next_level_word_count_m_list[[new_topic_ind]][i, j] =
          next_level_word_count_m_list[[new_topic_ind]][i, j] + 1
        total_level_word_count_m[i, new_topic_ind] = 
          total_level_word_count_m[i, new_topic_ind] + 1
        new_topic = doc_topic_assignment_m[i, new_topic_ind]
        topic_word_count_m[new_topic, j] = topic_word_count_m[new_topic, j] + 1
        total_topic_word_count[new_topic] = total_topic_word_count[new_topic] + 1
      }
    }
  }
  return(list("level_word_count_m_list" = next_level_word_count_m_list,
              "total_level_word_count_m" = total_level_word_count_m,
              "topic_word_count_m" = topic_word_count_m,
              "total_topic_word_count" = total_topic_word_count))
}

sample_z_for_doc_level <- function(d, sample_order_counts, sample_info_m, 
                                   next_level_word_count_m_list,
                                   total_level_word_count_m,
                                   topic_word_count_m, total_topic_word_count,
                                   doc_topic_assignment_m,
                                   a_beta, a_theta) {
  sample_order_ind_v <- sample(rep(1:nrow(sample_info_m), 
                                   times = sample_order_counts))
  for (ind in sample_order_ind_v) {
    level = sample_info_m[ind, 1]
    j = sample_info_m[ind, 2]
    curr_topic = doc_topic_assignment_m[d, level]
    next_level_word_count_m_list[d, j, level] = 
      next_level_word_count_m_list[d, j, level] - 1
    total_level_word_count_m[d, level] = total_level_word_count_m[d, level] - 1
    topic_word_count_m[curr_topic, j] = topic_word_count_m[curr_topic, j] - 1
    total_topic_word_count[curr_topic] = total_topic_word_count[curr_topic] - 1
    
    #Sample new topic assignment
    topic_word_denom <- total_topic_word_count[doc_topic_assignment_m[d,]] + sum(a_beta)
    topic_probs <- (total_level_word_count_m[d,] + a_theta) * 
      (topic_word_count_m[doc_topic_assignment_m[d,], j] + 
         a_beta[j]) / topic_word_denom
    new_topic_ind <- sample(ncol(total_level_word_count_m), 1, prob = topic_probs)
    
    #Update info
    next_level_word_count_m_list[d, j, new_topic_ind] =
      next_level_word_count_m_list[d, j, new_topic_ind] + 1
    total_level_word_count_m[d, new_topic_ind] = 
      total_level_word_count_m[d, new_topic_ind] + 1
    new_topic = doc_topic_assignment_m[d, new_topic_ind]
    topic_word_count_m[new_topic, j] = topic_word_count_m[new_topic, j] + 1
    total_topic_word_count[new_topic] = total_topic_word_count[new_topic] + 1
  }
  return(list("level_word_count_m_list" = next_level_word_count_m_list,
              "total_level_word_count_m" = total_level_word_count_m,
              "topic_word_count_m" = topic_word_count_m,
              "total_topic_word_count" = total_topic_word_count))
}

sample_z <- function(topic_word_count_m, topic_word_assignment_list, 
                     doc_path_assignment_m, a_beta, a_theta) {
  for (d in 1:nrow(doc_path_assignment_m)) {
    doc_update <- sample_z_for_doc(topic_word_count_m[doc_path_assignment_m[d,],],
                                   topic_word_assignment_list[[d]], 
                                   doc_path_assignment_m[d,], a_beta, a_theta)
    topic_word_count_m[doc_path_assignment_m[d,],] <- doc_update$topic_word_count
    topic_word_assignment_list[[d]] <- doc_update$topic_word_assignment
  }
  return(list("topic_word_count" = topic_word_count_m,
              "topic_word_assignment" = topic_word_assignment_list))
}

sample_z_for_all_levels <- function(level_word_count_m_list,
                                    total_level_word_count_m,
                                    topic_word_count_m, total_topic_word_count,
                                    doc_topic_assignment_m,
                                    a_beta, a_theta) {
  next_level_word_count_m_list <- level_word_count_m_list
  level_assignment_m <- expand.grid(1:dim(level_word_count_m_list)[3],
                                    1:ncol(topic_word_count_m))
  for (d in 1:nrow(doc_topic_assignment_m)) {
    word_assignment_counts <- t(level_word_count_m_list[d,,])
    level_sample_info <- sample_z_for_doc_level(d,
      word_assignment_counts, level_assignment_m, next_level_word_count_m_list,
      total_level_word_count_m, topic_word_count_m, total_topic_word_count,
      doc_topic_assignment_m, a_beta, a_theta)
    next_level_word_count_m_list <- level_sample_info[[1]]
    total_level_word_count_m <- level_sample_info[[2]]
    topic_word_count_m <- level_sample_info[[3]]
    total_topic_word_count <- level_sample_info[[4]]
  }
  return(list("level_word_count_m_list" = next_level_word_count_m_list,
              "total_level_word_count_m" = total_level_word_count_m,
              "topic_word_count_m" = topic_word_count_m,
              "total_topic_word_count" = total_topic_word_count))
}


#possible_paths_m matrix of topics that represent possible paths
#Assumes topics are in sorted order
sample_c <- function(topic_word_count_m, topic_word_assignment_list, 
                     doc_path_assignment_m, a_beta, possible_paths_m,
                     pooled = F) {
  for (d in 1:nrow(doc_path_assignment_m)) {
    #Remove existing information
    cur_path = doc_path_assignment_m[d,]
    word_counts <- table(topic_word_assignment_list[[d]][, c(2, 1)])
    word_assignment_counts <- matrix(0, nrow = ncol(doc_path_assignment_m),
                                     ncol = ncol(word_counts))
    word_assignment_counts[which(sort(cur_path) 
                           %in% rownames(word_counts)),] <- word_counts
    doc_possible_paths_m <- possible_paths_m
    if (pooled) {
      word_assignment_counts <- word_assignment_counts[cur_path,]
      subset_ind <- sort(unique(topic_word_assignment_list[[d]][, 3]))
      doc_possible_paths_m <- doc_possible_paths_m[
        !duplicated(doc_possible_paths_m[,subset_ind]), , drop = F]
    }
    word_list <- as.integer(colnames(word_counts))
    topic_word_count_m[cur_path, word_list] <-
      topic_word_count_m[cur_path, word_list] - word_assignment_counts
    
    #word_assignment_level_counts <- table(topic_word_assignment_list[[d]][, c(3, 1)])
    words_per_level <- rowSums(word_assignment_counts)
    path_probs = rep(0, nrow(doc_possible_paths_m))
    for (i in 1:length(path_probs)) {
      path = doc_possible_paths_m[i,]
      for (j in 1:length(path)) {
        if (words_per_level[j] == 0) {
          next
        }
        words_ind_for_topic <- 1:length(word_list)
        zero_word_ind = which(word_assignment_counts[j,] == 0)
        if (length(zero_word_ind) > 0) {
          words_ind_for_topic <- words_ind_for_topic[-zero_word_ind]
        }
        prob_numer = sum(sapply(words_ind_for_topic, function(w_i) {
          sum(log(topic_word_count_m[path[j], word_list[w_i]] + a_beta[word_list[w_i]] +
                 0:(word_assignment_counts[j, w_i] - 1)))
        }))
        prob_denom = sum(log(sum(topic_word_count_m[path[j],]) + sum(a_beta) + 
                       0:(words_per_level[j] - 1)))
        path_probs[i] = path_probs[i] + prob_numer - prob_denom
      }
    }
    new_path <- doc_possible_paths_m[sample(length(path_probs), 1, 
                                            prob = exp(path_probs - max(path_probs))),]
    if (pooled && (length(subset_ind) < (length(cur_path) - 1))) {
      new_path[-subset_ind] <- sample(new_path[-subset_ind])
    }

    #Update existing information
    topic_word_count_m[new_path, word_list] <- 
      topic_word_count_m[new_path, word_list] + word_assignment_counts
    topic_word_assignment_list[[d]][, 2] <- new_path[topic_word_assignment_list[[d]][, 3]]
    doc_path_assignment_m[d,] <- new_path
  }
  return(list("topic_word_count" = topic_word_count_m,
              "topic_word_assignment" = topic_word_assignment_list,
              "doc_path_assignment" = doc_path_assignment_m))
}

sample_c_for_all_level <- function(level_word_count_m_list,
                                   total_level_word_count_m,
                                   topic_word_count_m, total_topic_word_count,
                                   doc_path_assignment_m,
                                   a_beta, possible_paths_m, pooled = F) {
  for (d in 1:nrow(doc_path_assignment_m)) {
    #Remove existing information
    cur_path = doc_path_assignment_m[d,]
    word_assignment_counts <- t(level_word_count_m_list[d,,])

    #word_assignment_level_counts <- table(topic_word_assignment_list[[d]][, c(3, 1)])
    doc_possible_paths_m <- possible_paths_m
    # if (pooled) {
    #   word_assignment_counts <- word_assignment_counts[cur_path,]
    #   subset_ind <- sort(unique(topic_word_assignment_list[[d]][, 3]))
    #   doc_possible_paths_m <- doc_possible_paths_m[
    #     !duplicated(doc_possible_paths_m[,subset_ind]),, drop = F]
    # }
    topic_word_count_m[cur_path,] <-
      topic_word_count_m[cur_path,] - word_assignment_counts
    total_topic_word_count[cur_path] <- total_topic_word_count[cur_path] - 
      total_level_word_count_m[d,]
    words_per_level <- total_level_word_count_m[d,]
    path_probs = rep(0, nrow(doc_possible_paths_m))
    for (i in 1:length(path_probs)) {
      path = doc_possible_paths_m[i,]
      for (j in 1:length(path)) {
        if (words_per_level[j] == 0) {
          next
        }
        words_ind_for_topic <- 1:ncol(word_assignment_counts)
        zero_word_ind = which(word_assignment_counts[j,] == 0)
        if (length(zero_word_ind) > 0) {
          words_ind_for_topic <- words_ind_for_topic[-zero_word_ind]
        }
        prob_numer = sum(sapply(words_ind_for_topic, function(w_i) {
          sum(log(topic_word_count_m[path[j], w_i] + a_beta[w_i] +
                    0:(word_assignment_counts[j, w_i] - 1)))
        }))
        prob_denom = sum(log(sum(topic_word_count_m[path[j],]) + sum(a_beta) + 
                               0:(words_per_level[j] - 1)))
        path_probs[i] = path_probs[i] + prob_numer - prob_denom
      }
    }
    new_path <- doc_possible_paths_m[sample(length(path_probs), 1, 
                                            prob = exp(path_probs - max(path_probs))),]
    # if (pooled && (length(subset_ind) < (length(cur_path) - 1))) {
    #   new_path[-subset_ind] <- sample(new_path[-subset_ind])
    # }
    
    #Update existing information
    topic_word_count_m[new_path,] <- 
      topic_word_count_m[new_path,] + word_assignment_counts
    total_topic_word_count[new_path] <- total_topic_word_count[new_path] + 
      total_level_word_count_m[d,]
    doc_path_assignment_m[d,] <- new_path
  }
  return(list("topic_word_count_m" = topic_word_count_m,
              "total_topic_word_count" = total_topic_word_count,
              "doc_path_assignment" = doc_path_assignment_m))
}

sample_c_for_all_level_hier <- function(level_word_count_m_list,
                                        total_level_word_count_m,
                                        topic_word_count_m, total_topic_word_count,
                                        doc_path_assignment_m,
                                        a_beta, possible_paths_m, J, K) {
  for (d in 1:nrow(doc_path_assignment_m)) {
    #Remove existing information
    cur_path = doc_path_assignment_m[d,]
    word_assignment_counts <- t(level_word_count_m_list[d,,])
    
    #word_assignment_level_counts <- table(topic_word_assignment_list[[d]][, c(3, 1)])
    topic_word_count_m[cur_path,] <-
      topic_word_count_m[cur_path,] - word_assignment_counts
    total_topic_word_count[cur_path] <- total_topic_word_count[cur_path] - 
      total_level_word_count_m[d,]
    words_per_level <- total_level_word_count_m[d,]
    topic_start <- 1
    topic_end <- 1
    log_topic_probs <- rep(0, nrow(topic_word_count_m))
    for (j in 1:J) {
      for (k in topic_start:topic_end) {
        if (words_per_level[j] == 0) {
          next
        }
        words_ind_for_topic <- 1:ncol(word_assignment_counts)
        zero_word_ind = which(word_assignment_counts[j,] == 0)
        if (length(zero_word_ind) > 0) {
          words_ind_for_topic <- words_ind_for_topic[-zero_word_ind]
        }
        prob_numer = sum(sapply(words_ind_for_topic, function(w_i) {
          sum(log(topic_word_count_m[k, w_i] + a_beta[w_i] +
                    0:(word_assignment_counts[j, w_i] - 1)))
        }))
        prob_denom = sum(log(total_topic_word_count[k] + sum(a_beta) + 
                               0:(words_per_level[j] - 1)))
        log_topic_probs[k] = log_topic_probs[k] + prob_numer - prob_denom
      }
      topic_start = topic_end + 1
      topic_end = K^j + topic_end
    }
    log_topic_probs = log_topic_probs
    path_probs <- apply(possible_paths_m, 1, function(path) {
      sum(log_topic_probs[path])
    })
    new_path <- possible_paths_m[sample(length(path_probs), 1, 
                                 prob = exp(path_probs - max(path_probs))),]
    
    #Update existing information
    topic_word_count_m[new_path,] <- 
      topic_word_count_m[new_path,] + word_assignment_counts
    total_topic_word_count[new_path] <- total_topic_word_count[new_path] + 
      total_level_word_count_m[d,]
    doc_path_assignment_m[d,] <- new_path
  }
  return(list("topic_word_count_m" = topic_word_count_m,
              "total_topic_word_count" = total_topic_word_count,
              "doc_path_assignment" = doc_path_assignment_m))
}

initialize_data <- function(doc_word_list, num_words, J, K, a_beta, a_theta, 
                            pooled = F, num_topics = sum(K^(0:(J - 1)))) {
  #Create possible paths
  if (pooled) {
    possible_paths_m <- as.matrix(generate_topic_permutations(num_topics))
  } else {
    possible_paths_m <- create_possible_topic_path_m(J, K)
  }
  doc_path_assignment_m <- possible_paths_m[
    sample(nrow(possible_paths_m), length(doc_word_list), replace = T),]
  topic_word_assignment_list <- list()
  topic_word_count_m <- matrix(0, nrow = num_topics, ncol = num_words)
  for (d in 1:length(doc_word_list)) {
    levels <- sample(ncol(possible_paths_m), 
                     length(doc_word_list[[d]]), replace = T)
    topic_word_doc <- as.data.frame(cbind(doc_word_list[[d]], 
                                          doc_path_assignment_m[d, levels], levels))
    topic_word_assignment_list <- c(topic_word_assignment_list, list(topic_word_doc))
    
    topic_word_counts <- table(topic_word_doc[, c(2, 1)])
    topic_word_count_m[as.integer(rownames(topic_word_counts)), 
                       as.integer(colnames(topic_word_counts))] <-
      topic_word_count_m[as.integer(rownames(topic_word_counts)), 
                         as.integer(colnames(topic_word_counts))] + topic_word_counts
  }
  return(list("topic_word_count" = topic_word_count_m,
              "topic_word_assignment" = topic_word_assignment_list,
              "doc_path_assignment" = doc_path_assignment_m, 
              "possible_paths_m" = possible_paths_m))
}

initialize_data_all_levels <- function(doc_word_list, num_words, J, K, a_beta, 
                                       a_theta, pooled = F, 
                                       num_topics = sum(K^(0:(J - 1)))) {
  #Create possible paths
  if (pooled) {
    possible_paths_m <- as.matrix(generate_topic_permutations(num_topics))
    level_word_count_m_list <- lapply(1:num_topics, function(j) {
      matrix(0, nrow = length(doc_word_list), ncol = num_words)})
    total_level_word_count_m <- matrix(0, nrow = length(doc_word_list), 
                                       ncol = num_topics)
  } else {
    possible_paths_m <- create_possible_topic_path_m(J, K)
    level_word_count_m_list <- lapply(1:J, function(j) {
      matrix(0, nrow = length(doc_word_list), ncol = num_words)})
    total_level_word_count_m <- matrix(0, nrow = length(doc_word_list), ncol = J)
  }
  doc_path_assignment_m <- possible_paths_m[
    sample(nrow(possible_paths_m), length(doc_word_list), replace = T),]
  topic_word_count_m <- matrix(0, nrow = num_topics, ncol = num_words)
  total_topic_word_count <- rep(0, num_topics)
  for (d in 1:length(doc_word_list)) {
    levels <- sample(ncol(possible_paths_m), length(doc_word_list[[d]]), replace = T)
    level_word_doc <- as.data.frame(cbind(doc_word_list[[d]], levels))
    word_counts <- table(level_word_doc[, c(2, 1)])
    level_word_counts <- matrix(0, nrow = ncol(possible_paths_m), ncol = num_words)
    level_word_counts[which(1:ncol(possible_paths_m) %in% as.integer(rownames(word_counts))),
                      which(1:num_words %in% as.integer(colnames(word_counts)))] <- 
      word_counts
    for (l in 1:ncol(possible_paths_m)) {
      level_word_count_m_list[[l]][d,] <- level_word_counts[l,]
    }
    total_level_word_count_m[d,] <- rowSums(level_word_counts)
    
    topic_word_count_m[doc_path_assignment_m[d,],] <-
      topic_word_count_m[doc_path_assignment_m[d,],] + level_word_counts
    total_topic_word_count[doc_path_assignment_m[d,]] <-
      total_topic_word_count[doc_path_assignment_m[d,]] + total_level_word_count_m[d,]
  }
  return(list("level_word_count_m_list" = simplify2array(level_word_count_m_list),
              "total_level_word_count_m" = total_level_word_count_m,
              "topic_word_count_m" = topic_word_count_m,
              "total_topic_word_count" = total_topic_word_count,
              "doc_path_assignment" = doc_path_assignment_m, 
              "possible_paths_m" = possible_paths_m))
}

initialize_data_all_levels_vb <- function(vb_results, doc_word_list, num_words, 
                                          J, K, a_beta, a_theta, 
                                          pooled = F, num_topics = sum(K^(0:(J - 1)))) {
  D = nrow(vb_results)
  if (pooled) {
    possible_paths_m <- as.matrix(generate_topic_permutations(num_topics))
    level_word_count_m_list <- lapply(1:num_topics, function(j) {
      matrix(0, nrow = length(doc_word_list), ncol = num_words)})
    total_level_word_count_m <- matrix(0, nrow = length(doc_word_list), ncol = num_topics)
  } else {
    possible_paths_m <- create_possible_topic_path_m(J, K)
    level_word_count_m_list <- lapply(1:J, function(j) {
      matrix(0, nrow = length(doc_word_list), ncol = num_words)})
    total_level_word_count_m <- matrix(0, nrow = length(doc_word_list), ncol = J)
  }
  opt_path <- apply(vb_results$c_d_probs_m, 1, function(probs) {
    max_probs_ind <- which(probs == max(probs))
    if (length(max_probs_ind) == 1) {
      return(max_probs_ind)
    }
    sample(max_probs_ind, 1)
  })
  doc_path_assignment_m <- possible_paths_m[opt_path,]
  topic_word_count_m <- matrix(0, nrow = num_topics, ncol = num_words)
  total_topic_word_count <- rep(0, num_topics)
  for (d in 1:length(doc_word_list)) {
    levels <- apply(vb_results$z_d_n_probs_m[,1:length(doc_word_list[[d]]),d], 
                    2, function(probs) {
                      sample(J, 1, prob = probs)
                    })
    level_word_doc <- as.data.frame(cbind(doc_word_list[[d]], levels))
    word_counts <- table(level_word_doc[, c(2, 1)])
    level_word_counts <- matrix(0, nrow = ncol(possible_paths_m), ncol = num_words)
    level_word_counts[which(1:ncol(possible_paths_m) %in% as.integer(rownames(word_counts))),
                      which(1:num_words %in% as.integer(colnames(word_counts)))] <- 
      word_counts
    for (l in 1:ncol(possible_paths_m)) {
      level_word_count_m_list[[l]][d,] <- level_word_counts[l,]
    }
    total_level_word_count_m[d,] <- rowSums(level_word_counts)
    
    topic_word_count_m[doc_path_assignment_m[d,],] <-
      topic_word_count_m[doc_path_assignment_m[d,],] + level_word_counts
    total_topic_word_count[doc_path_assignment_m[d,]] <-
      total_topic_word_count[doc_path_assignment_m[d,]] + total_level_word_count_m[d,]
  }
  return(list("level_word_count_m_list" = simplify2array(level_word_count_m_list),
              "total_level_word_count_m" = total_level_word_count_m,
              "topic_word_count_m" = topic_word_count_m,
              "total_topic_word_count" = total_topic_word_count,
              "doc_path_assignment" = doc_path_assignment_m, 
              "possible_paths_m" = possible_paths_m))
}

sample_lda_model <- function(c_results, possible_paths_m, num_iter, a_beta, a_theta, 
                             pooled = F, keep_all = T) {
  sampling_results <- list()
  for (i in 1:num_iter) {
    print(i)
    z_results <- sample_z(c_results$topic_word_count, c_results$topic_word_assignment, 
                          c_results$doc_path_assignment, a_beta, a_theta)
    print("Done z")
    c_results <- sample_c(z_results$topic_word_count, 
                          z_results$topic_word_assignment, 
                          c_results$doc_path_assignment, a_beta, possible_paths_m,
                          pooled = pooled)
    print("Done c")
    if (keep_all || (i >= num_iter / 2)) {
      sampling_results <- c(sampling_results, list(c_results))
    }
  }
  return(sampling_results)
  
}

sample_lda_model_for_all_level <- function(results, possible_paths_m, 
                                           num_iter, a_beta, a_theta, 
                                           J = 1, K = 1, pooled = F, 
                                           keep_all = T) {
  sampling_results <- list()
  for (i in 1:num_iter) {
    print(i)
    z_results <- sample_z_for_all_levels(
       results$level_word_count_m_list,
       results$total_level_word_count_m,
       results$topic_word_count_m, results$total_topic_word_count,
       results$doc_path_assignment,
       a_beta, a_theta)
    results$level_word_count_m_list <- z_results$level_word_count_m_list
    results$total_level_word_count_m <- z_results$total_level_word_count_m 
    results$topic_word_count_m <- z_results$topic_word_count_m
    results$total_topic_word_count <- z_results$total_topic_word_count
    print("Done z")
    remove(z_results)
    
    if (pooled) {
      c_results <- sample_c_for_all_level(
        results$level_word_count_m_list, results$total_level_word_count_m,
        results$topic_word_count_m, results$total_topic_word_count,
        results$doc_path_assignment, a_beta, possible_paths_m, pooled = T)
    } else {
      c_results <- sample_c_for_all_level_hier(
        results$level_word_count_m_list, results$total_level_word_count_m,
        results$topic_word_count_m, results$total_topic_word_count,
        results$doc_path_assignment, a_beta, possible_paths_m, J, K)
    }
    results$topic_word_count_m <- c_results$topic_word_count_m
    results$total_topic_word_count <- c_results$total_topic_word_count
    results$doc_path_assignment <- c_results$doc_path_assignment
    print("Done c")
    remove(c_results)
    
    if (keep_all || (i >= num_iter / 2)) {
      sampling_results <- c(sampling_results, list(results))
    }
  }
  return(sampling_results)
}

sample_hier_lda <- function(doc_word_list, num_words, J, K, 
                            a_beta, a_theta, num_iter = 1000,
                            num_chains = 1, keep_all = T) {
  c_results <- initialize_data(doc_word_list, num_words, J, K, a_beta, a_theta)
  possible_paths_m <- c_results$possible_paths_m
  if (num_chains == 1) {
    return(sample_lda_model(c_results, possible_paths_m, num_iter, 
                            a_beta, a_theta, keep_all = keep_all))
  } else {
    return(mclapply(1:num_chains, function(i) {
      sample_lda_model(c_results, possible_paths_m, num_iter, 
                       a_beta, a_theta, keep_all = keep_all)
    }, mc.cores = num_chains))
  }
}

sample_hier_lda_levels <- function(doc_word_list, num_words, J, K, 
                                   a_beta, a_theta, num_iter = 1000,
                                   num_chains = 1, keep_all = T,
                                   vb_results = NULL) {
  if (is.null(vb_results)) {
    results <- initialize_data_all_levels(doc_word_list, num_words,
                                          J, K, a_beta, a_theta)
  } else {
    results <- initialize_data_all_levels_vb(vb_results, doc_word_list,
                                             num_words, J, K, a_beta, a_theta)
  }
  possible_paths_m <- results$possible_paths_m
  results <- results[names(results) != "possible_paths_m"]
  if (num_chains == 1) {
    return(sample_lda_model_for_all_level(
      results, possible_paths_m, num_iter, a_beta, a_theta, keep_all = keep_all,
      J = J, K = K))
  } else {
    return(mclapply(1:num_chains, function(i) {
      sample_lda_model_for_all_level(
        results, possible_paths_m, num_iter, a_beta, a_theta,
        keep_all = keep_all, J = J, K = K)
    }, mc.cores = num_chains))
  }
}

sample_pooled_lda <- function(doc_word_list, num_words, num_topics, 
                              a_beta, a_order, num_iter = 1000, 
                              num_chains = 1, keep_all = T) {
  c_results <- initialize_data(doc_word_list, num_words, 1, 1, a_beta, a_order, 
                               pooled = T, num_topics = num_topics)
  possible_paths_m <- c_results$possible_paths_m
  if (num_chains == 1) {
    return(sample_lda_model(c_results, possible_paths_m, num_iter, 
                            a_beta, a_order, pooled = T, keep_all = keep_all))
  } else {
    return(mclapply(1:num_chains, function(i) {
      sample_lda_model(c_results, possible_paths_m, num_iter, 
                       a_beta, a_order, pooled = T, keep_all = keep_all)
    }, mc.cores = num_chains))
  }
  return(sampling_results)
}


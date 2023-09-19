library(kernlab)
source("tbm_helper_functions.R")
source("vlad.R")
source("gdm.R")

identify_group_to_start <- function(groups_available, available_ind, dist_m) {
  test_inds <- unlist(available_ind[groups_available])
  max_min = Inf
  max_group_ind = -1
  for (group_ind in groups_available) {
    group_max_min_tmp = 
      max(dist_m[available_ind[[group_ind]], unlist(available_ind[
        groups_available[groups_available != group_ind]])])
    if (group_max_min_tmp < max_min) {
      max_min = group_max_min_tmp
      max_group_ind = group_ind
    }
  }
  return(max_group_ind)
}

find_closest_neighbors_for_level <- function(start_ind, num_members,
                                             available_ind_list, dist_m, J,
                                             used_groups = NULL) {
  start_group = (start_ind - 1) %/% J + 1
  closest_neighbors = c(start_ind)
  groups_involved = c(start_group)
  overall_max_dist = 0
  for (i in 1:(num_members - 1)) {
    choices <- unlist(available_ind_list[-c(groups_involved, used_groups)])
    min_dist = Inf
    interested_ind = -1
    for (j in choices) {
      new_min_dist = max(dist_m[closest_neighbors, j])
      if (new_min_dist < min_dist) {
        min_dist = new_min_dist
        interested_ind = j
      }
    }
    if (overall_max_dist < min_dist) {
      overall_max_dist = min_dist
    }
    closest_neighbors <- c(closest_neighbors, interested_ind)
    groups_involved <- c(groups_involved, (interested_ind - 1) %/% J + 1)
  }
  return(list("neighbors" = closest_neighbors,
              "groups_used" = groups_involved,
              "max_dist" = overall_max_dist))
}

find_closet_neighbors_for_subset <- function(gdm_results_m, dist_m, J, K, 
                                             available_ind, subset_groups, level) {
  num_paths = K^(J - 1)
  groups_available = subset_groups
  num_members = K^(J - level)
  if (length(subset_groups) != num_paths) {
    all_groups = 1:length(available_ind)
    used_groups = all_groups[!(all_groups %in% groups_available)]
  } else {
    used_groups = c()
  }
  merge_groups_info <- list()
  repeat {
    min_overall_dist = Inf
    neighbors_list <- c()
    groups_used_list <- c()
    group_start_ind <- 
      identify_group_to_start(groups_available, available_ind, dist_m)
    for (start_ind in available_ind[[group_start_ind]]) {
      closest_neighbor_info <- find_closest_neighbors_for_level(
        start_ind, num_members, available_ind, dist_m, J, used_groups)
      # print(closest_neighbor_info$max_dist)
      if (closest_neighbor_info$max_dist < min_overall_dist) {
        min_overall_dist = closest_neighbor_info$max_dist
        neighbors_list = closest_neighbor_info$neighbors
        groups_used_list = closest_neighbor_info$groups_used
      }
    }
    # print(groups_used_list)
    # print(neighbors_list)
    merge_groups_info <- c(merge_groups_info, list(neighbors_list))
    groups_available <- groups_available[!(groups_available %in% 
                                             groups_used_list)]
    used_groups <- c(used_groups, groups_used_list)
    for (j in 1:length(groups_used_list)) {
      tmp = available_ind[[groups_used_list[j]]]
      tmp <- tmp[tmp != neighbors_list[j]]
      available_ind[[groups_used_list[j]]] <- tmp
    }
    if (length(groups_available) == 0) {
      break
    }
  }
  return(list(merge_groups_info, available_ind))
  
} 

find_closest_neighbors_ordered <- function(gdm_results_m, J, K) {
  num_paths = K^(J - 1)
  dist_m <- as.matrix(dist(gdm_results_m))
  available_ind <- lapply(1:num_paths, function(i) {J * (i - 1) + 1:J})
    
  neighbors_info <- find_closet_neighbors_for_subset(
    gdm_results_m, dist_m, J, K, available_ind, 1:num_paths, 1)
  merge_groups_info <- neighbors_info[[1]]
  available_ind <- neighbors_info[[2]]
  subset_groups <- list(1:num_paths)
  for (j in 2:(J - 1)) {
    new_subset_groups <- list()
    for (k in 1:length(subset_groups)) {
      neighbors_info <- find_closet_neighbors_for_subset(
        gdm_results_m, dist_m, J, K, available_ind, subset_groups[[k]], j)
      merge_groups_info <- c(merge_groups_info, neighbors_info[[1]])
      available_ind <- neighbors_info[[2]]
      new_subset_groups <- 
        c(new_subset_groups, 
          lapply(neighbors_info[[1]], function(merge_ind) {
            (merge_ind - 1) %/% J + 1
          }))
    }
    subset_groups <- new_subset_groups
    print(subset_groups)
  }
  return(c(merge_groups_info, available_ind[unlist(subset_groups)]))
}

merge_topics_hclust <- function(gdm_results_m, J, K, num_words) {
  # merge_plan <- group_hier_clust(gdm_results_m, J, K)
  # merge_plan <- find_closest_neighbors(gdm_results_m, J, K)
  merge_plan <- find_closest_neighbors_ordered(gdm_results_m, J, K)
  print(merge_plan)
  gdm_topics_m <- matrix(0, nrow = sum(K^(0:(J - 1))), ncol = num_words)
  for (i in 1:length(merge_plan)) {
    gdm_topics_m[i,] <- colMeans(gdm_results_m[merge_plan[[i]], , drop = F])
  }
  return(gdm_topics_m)
}

get_topics_alg <- function(level_labels, prop_doc_word_m, num_words, J, alg) {
  
  unique_level_labels <- sort(unique(level_labels))
  if (alg == "GDM") {
    gdm_results <- lapply(unique_level_labels,
                          function(label) {
                            run_tuned_gdm(as.data.frame(
                              prop_doc_word_m[which(level_labels == label),, drop = F]),
                              1:num_words, J, tune = T)
                          })
    gdm_results_m <- t(do.call(cbind, gdm_results))
  }
  if (alg == "hybrid") {
    gdm_results <- lapply(unique_level_labels,
                          function(label) {
                            tune_vlad(prop_doc_word_m[which(level_labels == label),, drop = F], J)
                          })
    gdm_results_m <- do.call(rbind, gdm_results)
  }
  if (alg == "VLAD") {
    gdm_results <- lapply(unique_level_labels,
                          function(label) {
                            vlad(prop_doc_word_m[which(level_labels == label),, drop = F], J)
                          })
    gdm_results_m <- do.call(rbind, gdm_results)
  }
  return(gdm_results_m)
}         

hier_gdm <- function(prop_doc_word_m, num_words, J, K,
                     start_labels = NULL, alg = "hybrid") {
  dist_prev <- Inf
  level_labels <- start_labels
  if (is.null(level_labels)) {
    #level_labels <- kmeans(prop_doc_word_m, K^(J - 1))$cluster
    # level_labels <- sample(K^(J - 1), nrow(prop_doc_word_m), replace = T)
    # level_labels <- specClust(prop_doc_word_m, K^(J - 1))$cluster
    level_labels <- specc(prop_doc_word_m, K^(J - 1))@.Data
  }
  print(table(level_labels))
  topic_counts <- list()
  possible_paths_m <- create_possible_topic_path_m(J, K)
  i = 1
  gdm_topics_m_prev <- matrix(1, nrow = max(possible_paths_m), ncol = num_words)
  symplecial_list <- c(Inf)
  prev_label_count <- rep(-1, K^(J - 1))
  next_label_count <- prev_label_count
  repeat {
    print(i)
    unique_level_labels <- sort(unique(level_labels))

    #Get topics
    gdm_results_m <- get_topics_alg(level_labels, prop_doc_word_m, num_words, J, alg)
    group_centers <- t(sapply(unique_level_labels, function(label) {
      colMeans(prop_doc_word_m[which(level_labels == label),, drop = F])
    }))

    #Merge topics
    gdm_topics_m_next <- merge_topics_hclust(gdm_results_m, J, K, num_words)
    # gdm_topics_m_next <- merge_topics_hclust_intersection(gdm_results_m, group_centers, J, K, num_words)
  
    #Reassign topics
    level_dist <- t(apply(prop_doc_word_m, 1, function(doc) {
      apply(possible_paths_m, 1, function(path) {
        proj_on_s(t(gdm_topics_m_next[path,]), doc, 
                  length(path), 1:length(path),
                  dist = T)
      })
    }))
    level_labels <- apply(level_dist, 1, which.min)
    symplecial_dist <- mean(mapply(function(j, k) {level_dist[j, k]}, 
                                   1:nrow(level_dist), level_labels))
    print(table(level_labels))
    print(c("Symplecial distance: ", symplecial_dist))
    symplecial_list <- c(symplecial_list, symplecial_dist)
    i = i + 1
    next_label_count = table(level_labels)
    
    # #Randomly re-assign labels if any have less points than the depth
    # level_counts <- table(level_labels)
    # missing_levels <- which(level_counts < J)
    # if (length(missing_levels) != 0) {
    #   missing_levels
    # }
    

    #Calculate shift in topics
    # dist_next = find_closest_topics_dist(gdm_topics_m_prev, 
    #                                      gdm_topics_m_next,
    #                                      J, K)
    # print(dist_next)
    # print(dist_next - dist_prev)
    if (abs(symplecial_list[i - 1] - symplecial_list[i]) < 1e-4 || all(next_label_count == prev_label_count) || any(table(level_labels) <= J)) {
      break
    }
    prev_label_count = next_label_count
    prev_level_labels <- level_labels
    prev_level_dist <- level_dist
    # if (abs(dist_prev - dist_next) < 1e-4) {
    #   break
    # }
    # dist_prev = dist_next
    gdm_topics_m_prev = gdm_topics_m_next
  }
  if (symplecial_list[i - 1] < symplecial_list[i]) {
    return(list("level_labels" = prev_level_labels,
                "level_label_m" = possible_paths_m[prev_level_labels,],
                "topics_m" = gdm_topics_m_prev,
                "symplecial_list" = symplecial_list,
                "level_dist" = prev_level_dist))    
  }
  return(list("level_labels" = level_labels,
              "level_label_m" = possible_paths_m[level_labels,],
              "topics_m" = gdm_topics_m_next,
              "symplecial_list" = symplecial_list,
              "level_dist" = level_dist))
}


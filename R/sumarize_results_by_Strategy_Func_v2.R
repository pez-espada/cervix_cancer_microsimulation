summarize_results_by_Strategy_v2 <- function(strategy, results_list, numb_of_sims) {
  library(dplyr)
  library(tibble)
  
  result <- list()
  Strategy_name <- strategy
  result[[Strategy_name]] <- list()
  
  # Remove unwanted large matrices before stacking
  results_list <- lapply(results_list, function(df) {
    df$tc_disc <- NULL
    df$tc_undisc <- NULL
    df$te_disc <- NULL
    df$te_undisc <- NULL
    df$Tot_Trans_per_t <- NULL
    return(df)
  })
  
  sim <- results_list
  names_sim <- names(sim[[1]])
  
  for (name_level_of_sim in names_sim) {
    result[[Strategy_name]][[name_level_of_sim]] <- bind_rows(lapply(seq_along(sim), function(i) {
      df <- as.data.frame(sim[[i]][[name_level_of_sim]])
      
      # ✅Prevent duplicate 'sim' column conflicts
      if ("sim" %in% colnames(df)) {
        df <- df %>% dplyr::select(-sim)
      }
      df$sim <- i
      
      if (name_level_of_sim != "TR") {
        df$row_names <- as.numeric(rownames(df))
        rownames(df) <- NULL
        df <- df[, c("sim", "row_names", setdiff(colnames(df), c("sim", "row_names")))]
      } else {
        rownames(df) <- NULL
        df$cycle <- seq_len(nrow(df))
      }
      
      return(df)
    }))
    
    # ✅ Rename misnamed last columns if necessary
    if (any(names(result[[Strategy_name]][[name_level_of_sim]]) %in% c("V1", "sim[[i]][[name_level_of_sim]]"))) {
      colnames(result[[Strategy_name]][[name_level_of_sim]])[ncol(result[[Strategy_name]][[name_level_of_sim]])] <- name_level_of_sim
    }
    
    # ✅ Clean up row_names for transition matrices
    if (any(endsWith(names(result[[Strategy_name]][[name_level_of_sim]]), "_disc")) |
        any(endsWith(names(result[[Strategy_name]][[name_level_of_sim]]), "_undisc"))) {
      result[[Strategy_name]][[name_level_of_sim]]$row_names <- NULL
    }
  }
  
  # Optional cleanup wrapper for returning only data.frames
  new_result <- list()
  for (name_level_of_sim in names_sim) {
    new_result[[Strategy_name]][[name_level_of_sim]] <- as.data.frame(result[[Strategy_name]][[name_level_of_sim]])
  }
  
  return(new_result)
}

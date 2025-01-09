#########################
### Second Attempt ######


# summarize_results_by_Strategy <- function(results_list, numb_of_sims) {
#   library(dplyr)
#   
#   result <- list()
#   Strategy_name <- "No Intervention"
#   result[[Strategy_name]] <- list()
#   
#   # Remove unwanted columns from results_list
#   results_list <- lapply(results_list, function(df) {
#     df$tc_disc <- NULL
#     df$tc_undisc <- NULL
#     df$te_disc <- NULL
#     df$te_undisc <- NULL
#     df$Tot_Trans_per_t <- NULL
#     return(df)
#   })
#   
#   sim <- results_list
#   names_sim <- names(sim[[1]])
#   
#   for (name_level_of_sim in names_sim) {
#     result[[Strategy_name]][[name_level_of_sim]] <- bind_rows(lapply(seq_along(sim), function(i) {
#       df <- as.data.frame(sim[[i]][[name_level_of_sim]])
#       
#       # Remove pre-existing "sim" column if it exists
#       if ("sim" %in% colnames(df)) {
#         df <- df %>% select(-sim)
#       }
#       
#       # Add a new "sim" column
#       df <- df %>% mutate(sim = i)
#       
#       if (name_level_of_sim != "TR") {
#         df <- df %>% mutate(row_names = as.numeric(rownames(df)))
#         rownames(df) <- NULL 
#         df <- df[, c("sim", "row_names", setdiff(colnames(df), c("sim", "row_names")))]
#       } else {
#         rownames(df) <- NULL 
#         df <- df %>% mutate(cycle = seq(1:nrow(df)))
#       }
#       
#       return(df)
#     }))
#     
#     # Ensure column names are consistent
#     if (sum(names(result[[Strategy_name]][[name_level_of_sim]]) == "V1") > 0 || 
#         sum(names(result[[Strategy_name]][[name_level_of_sim]]) == "sim[[i]][[name_level_of_sim]]") > 0) {
#       colnames(result[[Strategy_name]][[name_level_of_sim]])[ncol(result[[Strategy_name]][[name_level_of_sim]])] <- name_level_of_sim
#     }
#     
#     if (sum(endsWith(names(result[[Strategy_name]][[name_level_of_sim]]), "_disc")) > 0 || 
#         sum(endsWith(names(result[[Strategy_name]][[name_level_of_sim]]), "_undisc")) > 0) {
#       result[[Strategy_name]][[name_level_of_sim]]["row_names"] <- NULL
#     }
#   }
#   
#   # Convert lists to data frames
#   new_result <- list()
#   for (name_level_of_sim in names_sim) {
#     new_result[[Strategy_name]][[name_level_of_sim]] <- as.data.frame(result[[Strategy_name]][[name_level_of_sim]])
#   }
#   
#   return(new_result)
# }




#########################
### Second Attempt ######
#summarize_results_by_Strategy <- function(results_list, numb_of_sims) {
#  library(dplyr)
#  
#  result <- list()
#  Strategy_name <- "No Intervention"
#  result[[Strategy_name]] <- list()
#  
#  # Clean up unwanted columns
#  results_list <- lapply(results_list, function(df) {
#    df$tc_disc <- NULL
#    df$tc_undisc <- NULL
#    df$te_disc <- NULL
#    df$te_undisc <- NULL
#    df$Tot_Trans_per_t <- NULL
#    return(df)
#  })
#  
#  sim <- results_list
#  names_sim <- names(sim[[1]])
#  
#  for (name_level_of_sim in names_sim) {
#    result[[Strategy_name]][[name_level_of_sim]] <- bind_rows(lapply(seq_along(sim), function(i) {
#      df <- as.data.frame(sim[[i]][[name_level_of_sim]])
#      
#      # Remove pre-existing "sim" column if it exists
#      if ("sim" %in% colnames(df)) {
#        df <- df %>% select(-sim)
#      }
#      
#      # Add a new "sim" column
#      df <- df %>% mutate(sim = i)
#      
#      if (name_level_of_sim != "TR") {
#        df <- df %>% mutate(row_names = as.numeric(rownames(df)))
#        rownames(df) <- NULL 
#        df <- df[, c("sim", "row_names", setdiff(colnames(df), c("sim", "row_names")))]
#      } else {
#        rownames(df) <- NULL 
#        df <- df %>% mutate(cycle = seq(1:nrow(df)))
#      }
#      
#      return(df)
#    }))
#    
#    # Ensure column names are consistent and avoid non-existent column renaming
#    column_names <- names(result[[Strategy_name]][[name_level_of_sim]])
#    
#    if ("V1" %in% column_names) {
#      colnames(result[[Strategy_name]][[name_level_of_sim]])[which(column_names == "V1")] <- name_level_of_sim
#    }
#  }
#  
#  # Convert lists to data frames
#  new_result <- list()
#  for (name_level_of_sim in names_sim) {
#    new_result[[Strategy_name]][[name_level_of_sim]] <- as.data.frame(result[[Strategy_name]][[name_level_of_sim]])
#  }
#  
#  return(new_result)
#}

##########################
##### Third Attempt ### 
#################################################################################
#summarize_results_by_Strategy <- function(results_list, numb_of_sims) {
#  library(dplyr)
#  
#  result <- list()
#  Strategy_name <- "No Intervention"
#  result[[Strategy_name]] <- list()
#  
#  # Clean up unwanted columns dynamically
#  results_list <- lapply(results_list, function(df) {
#    columns_to_remove <- c("tc_disc", "tc_undisc", "te_disc", "te_undisc", "Tot_Trans_per_t", "seed")
#    df <- df[, !(names(df) %in% columns_to_remove), drop = FALSE]
#    return(df)
#  })
#  
#  sim <- results_list
#  names_sim <- names(sim[[1]])
#  
#  for (name_level_of_sim in names_sim) {
#    result[[Strategy_name]][[name_level_of_sim]] <- bind_rows(lapply(seq_along(sim), function(i) {
#      df <- as.data.frame(sim[[i]][[name_level_of_sim]])
#      
#      # Remove pre-existing "sim" column if it exists
#      if ("sim" %in% colnames(df)) {
#        df <- df %>% select(-sim)
#      }
#      
#      # Add a new "sim" column
#      df <- df %>% mutate(sim = i)
#      
#      if (name_level_of_sim != "TR") {
#        # Add row names column if row names exist
#        if (!is.null(rownames(df))) {
#          df <- df %>% mutate(row_names = as.numeric(rownames(df)))
#        }
#        rownames(df) <- NULL
#        df <- df[, c("sim", setdiff(colnames(df), "row_names"))]
#      } else {
#        rownames(df) <- NULL
#        df <- df %>% mutate(cycle = seq(1:nrow(df)))
#      }
#      
#      return(df)
#    }))
#    
#    # Ensure column names are consistent
#    column_names <- names(result[[Strategy_name]][[name_level_of_sim]])
#    if ("V1" %in% column_names) {
#      colnames(result[[Strategy_name]][[name_level_of_sim]])[which(column_names == "V1")] <- name_level_of_sim
#    }
#  }
#  
#  # Convert lists to data frames
#  new_result <- list()
#  for (name_level_of_sim in names_sim) {
#    new_result[[Strategy_name]][[name_level_of_sim]] <- as.data.frame(result[[Strategy_name]][[name_level_of_sim]])
#  }
#  
#  return(new_result)
#}


#########################
#### Fourt Attempt ### 
################################################################################
summarize_results_by_Strategy <- function(results_list, numb_of_sims) {
  library(dplyr)
  
  result <- list()
  Strategy_name <- "No Intervention"
  result[[Strategy_name]] <- list()
  
  # Clean up unwanted columns dynamically
  results_list <- lapply(results_list, function(item) {
    # Ensure the item is a data frame
    if (!is.data.frame(item)) {
      stop("Each element of results_list must be a data frame.")
    }
    
    columns_to_remove <- c("tc_disc", "tc_undisc", "te_disc", "te_undisc", "Tot_Trans_per_t", "seed")
    item <- item[, !(names(item) %in% columns_to_remove), drop = FALSE]
    return(item)
  })
  
  sim <- results_list
  names_sim <- names(sim[[1]])
  
  for (name_level_of_sim in names_sim) {
    result[[Strategy_name]][[name_level_of_sim]] <- bind_rows(lapply(seq_along(sim), function(i) {
      df <- as.data.frame(sim[[i]][[name_level_of_sim]])
      
      # Remove pre-existing "sim" column if it exists
      if ("sim" %in% colnames(df)) {
        df <- df %>% select(-sim)
      }
      
      # Add a new "sim" column
      df <- df %>% mutate(sim = i)
      
      if (name_level_of_sim != "TR") {
        # Add row names column if row names exist
        if (!is.null(rownames(df))) {
          df <- df %>% mutate(row_names = as.numeric(rownames(df)))
        }
        rownames(df) <- NULL
        df <- df[, c("sim", setdiff(colnames(df), "row_names"))]
      } else {
        rownames(df) <- NULL
        df <- df %>% mutate(cycle = seq(1:nrow(df)))
      }
      
      return(df)
    }))
    
    # Ensure column names are consistent
    column_names <- names(result[[Strategy_name]][[name_level_of_sim]])
    if ("V1" %in% column_names) {
      colnames(result[[Strategy_name]][[name_level_of_sim]])[which(column_names == "V1")] <- name_level_of_sim
    }
  }
  
  # Convert lists to data frames
  new_result <- list()
  for (name_level_of_sim in names_sim) {
    new_result[[Strategy_name]][[name_level_of_sim]] <- as.data.frame(result[[Strategy_name]][[name_level_of_sim]])
  }
  
  return(new_result)
}


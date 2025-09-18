################################################################################
###################################
## Post-simulation Computations: ##
###################################
## This script compute post-processing for
## n_strat strategies times n_sim batch simulations
################################################################################
library(tidyverse)

# browser() # This pauses execution like a breakpoint

# Adding runtime execution time and some other parameters:
runtime <- comp.time %>% as_tibble() %>% `colnames<-`("runtime")
#sim_result[[1]]$runtime        <- runtime

# browser() # This pauses execution like a breakpoint

# For clarity sake, promote inner named elements to top level
sim_result <- setNames(
  lapply(sim_result, `[[`, 1),  # Get the inner object (first element of each sublist)
  sapply(sim_result, function(x) names(x)[1])  # Use the inner name (e.g., "STRATEGY A") as the new name
)

# browser() # This pauses execution like a breakpoint

sim_result$runtime  <- runtime #add runtime at the base lavel of the list

# browser() # This pauses execution like a breakpoint

################################################################################
################################################################################
for (strategy_name in names(sim_result[names(sim_result) != "runtime"])) {
  strategy_data <-  sim_result[[strategy_name]]
  cat("Processing:", strategy_name, "\n")
  
  # Processing code here
  
  #sim_result[[strategy_name]]$strategy       <- strategy
  sim_result[[strategy_name]]$strategy       <- strategy_name
  sim_result[[strategy_name]]$numb_of_sims   <- numb_of_sims
  sim_result[[strategy_name]]$numb_of_ind    <- n_i
  sim_result[[strategy_name]]$numb_of_cycles <- n_t
 
   
  
  sim_result[[strategy_name]]$seed <- sim_result[[strategy_name]]$seed %>% 
    dplyr::select(-c("seed", "row_names")) %>% 
    dplyr::rename("seed" = "sim[[i]][[name_level_of_sim]]")
  
  
  
  sim_result[[strategy_name]]$tc_hat_undisc <- sim_result[[strategy_name]]$tc_hat_undisc %>%
    dplyr::select(-c(tc_hat_undisc)) %>% 
    dplyr::rename("tc_hat_undisc" = "sim[[i]][[name_level_of_sim]]")
  
  sim_result[[strategy_name]]$tc_hat_disc <- sim_result[[strategy_name]]$tc_hat_disc %>%
    dplyr::select(-c(tc_hat_disc)) %>% 
    dplyr::rename("tc_hat_disc" = "sim[[i]][[name_level_of_sim]]")
  
  sim_result[[strategy_name]]$te_hat_undisc <- sim_result[[strategy_name]]$te_hat_undisc %>%
    dplyr::select(-c(te_hat_undisc)) %>% 
    dplyr::rename("te_hat_undisc" = "sim[[i]][[name_level_of_sim]]")
  
  sim_result[[strategy_name]]$te_hat_disc <- sim_result[[strategy_name]]$te_hat_disc %>%
    dplyr::select(-c(te_hat_disc)) %>% 
    dplyr::rename("te_hat_disc" = "sim[[i]][[name_level_of_sim]]")
  sim_result[[strategy_name]]$vacc_coverage <- vacc_coverage
  
  
  ################################################################################
  ## Prevalence is defined as number of infected divided by total alive 
  ## individuals for that cycle/time step
  ################################################################################
  mean_prevalence_func <- function(sim_stalked_result, my_Probs) {
    # Extract unique age intervals and ensure Larger doesn't exceed 84
    age_intervals <- my_Probs %>% 
      select(Lower, Larger) %>% 
      unique() %>% 
      arrange(Lower)
    
    # Create a vector of breaks for the intervals
    breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
    
    # Create labels for the intervals
    labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
    
    # Compute prevalence and average it by age intervals
    df <- sim_stalked_result$TR %>% 
      dplyr::select(everything()) %>% 
      dplyr::mutate(total_alive = H + HR.HPV.infection + CIN1 + CIN2 + CIN3 +
                      FIGO.I + FIGO.II + FIGO.III + FIGO.IV + Survival) %>%
      dplyr::mutate(prevalence = HR.HPV.infection / total_alive) %>%
      dplyr::mutate(age_interval = cut(age, breaks = breaks, 
                                       labels = labels, right = FALSE)) %>%
      dplyr::group_by(age_interval) %>%
      dplyr::summarise(prevalence = mean(prevalence, na.rm = TRUE)) %>%
      dplyr::ungroup()
    
    # Store the result in the list
    # sim_stalked_result[[strategy_name]]$mean_HPV_prevalence_per_age_interval <- df
    # return(sim_stalked_result)
     return(df)
  }
  ################################################################################
  
  sim_result[[strategy_name]]$mean_HPV_prevalence_per_age_interval <-
    mean_prevalence_func(sim_stalked_result = sim_result[[strategy_name]],
                         my_Probs = my_Probs)
 
   
  ################################################################################
  # Incidence is defined by the NEW number of individuals in the state of interest
  # in the time t divided by all the individuals in  the epidemiological "precedent" 
  # states at time t-1 (in the previous cycle). This is how is defined in the  the 
  # Markov model. We follow that definition to compare the micro sim and the Markov. 
  mean_incidence_func <- function(sim_stalked_result, state, my_Probs) {
    # Extract unique age intervals and ensure Larger doesn't exceed 84
    age_intervals <- my_Probs %>% 
      dplyr::select(Lower, Larger) %>% 
      unique() %>% 
      arrange(Lower)
    
    # Create a vector of the breaks for the intervals
    breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
    
    # Create labels for the intervals
    labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
    
    # Define all previous epi states
    all_previous_states <- v_n[1 : ( which(v_n == state) - 1)]
    
    my_incidence_df <- sim_stalked_result$TR
    
    # Define the new state for incidence calculation
    new_state <- paste0("new_", state)
    new_state_df <- sim_stalked_result[new_state][[1]]
    
    my_incidence_df <- 
      my_incidence_df %>% 
      dplyr::left_join(new_state_df %>%
                         dplyr::select(-c(sim.1, row_names)), 
                       by = c("age"), relationship = "many-to-many")
    
    # Find the column containing "->" and the state
    target_column <- grep(paste0("->", state), names(my_incidence_df), value = TRUE)
    
    # Calculate the incidence rate
    my_incidence_df <- my_incidence_df %>%
      dplyr::mutate(!!paste0("incidence_", state) := 
                      (get(target_column) /
                         dplyr::lag(rowSums(select(., all_of(all_previous_states))), 
                                    n=1, default = NA)) * 10^5 )
    
    # Ensure the 'rlang' library is available
    ensure_library("rlang")
    
    # Create the age intervals and compute mean incidence
    df <- my_incidence_df %>%
      dplyr::select(everything()) %>%
      dplyr::mutate(age_interval = 
                      cut(age, breaks = breaks, 
                          labels = labels, 
                          right = FALSE)) %>%
      group_by(age_interval) %>%
      summarise(!!paste0("mean_incidence_", state) := 
                  mean(!!sym(paste0("incidence_", state)), na.rm = TRUE)) %>%
      ungroup()
    
    ## Store the result in the list
    #sim_stalked_result[[1]][[paste0("mean_incidence_", state, "_per_age_interval")]] <- df
    #return(sim_stalked_result)
    return(df)
  }
  ##############################################################################
  
  # Computing incidences:
  incidence_states_to_compute <- c("CIN1", "CIN2", "CIN3") 
  
  # Initialize the result with the original structure
  #mean_incidence_result <- mean_prevalence_result
  
  # Apply the incidence function to each state and update the result structure
  for (my_state in incidence_states_to_compute) {
    #print(my_state)
   sim_result[[strategy_name]][[paste0("mean_incidence_", my_state, "_per_age_interval")]]<- 
      mean_incidence_func(sim_stalked_result = sim_result[[strategy_name]], 
                          state = my_state, my_Probs = my_Probs)
  }
  ##############################################################################
  
  ##############################################################################
  # Computing Cervix Cancer incidence:
  mean_CC_incidence_func <- function(sim_stalked_result, my_Probs) {
    
    age_intervals <- my_Probs %>% 
      dplyr::select(Lower, Larger) %>% 
      unique() %>% 
      arrange(Lower)
    
    # Create a vector of the breaks for the intervals
    breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
    
    # Create labels for the intervals
    labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
    sim_stalked_result$new_Cancer <- sim_stalked_result$new_Cancer #%>%
    # dplyr::select(-sim.1)
    
    # Compute prevalence and average it by age intervals
    df <- merge(sim_stalked_result$TR, 
                sim_stalked_result$new_Cancer, by = c("sim", "age")) %>%
      dplyr::rename(new_Cancer = `CIN3->FIGO.I`)
    
    df <- df %>% 
      #dplyr::select(sim, cycle, age, H, HR.HPV.infection) %>% 
      dplyr::select(everything()) %>% 
      # total_alive at the previous time step
      dplyr::mutate(total_alive_lagged = dplyr::lag(H + HR.HPV.infection + CIN1 + CIN2 + CIN3 +
                                                      FIGO.I + FIGO.II + FIGO.III + FIGO.IV + Survival, n=1)) %>%
      dplyr::mutate(CC_incidence = (new_Cancer / total_alive_lagged) * 10^5) %>% 
      dplyr::mutate(age_interval = cut(age, breaks = breaks, labels = labels, right = FALSE)) %>% 
      dplyr::group_by(age_interval) %>% 
      dplyr::summarise(CC_mean_incidence = mean(CC_incidence, na.rm = TRUE)) %>% 
      dplyr::ungroup()
    
    return(df)
    #sim_stalked_result[[1]]$mean_CC_incidence <- df
    #return(sim_stalked_result)
  }
  ##############################################################################
  
   sim_result[[strategy_name]]$mean_CC_incidence <-
    mean_CC_incidence_func(sim_stalked_result = sim_result[[strategy_name]],
                           my_Probs = my_Probs)  
  
  ############################################################################## 
  
  ##############################################################################
  # Computing Mortality:
  # A. Cancer-related deaths at certain age (cycle) / total alive at that age (cycle)
  # B. Cancer-unrelated deaths at certain age (cycle) / total alive at that age (cycle)
  
  ##############################################################################
  # A. Cancer-related Deaths per age
  mean_CC_mortality_func <- function(sim_stalked_result, my_Probs) {
    
    age_intervals <- my_Probs %>% 
      dplyr::select(Lower, Larger) %>% 
      unique() %>% 
      arrange(Lower)
    
    # Create a vector of the breaks for the intervals
    breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
    
    # Create labels for the intervals
    labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
    
    # Define the age range you want to keep
    age_range <- 10:84
    
    # Left join sim_result[[1]]$TR with sim_result[[1]]$new_CC_Death by age
    #df <- sim_result[[1]]$TR %>%
    df <- sim_stalked_result$TR %>%
      #left_join(sim_result[[1]]$new_CC_Death %>%
      left_join(sim_stalked_result$new_CC_Death %>%
                  dplyr::select(sim, age, CC_Death_per_t), 
                by = c("sim", "age"))  %>% #, relationship = "many-to-many") %>%
      
      # Filter for ages in the desired range
      dplyr::filter(age %in% age_range) %>%
      
      # Fill missing CC_Death_per_t with zeros (for cases where the age doesn't exist)
      dplyr::mutate(CC_Death_per_t = coalesce(CC_Death_per_t, 0)) %>%
      
      # Compute total_alive one previous time stpe / cylce:
      dplyr::mutate(total_alive_lagged =
                      dplyr::lag(H + HR.HPV.infection + CIN1 + CIN2 + CIN3 +
                                   FIGO.I + FIGO.II + FIGO.III + FIGO.IV + Survival), n=1) %>%
      
      # Compute CC_mortality based on CC_Death_per_t and total_alive
      dplyr::mutate(CC_mortality = (CC_Death_per_t / total_alive_lagged) * 10^5) %>%
      
      # Compute age intervals and average CC_mortality by age intervals
      dplyr::mutate(age_interval = cut(age, breaks = breaks, labels = labels, right = FALSE)) %>%
      dplyr::group_by(age_interval) %>%
      dplyr::summarise(CC_mean_mortality = mean(CC_mortality, na.rm = TRUE)) %>%
      dplyr::ungroup()
    
    #sim_stalked_result[[1]]$CC_mean_mortality <- df
    #return(sim_stalked_result)
    return(df)
  }
  ##############################################################################
  
  sim_result[[strategy_name]]$CC_mean_mortality <-
    mean_CC_mortality_func(sim_stalked_result = sim_result[[strategy_name]],
                           my_Probs = my_Probs)  
  

  
  ##############################################################################
  # A.2 Cancer-related Deaths (per differences) per age
  mean_CC_mortality_by_diff_func <- function(sim_stalked_result, my_Probs) {
    
    age_intervals <- my_Probs %>% 
      dplyr::select(Lower, Larger) %>% 
      unique() %>% 
      arrange(Lower)
    
    # Create a vector of the breaks for the intervals
    breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
    
    # Create labels for the intervals
    labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
    
    # Define the age range you want to keep
    age_range <- min(age_intervals$Lower):max(age_intervals$Larger) 
    
    # Left join sim_result[[1]]$TR with sim_result[[1]]$new_CC_Death by age
    #df <- sim_result[[1]]$TR %>%
    df <- sim_stalked_result$TR %>%
      left_join(sim_stalked_result$CC_Death_by_diff %>%
                  dplyr::select(sim, age, CC_Death_by_diff), 
                by = c("sim", "age", "CC_Death_by_diff"), 
                relationship = "many-to-many") %>%
      
      # Filter for ages in the desired range
      dplyr::filter(age %in% age_range) %>%
      
      # Fill missing CC_Death_per_t with zeros (for cases where the age doesn't exist)
      dplyr::mutate(CC_Death_by_diff_per_t = coalesce(CC_Death_by_diff, 0)) %>%
      
      # Compute total_alive
      dplyr::mutate(total_alive = H + HR.HPV.infection + CIN1 + CIN2 + CIN3 +
                      FIGO.I + FIGO.II + FIGO.III + FIGO.IV + Survival) %>%
      
      # Compute CC_mortality based on CC_Death_per_t and total_alive
      dplyr::mutate(CC_by_diff_mortality = (CC_Death_by_diff_per_t / total_alive) * 10^5) %>%
      
      # Compute age intervals and average CC_mortality by age intervals
      dplyr::mutate(age_interval = cut(age, breaks = breaks, labels = labels, right = FALSE)) %>%
      dplyr::group_by(age_interval) %>%
      dplyr::summarise(CC_by_diff_mean_mortality = mean(CC_by_diff_mortality, na.rm = TRUE)) %>%
      dplyr::ungroup()
    
    return(df)
    #sim_stalked_result[[1]]$CC_by_diff_mean_mortality <- df
    #return(sim_stalked_result)
  }
  ##############################################################################
  
  sim_result[[strategy_name]]$CC_by_diff_mean_mortality <-
    mean_CC_mortality_by_diff_func(sim_stalked_result = sim_result[[strategy_name]],
                                   my_Probs = my_Probs)  
  
  
  ##############################################################################
  # B. Cancer-unrelated Mortality
  other_mean_mortality_func <- function(sim_stalked_result, my_Probs) {
    
    age_intervals <- my_Probs %>% 
      dplyr::select(Lower, Larger) %>% 
      unique() %>% 
      arrange(Lower)
    
    # Create a vector of the breaks for the intervals
    breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
    
    # Create labels for the intervals
    labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
    
    # Compute prevalence and average it by age intervals
    df <- sim_stalked_result$TR %>% 
      #dplyr::select(sim, cycle, age, H, HR.HPV.infection) %>% 
      dplyr::select(everything()) %>% 
      dplyr::mutate(total_alive = H + HR.HPV.infection + CIN1 + CIN2 + CIN3 +
                      FIGO.I + FIGO.II + FIGO.III + FIGO.IV + Survival) %>%
      dplyr::mutate(other_mortality = (Other.Death / total_alive) * 10^5) %>% 
      dplyr::mutate(age_interval = cut(age, breaks = breaks, labels = labels, right = FALSE)) %>% 
      dplyr::group_by(age_interval) %>% 
      dplyr::summarise(other_mean_mortality = mean(other_mortality, na.rm = TRUE)) %>% 
      dplyr::ungroup()
    
    return(df)
    #sim_stalked_result[[1]]$other_mean_mortality <- df
    #return(sim_stalked_result)
  }
  ##############################################################################
  
  sim_result[[strategy_name]]$other_mean_mortality <-
    other_mean_mortality_func(sim_stalked_result = sim_result[[strategy_name]],
                              my_Probs = my_Probs)  
  
  
  ##############################################################################
  # Mean FIGO states across simulations by age interval
  mean_FIGO_prevalence_Func <- function(sim_stalked_result, my_Probs) {
    age_intervals <- my_Probs %>% 
      dplyr::select(Lower, Larger) %>% 
      unique() %>% 
      arrange(Lower)
    
    # Create a vector of the breaks for the intervals
    breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
    
    # Create labels for the intervals
    labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
    
    # Compute prevalence and average it by age intervals
    df <- sim_stalked_result$TR %>% 
      dplyr::select(everything()) %>% 
      dplyr::mutate(total_alive = H + HR.HPV.infection + CIN1 + CIN2 + CIN3 +
                      FIGO.I + FIGO.II + FIGO.III + FIGO.IV + Survival) %>%
      # Prevalence for each FIGO state
      dplyr::mutate(FIGO.I_prev = (FIGO.I / total_alive) * 10^5) %>% 
      dplyr::mutate(FIGO.II_prev = (FIGO.II / total_alive) * 10^5) %>% 
      dplyr::mutate(FIGO.III_prev = (FIGO.III / total_alive) * 10^5) %>% 
      dplyr::mutate(FIGO.IV_prev = (FIGO.IV / total_alive) * 10^5) %>% 
      # Assigning age intervals
      dplyr::mutate(age_interval = cut(age, breaks = breaks, 
                                       labels = labels, right = FALSE)) %>% 
      dplyr::group_by(age_interval) %>% 
      # Summarizing the mean prevalence for each FIGO state
      dplyr::summarise(mean_FIGO.I_prev = mean(FIGO.I_prev, na.rm = TRUE),
                       mean_FIGO.II_prev = mean(FIGO.II_prev, na.rm = TRUE),
                       mean_FIGO.III_prev = mean(FIGO.III_prev, na.rm = TRUE),
                       mean_FIGO.IV_prev = mean(FIGO.IV_prev, na.rm = TRUE)) %>% 
      dplyr::ungroup()
    
    # # Storing the results in the simulation object
    # sim_stalked_result[[1]]$mean_FIGO_prevalence <- df
    # return(sim_stalked_result) 
    return(df)
  }
  ##############################################################################
  
  
  sim_result[[strategy_name]]$mean_FIGO_prevalence <-
    mean_FIGO_prevalence_Func(sim_stalked_result = 
                                sim_result[[strategy_name]], my_Probs = my_Probs)  
  ##############################################################################
  
  ##############################################################################
  # Mean (accross simulations) of Cancer (FIGO.I-.IV)
  mean_Figo_Func  <- function (sim_stalked_result, my_Probs) {
    age_intervals <- my_Probs %>% 
      dplyr::select(Lower, Larger) %>% 
      unique() %>% 
      arrange(Lower)
    
    # Create a vector of the breaks for the intervals
    breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
    
    # Create labels for the intervals
    labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
    
    # Compute prevalence and average it by age intervals
    df <- sim_stalked_result$TR %>% 
      dplyr::select(everything()) %>% 
      #dplyr::mutate(total_alive = H + HR.HPV.infection + CIN1 + CIN2 + CIN3 +
      #                FIGO.I + FIGO.II + FIGO.III + FIGO.IV + Survival) %>%
      # Prevalence for each FIGO state
      dplyr::mutate(FIGO.I   = (FIGO.I)) %>%  
      dplyr::mutate(FIGO.II  = (FIGO.II)) %>% 
      dplyr::mutate(FIGO.III = (FIGO.III)) %>% 
      dplyr::mutate(FIGO.IV  = (FIGO.IV)) %>% 
      # Assigning age intervals
      dplyr::mutate(age_interval = 
                      cut(age, breaks = 
                            breaks, labels = labels, right = FALSE)) %>% 
      dplyr::group_by(age_interval) %>% 
      # Summarizing the mean prevalence for each FIGO state
      dplyr::summarise(mean_FIGO.I = mean(FIGO.I, na.rm = TRUE),
                       mean_FIGO.II = mean(FIGO.II, na.rm = TRUE),
                       mean_FIGO.III = mean(FIGO.III, na.rm = TRUE),
                       mean_FIGO.IV = mean(FIGO.IV, na.rm = TRUE)) %>% 
      dplyr::ungroup()
    
    ## Storing the results in the simulation object
    #sim_stalked_result[[1]]$mean_FIGO <- df
    #return(sim_stalked_result) 
    return(df)
  }
  ##############################################################################
  
  # Concatenate the prevalence to the sim result 
  sim_result[[strategy_name]]$mean_FIGO <-
    mean_Figo_Func(sim_stalked_result = 
                     sim_result[[strategy_name]], my_Probs = my_Probs)  
  ##############################################################################
  
  ##############################################################################
  # Mean diagnosed of Cancer averaged by age intervals (FIGO.I-.IV) and by sims
  mean_Diagnosed_Per_Symp_Func  <- function (sim_stacked_result, my_Probs) {
    age_intervals <- my_Probs %>% 
      dplyr::select(Lower, Larger) %>% 
      unique() %>% 
      arrange(Lower)
    
    # Create a vector of the breaks for the intervals
    breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
    
    # Create labels for the intervals
    labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
    
    # extracting diagnosed:
    sympt <- sim_stacked_result$symptomatics
    
    # add age column:
    sympt <- sympt %>% dplyr::mutate(age = TimeStep + 9)
    
    # add age interval column:
    sympt <- sympt %>%
      mutate(age_interval = cut(age, 
                                breaks = seq(min(breaks), max(breaks), by = 5), 
                                labels = labels, 
                                right = FALSE))
    
    # Calculate the maximum simulation count
    max_sim <- max(sympt$sim)
    
    # Create a complete data frame with all combinations of age intervals and DiagnosedStates
    complete_data <- expand.grid(
      age_interval = labels,
      DiagnosedState = c("FIGO.I", "FIGO.II", "FIGO.III", "FIGO.IV")
    )
    
    # Summarize data in the desired format
    df <- sympt %>%
      group_by(age_interval, DiagnosedState) %>%
      summarise(mean_count = n() / max_sim, .groups = "drop") %>%
      right_join(complete_data, by = c("age_interval", "DiagnosedState")) %>%
      replace_na(list(mean_count = 0)) %>%  # Replace NA values with 0
      pivot_wider(names_from = DiagnosedState, 
                  values_from = mean_count, 
                  names_prefix = "mean_Diagnosed_") %>%
      arrange(age_interval)
    
    ## View the result
    #print(df)
    
    ## Storing the results in the simulation object
    #sim_stacked_result[[1]]$mean_Diagnosed <- df 
    #return(sim_stacked_result) 
    return(df)
  }
  ##############################################################################
  
  sim_result[[strategy_name]]$mean_Diagnosed_per_Symp  <-
    mean_Diagnosed_Per_Symp_Func(sim_stacked_result =
                          sim_result[[strategy_name]], my_Probs = my_Probs)  
  ##############################################################################

  ##############################################################################
  ## Compute averaged new_cases per age interval:
  ##############################################################################
  mean_new_cases_func <- function(sim_result, new_cases_name, my_Probs) {
    # Extract age intervals from my_Probs
    age_intervals <- my_Probs %>%
      dplyr::select(Lower, Larger) %>%
      unique() %>%
      dplyr::arrange(Lower)
    
    # Define breaks and labels
    breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
    labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
    
    # Get the new cases data.frame by name
    #df <- sim_result[[1]][[new_cases_name]] %>%
    df <- sim_result[[new_cases_name]] %>% dplyr::select(-row_names) %>%
      dplyr::rename(new_cases = 2) %>%  # assumes second column is the one with cases
      dplyr::mutate(age_interval = cut(age, breaks = breaks, labels = labels, right = FALSE)) %>%
      dplyr::group_by(sim, age_interval) %>%
      dplyr::summarise(mean_new_cases = mean(new_cases, na.rm = TRUE), .groups = "drop") %>%
      dplyr::group_by(age_interval) %>%
      dplyr::summarise(mean_new_cases = mean(mean_new_cases, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(mean_new_cases = round(mean_new_cases, 2))  # keep up to two decimals
    
    # Return summary table
    return(df)
  }
  ##############################################################################
  new_averaged_CIN1_per_age_interval <- 
    mean_new_cases_func(sim_result = sim_result[[strategy_name]], 
                        new_cases_name = "new_CIN1",
                        my_Probs = my_Probs)
  
  new_averaged_CIN2_per_age_interval <- 
    mean_new_cases_func(sim_result = sim_result[[strategy_name]], 
                        new_cases_name = "new_CIN2",
                        my_Probs = my_Probs)
  
  new_averaged_CIN3_per_age_interval <- 
    mean_new_cases_func(sim_result = sim_result[[strategy_name]], 
                        new_cases_name = "new_CIN3",
                        my_Probs = my_Probs)
  
  new_averaged_Cancer_per_age_interval <- 
    mean_new_cases_func(sim_result = sim_result[[strategy_name]], 
                        new_cases_name = "new_Cancer",
                        my_Probs = my_Probs)
  
  new_averaged_CC_Death_per_age_interval <- 
    mean_new_cases_func(sim_result = sim_result[[strategy_name]], 
                        new_cases_name = "new_CC_Death",
                        my_Probs = my_Probs)
  
  sim_result[[strategy_name]]$new_averaged_CIN1_per_age_interval <- 
    new_averaged_CIN1_per_age_interval
  
  sim_result[[strategy_name]]$new_averaged_CIN2_per_age_interval <- 
    new_averaged_CIN2_per_age_interval
  
  sim_result[[strategy_name]]$new_averaged_CIN3_per_age_interval <- 
    new_averaged_CIN3_per_age_interval
  
  sim_result[[strategy_name]]$new_averaged_Cancer_per_age_interval <- 
    new_averaged_Cancer_per_age_interval
  
  sim_result[[strategy_name]]$new_averaged_CC_Death_per_age_interval <- 
    new_averaged_CC_Death_per_age_interval
  
  ##############################################################################
 
  ##### Adding mean diagnosed per symptoms FIGOs and total cancer:
  sim_result[[strategy_name]]$mean_Diagnosed_per_symp_FIGOI <-
    sim_result[[strategy_name]]$symptomatics %>%
    dplyr::filter(DiagnosedState == 'FIGO.I') %>%
    nrow() %>%
    `/`(numb_of_sims)
  
  sim_result[[strategy_name]]$mean_Diagnosed_Per_Symp_FIGOII <-
    sim_result[[strategy_name]]$symptomatics %>%
    dplyr::filter(DiagnosedState == 'FIGO.II') %>%
    nrow() %>%
    `/`(numb_of_sims)
  
  sim_result[[strategy_name]]$mean_Diagnosed_Per_Symp_FIGOIII <-
    sim_result[[strategy_name]]$symptomatics %>%
    dplyr::filter(DiagnosedState == 'FIGO.III') %>%
    nrow() %>%
    `/`(numb_of_sims)
  
  sim_result[[strategy_name]]$mean_Diagnosed_Per_Symp_FIGOIV <-
    sim_result[[strategy_name]]$symptomatics %>%
    dplyr::filter(DiagnosedState == 'FIGO.IV') %>%
    nrow() %>%
    `/`(numb_of_sims)
  
  sim_result[[strategy_name]]$mean_Diagnosed_Per_Symp_Cancer <-
    sim_result[[strategy_name]]$symptomatics %>%
    nrow() %>%
    `/`(numb_of_sims)
  
  ##############################################################################
  ##############################################################################
  
  ## Cleaning
  source("./R/Remove_columnS_from_list_Func.R")
  sim_result[[strategy_name]] <- 
    remove_columns_from_list(complex_list = sim_result[[strategy_name]], 
                             ... = "sim.1", "row_names")
  
  
  ##############################################################################
  
} # endfor strategy_name

## Further Cleaning
remove_symptomatics <- function(sim_result, drop = TRUE) {
  if (drop) {
    sim_result <- lapply(sim_result, function(strategy) {
      strategy[setdiff(names(strategy), "symptomatics")]
    })
  }
  sim_result
}


# removing symptomatics object:
sim_result <- remove_symptomatics(sim_result, drop = TRUE)

################################################################################
################################################################################
##### END POST-PROCESSING ######################################################
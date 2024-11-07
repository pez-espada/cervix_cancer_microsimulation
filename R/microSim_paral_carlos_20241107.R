
MicroSim <- function(strategy="natural_history", numb_of_sims = 30,
                     v_M_1, n_i, n_t, v_n, d_c, d_e, TR_out = TRUE, 
                     TS_out = TRUE, Trt = FALSE,  seed = 1, Pmatrix) 
{
  seeds <- sample(1:10000, numb_of_sims, replace = FALSE)  # Generate random seeds
  seeds <- sample(1:100000, numb_of_sims, replace = FALSE)  # Generate random seeds
  
  
  # Register the parallel backend
  cl <- makeCluster(n_cores, timeout = 6*60*60) # 6-hours timeout to prevent socket drop issues
  clusterExport(cl, c("Costs_per_Cancer_Diag", "Effs", "trans_prb", "Probs",
                      "my_Probs", "utilityCoefs", "v_n", "samplev", 
                      "my_age_prob_matrix_func","diagnose_column", 
                      "update_column", "states_to_check", "symptom_prob_vec",
                      "survival_prob_vec", "global_diagnosed", 
                      "cost_Vec", "new_cases_2"))
  registerDoParallel(cl)
  #registerDoSEQ()
  
  simulation_results <- list() 
  
  #for(sim in 1:numb_of_sims) {
  # Parallel processing using foreach
  simulation_results <- 
    foreach(sim = 1:numb_of_sims, .packages = c("dplyr", "tidyr", "purrr") ) %dopar% { 
      
      cat("Running simulation", sim, "with seed", seeds[sim], "\n")
      symptomatics <-
        data.frame(ID = integer(), TimeStep = integer(), 
                   DiagnosedState = character(), 
                   RecoveredFromState = logical(), stringsAsFactors = FALSE)
      
      v_dwc <- 1 / (1 + d_c) ^ (0:(n_t-1))   # calculate the cost discount weight based
                                             # on the discount rate d_c 
      v_dwe <- 1 / (1 + d_e) ^ (0:(n_t-1))   # calculate the QALY discount weight based 
                                             # on the discount rate d.e
      
      # Create the matrix capturing the state name/costs/health outcomes 
      # for all individuals at each time point:
      #m_M <- m_C <- m_E <-  matrix(nrow = n_i, ncol = (n_t + 1), 
      m_M <- m_C <- m_E <- 
        matrix(nrow = n_i, ncol = (n_t), 
               dimnames = list( 1:n_i, 
                                #paste0("cycle_", 1:(n_t + 1), sep = "")))  
                                paste0("cycle_", 1:(n_t), sep = "")))  
      
      m_M[, 1] <- v_M_1  # indicate the initial health state   
      
      seed <- seeds[sim]
      #seed <- 17
      set.seed(seed) # set the seed for every individual for the random number generator
      
      
      m_C[, 1] <- Costs_per_Cancer_Diag(M_it = m_M[, 1], # estimate costs per individual for the 
                                        symptomatics = symptomatics,
                                        time_iteration = 1,
                                        cost_Vec = cost_Vec, # initial health state
                                        Trt)             
      
      m_E[, 1] <- Effs(m_M[, 1], Trt, utilityCoefs = utilityCoefs) # estimate QALYs
                                                                   # per individual 
                                                                   # for the initial
                                                                   # health state  
      stored_list <- list()
      ###################### run over all the cycles ########################### 
      #for (t in 1:(n_t)) {
      for (t in 1:(n_t-1)) {
        ########################################################################
        # Select the transition matrix based on the cycle `n_t`:
        # Since our age intervals start at 10 years old,
        age_in_loop <- t + 9
        ########################################################################
        
        # update/correct n_s (<<- let change variable from inside a function):
        n_s  <<- length(v_n)  
        
        ######################################################################## 
        #new code:
        new_entries <- diagnose_column(m_M[, t], t)
        
        if (!is.null(new_entries)) {
          stored_list[[t]] <- new_entries
        }
        if (nrow(new_entries) > 0) {
          symptomatics <- bind_rows(symptomatics, new_entries)
        }
        ######################################################################## 
        
        ########################################################################
        # NOTE: if my_age_in_loop = age_in_loop (without adding 1), then the 
        # microsim does not agree with the markov (for whatever reason)
        my_age_prob_matrix <- 
          my_age_prob_matrix_func(my_Prob_matrix = my_Probs, 
                                  my_age_in_loop = (age_in_loop + 1))
        # Add colnames and update `v_n`:
        rownames(my_age_prob_matrix) <- v_n <<- 
          my_age_prob_matrix %>%
          dplyr::select(-c(Age.group, Lower, Larger)) %>% 
          colnames()
        
        # Extract the transition probabilities of each individuals at cycle t
        # given the individual current state and the corresponding 
        # transition probability matrix that depends on age:
        # Next time (t+1) transition
        m_P <- Probs(M_it =  m_M[, t], my_Probs = my_age_prob_matrix)
        
        m_M[, t + 1] <- samplev(probs = m_P, m = 1)  # sample the next health state 
        # and store that state in  
        # matrix m_M 
        ########################################################################    
        
        # m_M[, t + 1] <- update_column(m_M[, t], new_entries)
        next_col <- m_M[, t + 1]
        next_col <- update_column(m_M[, t], new_entries, next_col)
        
        # Ensure next_col updates are preserved after sampling
        m_M[, t + 1] <- ifelse(next_col == "Survival", "Survival", m_M[, t + 1])
        
        ########################################################################    
        ## Costs per CC diagnose at time t + 1.
        # Estimate costs per individual during cycle t + 1 conditional on treatment:
        # Debugging:
        #m_C[, t] <-                              
        m_C[, t + 1] <-                              
          Costs_per_Cancer_Diag(M_it = m_M[, t + 1],  
                                symptomatics = symptomatics,
                                time_iteration = t,
                                cost_Vec = cost_Vec,    
                                Trt)            
      
        m_E[, t + 1] <- # estimate QALYs per individual during cycle t + 1
          Effs( m_M[, t + 1], Trt, 
                utilityCoefs = utilityCoefs)                   
        ############################################################################    
        cat('\r', paste(round(t/n_t * 100),          # display the 
                        "% done\n", sep = " "))        # progress of  the simulation                    
        
      }  
      ######################## close loop for cycles ############################### 
      
      # Combine stored entries in a single data frame()
      symptomatics <- bind_rows(stored_list)
      tc_disc <- m_C[,1:n_t] %*% v_dwc       # total (discounted) cost per individual
      te_disc <- m_E[,1:n_t] %*% v_dwe       # total (discounted) QALYs per individual 
      
      tc_undisc <- m_C[,1:n_t] %*% rep(1, n_t)       # total (discounted) cost per individual
      te_undisc <- m_E[,1:n_t] %*% rep(1, n_t)       # total (discounted) QALYs per individual 
      
      tc_hat_disc <- mean(tc_disc)        # average (discounted) cost 
      te_hat_disc <- mean(te_disc)        # average (discounted) QALYs
      tc_hat_undisc <- mean(tc_undisc)    # average (discounted) cost 
      te_hat_undisc <- mean(te_undisc)    # average (discounted) QALYs
      
      # Create a matrix of transitions across states transitions from one state to the other:
      if (TS_out == TRUE) {  
        TS <- paste(m_M, cbind(m_M[, -1], NA), sep = "->")    
        
        TS <- matrix(TS, nrow = n_i)
        rownames(TS) <- paste("Ind",   1:n_i, sep = " ")   # name the rows 
        #colnames(TS) <- paste0("cycle_", 1:(n_t + 1), sep = "")   # name the columns 
        colnames(TS) <- paste0("cycle_", 1:(n_t), sep = "")   # name the columns 
      } else {
        TS <- NULL
      }
      
      if (TR_out == TRUE) {
        TR <- t(apply(m_M, 2, 
                      function(x) table(factor(x, levels = v_n, ordered = TRUE))))
        #TR <- TR / n_i                                   # create a distribution 
        # trace
        
        #rownames(TR) <- paste("cycle", 1:(n_t + 1), sep = "_") # name the rows 
        rownames(TR) <- paste("cycle", 1:(n_t), sep = "_") # name the rows 
        colnames(TR) <- v_n                              # name the columns 
      } else {
        TR <- NULL
      }
      
      # If TS_out == TRUE we can then compute the number of new cases for each type
      # of cancer state per time (cycle). A new case of cancer state X in time t
      # is defined as an individual transition to this state X provided the
      # individual was not in that state X a time t-1.
      # NOTE that the TR output display individual transitions at each cycle t
      # that are going to occur at t + 1. That is, "XX->YY" in cycle t meant that the
      # corresponding individual is in state "XX" in t and is transiting to state
      # "YY" in t + 1.
      # A character with all transitions:
      transitions <- 
        TS %>% 
        as_tibble() %>% 
        pivot_longer(everything(), names_to = "column") %>% 
        distinct(value) %>%
        unique() %>% 
        as.list() %>%
        unlist()
      
      if(TS_out == TRUE){
        Tot_Trans_per_t <- 
          t(apply(TS, 2, 
                  function(x) 
                    table(factor(x, levels 
                                 = transitions, 
                                 ordered = TRUE))))
        # trace
        #rownames(Tot_Trans_per_t) <- paste0("cycle_", 1:(n_t + 1), sep = "") # name the rows 
        rownames(Tot_Trans_per_t) <- paste0("cycle_", 1:(n_t), sep = "") # name the rows 
      } else {
        Tot_Trans_per_t <- NULL
      }
      
      # New cases:
      new_CIN1 <- new_cases_2(state1 = "HR.HPV.infection", state2 = "CIN1", 
                              Tot_Trans_per_t = Tot_Trans_per_t)
      
      new_CIN2 <- new_cases_2(state1 = "CIN1", state2 = "CIN2", 
                              Tot_Trans_per_t = Tot_Trans_per_t)
      
      new_CIN3 <- new_cases_2(state1 = "CIN2", state2 = "CIN3", 
                              Tot_Trans_per_t = Tot_Trans_per_t)
      
      new_Cancer <- new_cases_2(state1 = "CIN3", state2 = "FIGO.I", 
                                Tot_Trans_per_t = Tot_Trans_per_t)
      
      new_CC_Death <- new_cases_2(state1 = c("CIN1", "CIN2","CIN3","FIGO.I", 
                                             "FIGO.II", "FIGO.III", "FIGO.IV"),
                                  state2 = "CC_Death", 
                                  Tot_Trans_per_t = Tot_Trans_per_t)
      
      
      # Before sending back, some cleaning regarding cycle `n_t+1` which is 
      # computed but no needed as a result:
      m_M <-m_M[, 1:n_t]
      m_C <-m_C[, 1:n_t]
      m_E <-m_E[, 1:n_t]
      new_CIN1 <- new_CIN1 %>% dplyr::slice(c(1:n_t))
      new_CIN2 <- new_CIN2 %>% dplyr::slice(c(1:n_t))
      new_CIN3 <- new_CIN3 %>% dplyr::slice(c(1:n_t))
      new_Cancer <- new_Cancer %>% dplyr::slice(c(1:n_t))
      new_CC_Death <- new_CC_Death %>% 
        dplyr::select(CC_Death_per_t, age, cycle) %>%
        dplyr::slice(c(1:n_t))
      
      # Removing no needed extra row from TR:
      row_to_remove <- n_t + 1
      TR <- TR[-row_to_remove, ]
      
      rm(row_to_remove)
      
      # Removing extra column no needed in TS
      TS <- TS[, -(n_t + 1)]
      
      ### add age to TR:
      TR <- as.data.frame(TR)
      #TR <- TR %>% mutate(age = row_number() + 10)
      #TR <- TR %>% mutate(age = row_number() + 8)
      TR <- TR %>% mutate(age = row_number() + 9)
      TR$sim <- sim
      
      #Remove large objects: 
      rm(m_M, m_C, m_E)
      
      # Computing new cancer cases pert cycle using diff() function:
      CC_Death_by_diff <- c(0, TR %>% 
                              select(CC_Death) %>% 
                              as_vector() %>% 
                              diff())
      TR$CC_Death_by_diff <- CC_Death_by_diff 
      TR$CC_Death_by_diff <- ifelse( TR$age==10, 0, TR$CC_Death_by_diff)
      
      CC_Death_by_diff <- TR %>% 
        dplyr::select(sim, age, CC_Death_by_diff) %>% 
        dplyr::as_tibble()
      
      
      # Store the results from the simulation in a list
      results <- list(strategy = strategy,
                      #seed = seeds[sim],
                      seed = seed,
                      sim_numb = sim, 
                      #m_M = m_M, 
                      #m_C = m_C, 
                      #m_E = m_E, 
                      #tc_disc = tc_disc, 
                      #tc_undisc = tc_undisc,
                      #te_disc = te_disc,
                      #te_undisc = te_undisc,
                      tc_hat_disc = tc_hat_disc,
                      tc_hat_undisc = tc_hat_undisc,
                      te_hat_disc = te_hat_undisc, 
                      te_hat_undisc = te_hat_undisc, 
                      #TS = TS,
                      TR = TR, 
                      #Tot_Trans_per_t = Tot_Trans_per_t, 
                      #symptomatics = symptomatics,
                      new_CIN1 = new_CIN1,
                      new_CIN2 = new_CIN2,
                      new_CIN3 = new_CIN3,
                      new_Cancer = new_Cancer,
                      new_CC_Death = new_CC_Death,
                      CC_Death_by_diff = CC_Death_by_diff)  
      
      #results$seed <- seeds[sim]
      #simulation_results[sim] <- list(results)
      cat("At sim number:", sim,  " tc_hat_undisc is ", tc_hat_undisc, "\n")
      rm(symptomatics)
      
      return(results)
     #gc() #Force memory cleanup after each sim/batch 
      
    } # end of `foreach/dopar` loop
  
  
  #stopCluster(cl)  # Stop the cluster when done
  
  #return(simulation_results)
  
  # stack results
  #source("./R/Sumarize_results_by_Strategy_Func.R")
  stacked_results <- 
    summarize_results_by_Strategy(results_list = simulation_results, 
                                  numb_of_sims = numb_of_sims)
  
  stopCluster(cl)  # Stop the cluster when done
  return(stacked_results)
} # end of MicroSim function

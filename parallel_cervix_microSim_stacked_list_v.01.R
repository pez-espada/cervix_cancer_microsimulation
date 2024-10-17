knitr::opts_chunk$set(echo = TRUE)


################################################################################
# This code is a modified version of the original code from:
# [https://github.com/DARTH-git/Microsimulation-tutorial] (Krijkamp et al 2018 
# Sick-Sicker model).
# Modifications by: Carlos Dommar D'Lima - carlos.dommar@gmail.com
# This code extends the "sick-sicker" model of the original authors to a
# multi-state cervix cancer model
# Parallelized version
################################################################################
rm(list = ls())
library(tidyverse)

# Sources:
# Sandra's function:
source("./R/sumarize_results_by_Strategy_Func.R")

# Ensure necessary libraries are loaded
ensure_library <- function(...) {
  pkgs <- unlist(list(...))
  pkgs <- gsub("[\"']", "", pkgs) # Remove quotes
  sapply(pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  })
}
################################################################################

my_Probs <- readRDS(file = "./data/probs.rds")

my_Probs <- # transition matrix (for all sim cycles) 
  my_Probs %>%
  as_tibble() # I need a tibble to use 'rename' function down there:

# tidying up a bit the transition matrix:
my_Probs <- my_Probs %>% dplyr::rename("H" = "Well")

my_Probs <- my_Probs %>% as.data.frame() # convert back to data.frame (no needed?)

###############################################################
# Function to extract and convert numbers from factor levels
extract_numbers <- function(range_factor) {
  range_string <- as.character(range_factor)
  numbers <- as.numeric(unlist(strsplit(range_string, "-")))
  return(numbers)
}
###############################################################

# Apply the function to the Range column and create new columns
my_Probs$Lower  <- sapply(my_Probs$Age.group, function(x) extract_numbers(x)[1])
my_Probs$Larger <- sapply(my_Probs$Age.group, function(x) extract_numbers(x)[2])
# For the last cycle/iteration we need to adjust the last transition matrix:
my_Probs$Larger <- 
  ifelse(my_Probs$Larger == max(my_Probs$Larger), my_Probs$Larger + 1, my_Probs$Larger) 


## ----model parameters-------------------------------------------------------------------------------------------------------------------------------------------------------------------
n_i <- 10^4                 # number of simulated individuals
n_t <- 75                   # time horizon, 75 cycles (it starts from 1)

################################################################################
### (THIS IS WORK IN PROGRESS):
# cycle_period can go from one month to one year. that is
# I think a sensible way is to offer the following frequencies
cycle_period <- "1mth"
cycle_period <- "2mth"
cycle_period <- "3mth"
cycle_period <- "4mth"
cycle_period <- "6mth"
cycle_period <- "1yr" # i.e. 12mth
if (cycle_period == "1mth"){
  n_t <- n_t * 12
} else if (cycle_period == "6mth") {
  n_t <- n_t * 2
} else if (cycle_period == "1yr") {
  n_t <- n_t * 1
}
################################################################################

v_n <- rownames(my_Probs)
v_n <- colnames(my_Probs)
v_n <- v_n[-c(1,14,15)]
n_s   <- length(v_n)                # the number of health states
v_M_1 <- rep("H", n_i)              # everyone begins in the healthy state 
#v_M_1 <- rep("Well", n_i)           # everyone begins in the healthy state 
d_c   <- d_e <- 0.03                # equal discounting of costs and QALYs by 3%
v_Trt <-
  c("No Treatment", "Treatment")    # store the strategy names
################################################################################

# Cost and utility inputs 
# From our Markov cervix model (CC's natural history?):
cost_Vec = c(0, 39.54, 288.91, 1552.27, 1552.27, 
             5759.81, 12903.63, 23032.41, 35323.14, 0, 0, 0)
utilityCoefs = c(1, 1, 0.987, 0.87, 0.87, 0.76, 0.67, 0.67, 0.67, 0.938, 0, 0)


################################################################################
#### For extracting the probabilities of transitions given the transition matrix:
########### Probably the following function is not needed ######################
#' Extract transition probability from Transition Matrix
#'
#' @param P 
#' @param state1 
#' @param state2 
#'
#' @return a numeric scalar corresponding to the asked probability of transition
#' @export
#'
#' @examples
#' trans_prb(P = my_Probs, state1 = "Well", state2 = "HR.HPV.infection") 
#' trans_prb(P = my_Probs, state1 = "CIN1", state2 = "CIN2") 
trans_prb <- function(P, state1, state2) {
  # If the matrix of transition, P, is given:
  # the probability of an individual to go to state 'state2' the next time
  # step given the individual is currently in state 'state1' is computed by:
  tryCatch(
    transition_prob <- P %>%  
      filter(row.names(P) %in% c(state1)) %>% # filter state1 row
      dplyr::select(all_of(state2)) %>%   # select state2 column
      as.numeric(),
    error = function(e){
      message("An error occurred:\n", e)
      print("Remember the valid states are:")
      P %>% rownames() %>% print()
    },
    warning = function(w){
      message("A warning occured:\n", w)
    }
  )
  return(transition_prob)
}
################################################################################



################################################################################
############################# Sampling function ################################
# Efficient implementation of the rMultinom() function of the Hmisc package #### 
samplev <- function (probs, m) {
  d <- dim(probs) # i.e. number of individuals times number of states: n_i x n_s
  n <- d[1]       # number of individuals n_s
  k <- d[2]       # number of states
  lev <- dimnames(probs)[[2]] # vector with  names of health states
  if (!length(lev)) # checks if `lev` vector (states names) is empty 
                    # or has length 0 
    lev <- 1:k # if empty (evaluates to `TRUE`), it assigns numeric state labels
               # (1:k) to `lev`
  ran <- 
    matrix(lev[1], ncol = m, nrow = n) # create array n_s x m (m=1) 
                                       # consisting in of health-state stored in
                                       # `lev[1]`, "H" in our case.
  
  ########## Creating the matrix of cumulative distributions U #################
  U <- t(probs)    # transpose probs from (`n_i*n_s`) to (`n_s*n_i`)
  for(i in 2:k) {  
    # This loop fills U with the cumulative probabilities of each individual
    # across all its possible transitions (`v_s`or `lev` within thus function).
    # That is each column of `U` represents the cumulative distribution for each
    # individual across its corresponding transitions. 
    # The last element of each column must sum 1 (or close enough:)
    U[i, ] <- U[i, ] + U[i - 1, ]
  }
  if (any((U[k, ] - 1) > 1e-04))
    stop("error in multinom: probabilities do not sum to 1")
  ##############################################################################
  
  ### Random sampling, binning, and moving states: 
  for (j in 1:m) {
    un <- rep(runif(n), rep(k, n)) # repeat `runif(n)` `rep(k,n)`times
    # this create a numeric of `n_i x n_s` that 
    # sample  an uniformed distributed number 
    # between 0 and 1. The generated random number
    # repeats itself `n_s` times and then another 
    # rand unif number is drawn. This process is 
    # carried out `n_i` times. NOTE: every time
    # runif() is run it produces a new random sample
    # i.e. it does not seem dependent on the seed
    
    # Here's where we choose the individuals' next states:
    ran[, j] <- lev[1 + colSums(un > U)]
  }
  ran
}
################################################################################

################################################################################
######################### Probability function #################################
## The Probs function that updates the transition probabilities of every cycle:
Probs <- function(M_it, my_Probs) {
  n_s <- length(v_n)
  n_i <- length(M_it)
  m_P_it <- matrix(NA, n_s, n_i) 
  rownames(m_P_it) <- v_n
  for (i in 1:length(v_n)) {
    state_mask <- !is.na(M_it) & M_it == v_n[i]
    
    if (sum(state_mask) > 0) {
      m_P_it[, state_mask] <- 
        lapply(X = v_n, function(x) trans_prb(P = my_Probs, state1 =
                                                v_n[i], state2 = x)) %>%
        unlist()
    } else {
      ## Debugging:
      #cat("State", v_n[i], "is not present in M_it at this time step\n")
    }
  }
  if (any(is.na(m_P_it))) {
    # Diagnostic message
    cat("Transition probabilities contain NA values\n")
  }
  ifelse(colSums(m_P_it, na.rm = TRUE) >= .991, 
         return(t(m_P_it)), 
         stop("Probabilities do not sum to 1"))
}
################################################################################


################################################################################
### Costs function
# The `Costs_per_Cancer_Diag` function estimates the costs of a diagnose 
# individual due to cancer symptoms (FIGO.I-IV) at every cycle. 
# This cost is only charged once in the patient's lifetime.
# NOTE: need to decide if the cost is applied on current time `t` or `t+1` as it is now.
Costs_per_Cancer_Diag <- function(M_it, cost_Vec, symptomatics,
                                  time_iteration, Trt = FALSE) {
  c_it <- rep(0, length(M_it))
  ci_t <- 0
  if(nrow(symptomatics) > 0 ) {
    c_it[symptomatics %>% 
           filter(DiagnosedState == "FIGO.I" & TimeStep == time_iteration) %>% 
           select(ID) %>% as.list() %>% 
           unlist()] <- cost_Vec[which(v_n %in% "FIGO.I")]
    c_it[symptomatics %>% 
           filter(DiagnosedState == "FIGO.II" & TimeStep == time_iteration) %>% 
           select(ID) %>% as.list() %>% 
           unlist()] <- cost_Vec[which(v_n %in% "FIGO.II")]
    c_it[symptomatics %>%
           filter(DiagnosedState == "FIGO.III" & TimeStep == time_iteration) %>% 
           select(ID) %>% as.list() %>% 
           unlist()] <- cost_Vec[which(v_n %in% "FIGO.III")]
    c_it[symptomatics %>% 
           filter(DiagnosedState == "FIGO.IV" & TimeStep == time_iteration) %>% 
           select(ID) %>% as.list() %>% 
           unlist()] <- cost_Vec[which(v_n %in% "FIGO.IV")]
  }
  return(c_it)              		                           # return the costs
}
################################################################################


################################################################################
### Health outcome function QALYs
Effs <- function (M_it, Trt = FALSE, cl = 1, utilityCoefs) {
  # check length of vector of states and vector of utility/QALYs are the same:
  u_it <- 0                   # by default the utility for everyone is zero
  tryCatch(
    for (i in 1:length(utilityCoefs)) {
      u_it[M_it == v_n[i]] <- utilityCoefs[i]   # update the utility if healthy
    },
    error = function(e){
      message("An error occurred:\n", e)
      print("Check state vector and utility vector have the same dimensions:")
      P %>% rownames() %>% print()
    },
    warning = function(w){
      message("A warning occured:\n", w)
    }
  )
  return(u_it)
}
################################################################################


################################################################################
## ----time period related functions, echo=TRUE
########### WORK IN PROGRESS #########################
age_factor <- function(my_period) {
  # it receives a string with the period of the cycle, and it can be:
  #  - "1mth"
  #  - "3mth"
  #  - "4mth"
  #  - "6mth"
  #  - "1yr" # i.e. 12 months
  # and it gives back an age factor for scaling cycle period.
  if (my_period == "1yr") {
    my_factor <- 1
  } else if (my_period == "6mth") {
    my_factor <- 2
  } else if (my_period == "4mth") {
    my_factor <- 3
  } else if (my_period == "3mth") {
    my_factor <- 4
  } else if (my_period == "1mth") {
    my_factor <- 12
  } else {print("Cycle period can only be: '1yr', '6mth','4mth', '3mth' and '12mth'")}
  return(my_factor)
}
################################################################################
 

################################################################################
#### ! NOT USED ! ############################
convert_matrix_to_proper_transition <- 
  function(my_age_prob_matrix, cycle_period) {
    my_age_prob_matrix %>% head(3)
    ensure_library(c("expm", "pracma", "ctmcd"))
    trans_matrix <- my_age_prob_matrix %>% 
      select(-c("Age.group", "Lower", "Larger")) %>% 
      as.matrix()
    # Referenece: https://rpubs.com/crossxwill/transition_matrix
    ## method 1: (not working atm)
    #ensure_library(expm)
    #TM.exp  <- expm::expm((1 / age_factor(cycle_period))) * log(trans_matrix) 
    
    #method 2 ;
    #ensure_library("pracma")
    TM_pracma <- 
      pracma::rootm(trans_matrix, p=age_factor(cycle_period), 
                    kmax = 20, tol = 1e-10)
    round(TM_pracma$B, 5)
    # Regularization with the `ctmcd` package, The code below uses the 
    # quasi-optimization of the generator (QOG) approach from 
    # Kreinin and Sidelnikova (2001).:
    ensure_library("ctmcd")
    TM_qo <- ctmcd::gm(TM_pracma$B, te=1, method = "QO") 
  }
#### ! NOT USED ! ############################
################################################################################

## symptoms
# An individual can be in cancer states, i.e. FIGO.I, FIGO.II. FIGO.III and FIGO.IV
# (in the model) and yet no develop symptoms. Form th Markov cohort model we have
# that the probability of developing symptoms are 0.11, 0.23, 0.66, and 0.9 for
# FIGO1...4 respectively. Symptoms are important for the cost-effectiveness analysis
# I build four n_i x (n_t + 1) matrices each with the actual individuals who developed
# symptoms according the aforementioned probabilities.

ensure_library("dplyr", "tidyverse", "purrr")
# Function to process each column version 3:
figoSymProb <- c(0.11, 0.23, 0.66, 0.9) 
screeProbs <- c(0, 0, 1, 1, 1, 0.9688, 0.9066, 0.7064, 0.3986, 0, 0, 0)
symptom_prob_vec <- figoSymProb
survival_prob_vec <- screeProbs[6:9]
states_to_check <- c("FIGO.I", "FIGO.II", "FIGO.III", "FIGO.IV")

stored_list <- vector("list", n_t)

# Initialize a global vector to store all diagnosed individuals
global_diagnosed <- integer()

################################################################################
# Function receives a column with current state of `n_i`individuals and gives
# a dataframe with `ID, TimeStep`, `state`, and `RecoveredFromState` 
diagnose_column <- function(col, time_step) {
  new_entries <- data.frame(ID = integer(), 
                            TimeStep = integer(), DiagnosedState = character(),
                            RecoveredFromState = logical())
  for (state_idx in seq_along(states_to_check)) {
    state <- states_to_check[state_idx]
    prob_symptom <- symptom_prob_vec[state_idx]
    prob_survival <- survival_prob_vec[state_idx]
    in_state <- which(col == state)
    if (length(in_state) > 0) {
      # Remove individuals who have already been diagnosed
      in_state <- setdiff(in_state, global_diagnosed)
      if (length(in_state) > 0) {
        # Store based on diagnose probability
        to_store <- in_state[runif(length(in_state)) < prob_symptom]
        if (length(to_store) > 0) {
          # Add these individuals to the global diagnosed list
          global_diagnosed <<- c(global_diagnosed, to_store)
          # Check another probability to potentially change their state to "Survival"
          recovered <- to_store[runif(length(to_store)) < prob_survival]
          # Store the individuals' IDs, time steps, diagnosed states, and recovery status
          new_entries <- rbind(new_entries, data.frame(
            ID = to_store, 
            TimeStep = time_step, 
            DiagnosedState = state, 
            RecoveredFromState = to_store %in% recovered))
        }
      }
    }
  }
  rownames(new_entries) <- NULL
  return(new_entries)
}
################################################################################


#################################################################################
update_column <- function(col, new_entries, next_col) {
  if (nrow(new_entries) > 0) {
    diagnosed_ids <- new_entries$ID
    recovered_ids <- new_entries$ID[new_entries$RecoveredFromState]
    
    # Update the states in the next column for recovered individuals
    next_col[recovered_ids] <- "Survival"
    
    # Ensure that individuals who were diagnosed but not recovered retain their diagnosed state
    non_recovered_ids <- diagnosed_ids[!diagnosed_ids %in% recovered_ids]
    next_col[non_recovered_ids] <-
      new_entries$DiagnosedState[!diagnosed_ids %in% recovered_ids]
  }
  return(next_col)
}
#################################################################################


################################################################################
new_cases_2 <- function(state1, state2, Tot_Trans_per_t) {
  # Convert the data to a tibble for easier manipulation
  Tot_Trans_per_t_tbl <- as_tibble(Tot_Trans_per_t)
  
  # Case when state1 is a single string
  if (length(state1) == 1) {
    transition_column <- paste0(state1, "->", state2)  # Create the transition name
    
    # If the transition column exists, select it
    if (transition_column %in% colnames(Tot_Trans_per_t_tbl)) {
      transition_cases <- Tot_Trans_per_t_tbl %>%
        dplyr::select(all_of(transition_column)) %>%  
        dplyr::mutate(age = row_number() + 10,        
                      cycle = age - 9) 
      
      # Modify the dataframe: Add a new row with age = 10 and transition column = 0, and delete the last row (age = 85)
      transition_cases <- transition_cases %>%
        add_row(!!transition_column := 0, age = 10, cycle = 1, .before = 1) %>%  # Add row at the beginning
        slice(-n()) %>%  # Remove the last row
        mutate(age = 10:(10 + n() - 1),  # Adjust age to start from 10
               cycle = age - 9)  # Adjust cycle
    } else {
      # Handle missing transition columns
      warning(paste0("Transition '", transition_column, 
                     "' not found! Using a column of zeros."))
      transition_cases <- tibble(
        !!transition_column := rep(0, nrow(Tot_Trans_per_t_tbl)),  
        age = row_number() + 10,
        cycle = age - 9
      ) %>%
        add_row(!!transition_column := 0, age = 10, cycle = 1, .before = 1) %>%  
        slice(-n()) %>%
        mutate(age = 10:(10 + n() - 1), 
               cycle = age - 9)
    }
    return(transition_cases)
    
    # Case when state1 is a vector (length > 1)
  } else if (length(state1) > 1) {
    transition_columns <- paste0(state1, "->", state2)  
    
    # Handle missing columns and replace them with zeros
    existing_cols <- intersect(transition_columns, colnames(Tot_Trans_per_t_tbl))
    missing_cols <- setdiff(transition_columns, colnames(Tot_Trans_per_t_tbl))
    
    if (length(missing_cols) > 0) {
      warning(paste0("Some transitions not found: ", 
                     paste(missing_cols, collapse = ", "), 
                     ". Using columns of zeros for these."))
    }
    
    # Create missing columns (zeros)
    # Create a tibble for the missing columns (zeros).
    # The operator ' unquote-splice` ("!!!") splices or unpack (corte y pega) 
    # a list or vector into multiple arguments (used with functions of `rlang`).
    # in our case the !!! is used to unpack the list returned by setNames() 
    # and pass it as individual arguments to tibble(). This way, each item in 
    # the list becomes a separate column in the tibble, with the names provided
    # by missing_cols.
    missing_df <- tibble(
      !!!setNames(lapply(missing_cols, function(x) rep(0, nrow(Tot_Trans_per_t_tbl))), missing_cols)
    )
    
    # Combine and process
    # The "unquote" operator unquotes a value or an expression, rather than 
    # treating it as a literal symbol or character string.
    # a) !! (Unquote): Injects a single value or expression into a function. 
    # It is typically used when you want to reference or compute something based
    # on a single variable or expression.
    # b) !!! (Unquote-splice): Injects or "splices" multiple values or elements 
    #from a list or vector into a function. It is used when you need to spread 
    # a list of arguments across multiple positions or inputs.
    transition_cases <- Tot_Trans_per_t_tbl %>%
      dplyr::select(all_of(existing_cols)) %>%
      bind_cols(missing_df) %>%  
      rowwise() %>%
      mutate(!!paste0(state2, "_per_t") := sum(c_across(everything()), na.rm = TRUE)) %>%
      ungroup() %>%
      dplyr::mutate(age = row_number() + 10,  
                    cycle = age - 9) %>%
      add_row(!!paste0(state2, "_per_t") := 0, age = 10, cycle = 1, .before = 1) %>%
      slice(-n()) %>%
      mutate(age = 10:(10 + n() - 1), 
             cycle = age - 9)
    return(transition_cases)
  }
}
################################################################################


################################################################################
## THE MICROSIMULATION MAIN FUNCTION
# Mod: incorporate loop over simulations:
# This version stacks solution of simulations but produces a list with stacked elements
# Parallelize!
ensure_library(parallel)

MicroSim <- function(strategy="natural_history", numb_of_sims = 1,
                     v_M_1, n_i, n_t, v_n, d_c, d_e, TR_out = TRUE, 
                     TS_out = TRUE, Trt = FALSE,  seed = 1, Pmatrix) {
   
  seeds <- sample(1:10000, numb_of_sims, replace = FALSE)  # Generate random seeds
  #simulation_results <- vector("list", numb_of_sims)
  simulation_results <- list() 
  ## Debuging:
  #TR_out = TRUE; TS_out = TRUE; Trt = FALSE; seed = 1
  
  #Arguments:
  # v_M_1:   vector of initial states for individuals
  # n_i:     number of individuals
  # n_t:     total number of cycles to run the model
  # v_n:     vector of health state names
  # d_c:     discount rate for costs
  # d_e:     discount rate for health outcome (QALYs)
  # TR_out:  should the output include a Microsimulation trace? 
  #          (default is TRUE)
  # TS_out:  should the output include a matrix of transitions between states? 
  #          (default is TRUE)
  # Trt:     are the n.i individuals receiving treatment? (scalar with a Boolean
  #          value, default is FALSE)
  # seed:    starting seed number for random number generator (default is 1)
  # Makes use of:
  # Probs:   function for the estimation of transition probabilities
  # Costs:   function for the estimation of cost state vamatrix: Matrix of 
  # tranistion probabilities for each sim cycle.
  # Effs:    function for the estimation of state specific health outcomes (QALYs)
  # Pmatrix: Matrix of transition probabilities for each sim cycle.
  
  # Symptomatic individuals are those who, while in a cancer state (FIGO I-IV),
  # develop symptoms according to the probability vector `figoSymProb`. It is 
  # assumed that all individuals who develop symptoms will visit a doctor. This 
  # event incurs a one-time, lifetime cost (applicable only once, upon diagnosis).
  #symptomatics <- data.frame()
  
  ## Initialize an empty list to store results from each simulation
  #simulation_results <- list()
  #seeds <- sample(1:10000, numb_of_sims, replace = FALSE)  # Generate random seeds
  
  ############################################################################
  my_age_prob_matrix_func <- function(my_Prob_matrix, my_age_in_loop) {
    my_age_prob_matrix <- my_Prob_matrix %>% 
      dplyr::filter(Lower <= my_age_in_loop  &
                      Larger >= my_age_in_loop) 
  }
  ############################################################################
 
  
  # Function for a single simulation
  run_single_simulation <- function(sim, seed, strategy, v_M_1, n_i, n_t, v_n, 
                                    d_c, d_e, TR_out, TS_out, Trt, Pmatrix) {
    set.seed(seed) 
    
    #cat("Running simulation", sim, "with seed", seeds[sim], "\n")
    
    # Initialize objects for storing simulation results
    symptomatics <-
      data.frame(ID = integer(), TimeStep = integer(), 
                 DiagnosedState = character(), 
                 RecoveredFromState = logical(), stringsAsFactors = FALSE)
    
    v_dwc <- 1 / (1 + d_c) ^ (0:(n_t-1)) # calculate the cost discount weight based
    v_dwe <- 1 / (1 + d_e) ^ (0:(n_t-1)) # calculate the QALY discount weight based 
    #v_dwe <- 1 / (1 + d_e) ^ (1:(n_t))   # calculate the QALY discount weight based 
    
    # Create the matrix capturing the state name/costs/health outcomes 
    # for all individuals at each time point:
    #m_M <- m_C <- m_E <-  matrix(nrow = n_i, ncol = (n_t + 1), 
    m_M <- m_C <- m_E <- 
      matrix(nrow = n_i, ncol = (n_t), 
             dimnames = list( 1:n_i, 
                              #paste0("cycle_", 1:(n_t + 1), sep = "")))  
                              paste0("cycle_", 1:(n_t), sep = "")))  
    
    m_M[, 1] <- v_M_1  # Set initial health state   
    
    seed <- seeds[sim]
    #seed <- 17
    cat ("This is simulation's seed:  ", seed, "\n")
    set.seed(seed) # set the seed for every individual for the random number generator
    
    
    m_C[, 1] <- Costs_per_Cancer_Diag(M_it = m_M[, 1], # estimate costs per individual for the 
                                      symptomatics = symptomatics,
                                      time_iteration = 1,
                                      cost_Vec = cost_Vec, # initial health state
                                      Trt)             
    
    m_E[, 1] <- Effs(m_M[, 1], Trt, utilityCoefs = utilityCoefs)  # estimate QALYs
                                                                  # per individual 
                                                                  # for the initial
                                                                  # health state  
    stored_list <- list()
    ######################## run over all the cycles ############################# 
    #for (t in 1:(n_t)) {
    # Debugging:
    #for (t in 2:(n_t)) {
    for (t in 1:(n_t - 1)) {
      ############################################################################
      # Select the transition matrix based on the cycle `n_t`:
      # Since our age intervals start at 10 years old,
      # `age_in_loop` is calculated as `t + 10 * age_factor(cycle_period)`.
      #age_in_loop <- t + 10
      #age_in_loop <- t + 8
      age_in_loop <- t + 9
      ########################################################################
      
      # update/correct n_s (<<- let change variable from inside a function):
      n_s  <<- length(v_n)  
      
      ######################################################################### 
      ## This piece of code has been moved to run after the transition has 
      ## been performed, and to have entries for states of the next time step, 
      ## i.e. a t + 1 as the computation of the cost is performed for t+1 during
      # the loop t.
      ##new code:
      ## Diagnose (or appearance of symptomatics):
      #new_entries <- diagnose_column(m_M[, t], t)
      #
      #if (!is.null(new_entries)) {
      #  stored_list[[t]] <- new_entries
      #}
      #if (nrow(new_entries) > 0) {
      #  symptomatics <- bind_rows(symptomatics, new_entries)
      #}
      ######################################################################## 
      
      
      ########################################################################    
      my_age_prob_matrix <- 
        my_age_prob_matrix_func(my_Prob_matrix = my_Probs, 
                                ## WHY add 1 to age_in_loop??
                                #my_age_in_loop = (age_in_loop + 1))
                                my_age_in_loop = (age_in_loop))
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
      ##########################################################################
      
      
      ##########################################################################
      ## This piece of coding is not needed as long as the diagnose/symptomatics
      ## determination if done after the computation of  transitions
      ## m_M[, t + 1] <- update_column(m_M[, t], new_entries)
      #next_col <- m_M[, t + 1]
      #next_col <- update_column(m_M[, t], new_entries, next_col)
      ## Ensure next_col updates are preserved after sampling
      #m_M[, t + 1] <- ifelse(next_col == "Survival", "Survival", m_M[, t + 1])
      ##########################################################################
      
      
      
      ##########################################################################
      #new code:
      # Diagnose (or appearance of symptomatics):
      new_entries <- diagnose_column(m_M[, t + 1], (t + 1))
      
      if (!is.null(new_entries)) {
        stored_list[[t+1]] <- new_entries
      }
      if (nrow(new_entries) > 0) {
        symptomatics <- bind_rows(symptomatics, new_entries)
      }
      ########################################################################## 
      
      
      ##########################################################################    
      ## Costs per CC diagnose at time t + 1.
      
      # Estimate costs per individual during cycle t + 1 conditional on treatment:
      m_C[, t + 1] <-                              
        Costs_per_Cancer_Diag(M_it = m_M[, t + 1],  
                              symptomatics = symptomatics,
                              #time_iteration = t,
                              time_iteration = (t+1),
                              cost_Vec = cost_Vec,    
                              Trt)            
      
      m_E[, t + 1] <- # estimate QALYs per individual during cycle t + 1
        Effs( m_M[, t + 1], Trt, 
              utilityCoefs = utilityCoefs)                   
    ############################################################################    
      
      # Conditional on treatment
      cat('\r', paste(round(t/n_t * 100),          # display the 
                      "% done", sep = " "))        # progress of  the simulation                    
      
    }##################### close loop for cycles ###############################  
    
    # Combine stored entries into a single data frame
    symptomatics <- bind_rows(stored_list)
    tc_disc <- m_C[,1:n_t] %*% v_dwc       # total (discounted) cost per individual
    te_disc <- m_E[,1:n_t] %*% v_dwe       # total (discounted) QALYs per individual 
    
    tc_undisc <- m_C[,1:n_t] %*% rep(1, n_t) # total (discounted) cost per individual
    te_undisc <- m_E[,1:n_t] %*% rep(1, n_t) # total (discounted) QALYs per individual 
    
    tc_hat_disc <- mean(tc_disc)        # average (discounted) cost 
    te_hat_disc <- mean(te_disc)        # average (discounted) QALYs
    tc_hat_undisc <- mean(tc_undisc)    # average (discounted) cost 
    te_hat_undisc <- mean(te_undisc)    # average (discounted) QALYs
    
    ##cat("I'm still here, debugging! \n")
    #browser()
    
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
    #rm(m_M, m_C, m_E)
    
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
    simulation_results[sim] <- list(results)
  } # end of `numb_of sims` loop
  
  #}  # end of `numb_of_sims` loop
  
  #return(simulation_results)
  
  # stack results
  #source("./R/Sumarize_results_by_Strategy_Func.R")
  stacked_results <- 
    summarize_results_by_Strategy(results_list = simulation_results, 
                                  numb_of_sims = numb_of_sims)
  return(stacked_results)
} # end of MicroSim function
################################################################################


################################################################################
########################## Run the simulation ##################################
## START SIMULATION
p = Sys.time()
# run for no treatment
sim_no_trt  <- MicroSim(strategy = "natural_history",numb_of_sims = 50, 
                        v_M_1 = v_M_1, n_i = n_i, n_t = n_t, v_n = v_n, 
                        d_c = d_c, d_e = d_e, TR_out = TRUE, TS_out = TRUE, 
                        Trt = FALSE, seed = 1, Pmatrix = Pmatrix)

# Load computed simulation if needed here:
#sim_no_trt <- readRDS(file = "./data/stacked_sims_100x10E6x75_FULL_IMPLEMENTATION.rds")
#sim_no_trt <- readRDS(file = "./data/stacked_sims_100x10E6x75.rds")
#sim_no_trt <- readRDS(file = "./data/stacked_sims_20x10E6x75_20241009.rds")
#sim_no_trt <- readRDS(file = "./data/stacked_sims_10x10E6x75_20241016.rds")

comp.time = Sys.time() - p
comp.time %>% print()
################################################################################
################################################################################



################################################################################
################################################################################
###################################
## Post-simulation Computations: ##
###################################
################################################################################
################################################################################
# # Prevalence is defined as number of infected divided by total alive individuals
# # for that cycle/time step
mean_prevalence_func <- function(sim_stalked_result, my_Probs) {
  # Extract unique age intervals and ensure Larger doesn't exceed 84
  age_intervals <- my_Probs %>% 
    select(Lower, Larger) %>% 
    unique() %>% 
    dplyr::mutate(Larger = ifelse(Larger == 85, 84, Larger)) %>% # Ensure Larger is capped at 84
    arrange(Lower)
  
  # Create a vector of breaks for the intervals
  breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
  
  # Create labels for the intervals
  labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
  
  # Compute prevalence and average it by age intervals
  df <- sim_stalked_result[[1]]$TR %>% 
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
  sim_stalked_result[[1]]$mean_HPV_prevalence_per_age_interval <- df
  return(sim_stalked_result)
}

# Concatenate the prevalence to the sim result 
mean_prevalence_result <-
  mean_prevalence_func(sim_stalked_result = sim_no_trt, my_Probs = my_Probs)  


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
    dplyr::mutate(Larger = ifelse(Larger == 85, 84, Larger)) %>%  # Cap at 84
    arrange(Lower)
  
  # Create a vector of the breaks for the intervals
  breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
  
  # Create labels for the intervals
  labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
  
  # Define all previous epi states
  all_previous_states <- v_n[1 : ( which(v_n == state) - 1)]
  
  my_incidence_df <- sim_stalked_result[[1]]$TR
  
  # Define the new state for incidence calculation
  new_state <- paste0("new_", state)
  new_state_df <- sim_stalked_result[[1]][new_state][[1]]
  
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
  
  # Store the result in the list
  sim_stalked_result[[1]][[paste0("mean_incidence_", state, "_per_age_interval")]] <- df
  return(sim_stalked_result)
}
################################################################################

# Computing incidences:
incidence_states_to_compute <- c("CIN1", "CIN2", "CIN3") 

# Initialize the result with the original structure
mean_incidence_result <- mean_prevalence_result

# Apply the incidence function to each state and update the result structure
for (my_state in incidence_states_to_compute) {
  #print(my_state)
  mean_incidence_result <- 
    mean_incidence_func(sim_stalked_result = mean_incidence_result, 
                        state = my_state, my_Probs = my_Probs)
}
################################################################################


################################################################################
# Computing Cervix Cancer incidence:
mean_CC_incidence_func <- function(sim_stalked_result, my_Probs) {
  # Extract unique age intervals
  #age_intervals <- my_Probs %>% 
  #  select(Lower, Larger) %>% 
  #  unique() %>% 
  #  ## making larger interval equals 84
  #  #dplyr::mutate(Larger = ifelse(Larger==85, 84, Larger)) %>%
  #  arrange(Lower)
  
  age_intervals <- my_Probs %>% 
    dplyr::select(Lower, Larger) %>% 
    unique() %>% 
    dplyr::mutate(Larger = ifelse(Larger == 85, 84, Larger)) %>%  # Cap at 84
    arrange(Lower)
  
  # Create a vector of the breaks for the intervals
  breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
  
  # Create labels for the intervals
  labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
  sim_stalked_result[[1]]$new_Cancer <- sim_stalked_result[[1]]$new_Cancer #%>%
  # dplyr::select(-sim.1)
  
  # Compute prevalence and average it by age intervals
  df <- merge(sim_stalked_result[[1]]$TR, 
              #sim_stalked_result[[1]]$new_Cancer, by = c("sim", "age", "cycle")) %>%
              sim_stalked_result[[1]]$new_Cancer, by = c("sim", "age")) %>%
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
  
  #return(df)
  sim_stalked_result[[1]]$mean_CC_incidence <- df
  return(sim_stalked_result)
}
################################################################################

# Initialize the result with the original structure
mean_CC_incidence_result <- mean_incidence_result

# Concatenate the prevalence to the sim result 
mean_CC_incidence_result <-
  mean_CC_incidence_func(sim_stalked_result = mean_CC_incidence_result,
                         my_Probs = my_Probs)  

################################################################################
# Computing Mortality:
# A. Cancer-related deaths at certain age (cycle) / total alive at that age (cycle)
# B. Cancer-unrelated deaths at certain age (cycle) / total alive at that age (cycle)

################################################################################
# A. Cancer-related Deaths per age
mean_CC_mortality_func <- function(sim_stalked_result, my_Probs) {
  # Extract unique age intervals
  #age_intervals <- my_Probs %>% 
  #  select(Lower, Larger) %>% 
  #  unique() %>% 
  #  ## making larger interval equals 84
  #  #dplyr::mutate(Larger = ifelse(Larger==85, 84, Larger)) %>%
  #  arrange(Lower)
  
  age_intervals <- my_Probs %>% 
    dplyr::select(Lower, Larger) %>% 
    unique() %>% 
    dplyr::mutate(Larger = ifelse(Larger == 85, 84, Larger)) %>%  # Cap at 84
    arrange(Lower)
  
  # Create a vector of the breaks for the intervals
  breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
  
  # Create labels for the intervals
  labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
  
  ## Compute prevalence and average it by age intervals
  #df <- sim_stalked_result[[1]]$TR %>% 
  #  #dplyr::select(sim, cycle, age, H, HR.HPV.infection) %>% 
  #  dplyr::select(everything()) %>% 
  #  dplyr::mutate(total_alive = H + HR.HPV.infection + CIN1 + CIN2 + CIN3 +
  #                  FIGO.I + FIGO.II + FIGO.III + FIGO.IV + Survival) %>%
  #  dplyr::mutate(CC_mortality = (CC_Death / total_alive) * 10^5) %>% 
  #  dplyr::mutate(age_interval = cut(age, breaks = breaks, labels = labels, right = FALSE)) %>% 
  #  dplyr::group_by(age_interval) %>% 
  #  dplyr::summarise(CC_mean_mortality = mean(CC_mortality, na.rm = TRUE)) %>% 
  #  dplyr::ungroup()
  
  ################# TESTING #################################################
  # Define the age range you want to keep
  age_range <- 10:84
  
  # Left join sim_no_trt[[1]]$TR with sim_no_trt[[1]]$new_CC_Death by age
  df <- sim_no_trt[[1]]$TR %>%
    left_join(sim_no_trt[[1]]$new_CC_Death %>%
                dplyr::select(sim, age, CC_Death_per_t), 
              by = c("sim", "age"))  %>% #, relationship = "many-to-many") %>%
    
    # Filter for ages in the desired range
    dplyr::filter(age %in% age_range) %>%
    
    # Fill missing CC_Death_per_t with zeros (for cases where the age doesn't exist)
    dplyr::mutate(CC_Death_per_t = coalesce(CC_Death_per_t, 0)) %>%
    
    # Compute total_alive one previous time stpe / cylce:
    dplyr::mutate(total_alive_lagged = dplyr::lag(H + HR.HPV.infection + CIN1 + CIN2 + CIN3 +
                                                    FIGO.I + FIGO.II + FIGO.III + FIGO.IV + Survival), n=1) %>%
    
    # Compute CC_mortality based on CC_Death_per_t and total_alive
    dplyr::mutate(CC_mortality = (CC_Death_per_t / total_alive_lagged) * 10^5) %>%
    
    # Compute age intervals and average CC_mortality by age intervals
    dplyr::mutate(age_interval = cut(age, breaks = breaks, labels = labels, right = FALSE)) %>%
    dplyr::group_by(age_interval) %>%
    dplyr::summarise(CC_mean_mortality = mean(CC_mortality, na.rm = TRUE)) %>%
    dplyr::ungroup()
  ################# TESTING #################################################
  
  #return(df)
  sim_stalked_result[[1]]$CC_mean_mortality <- df
  return(sim_stalked_result)
}
################################################################################

# Initialize the result with the original structure
mean_CC_mortality_result <- mean_CC_incidence_result 
# Concatenate the prevalence to the sim result 
mean_CC_mortality_result <-
  mean_CC_mortality_func(sim_stalked_result = mean_CC_mortality_result, my_Probs = my_Probs)  

################################################################################
# A.2 Cancer-related Deaths (per differences) per age
mean_CC_mortality_by_diff_func <- function(sim_stalked_result, my_Probs) {
  # Extract unique age intervals
  #age_intervals <- my_Probs %>% 
  #  select(Lower, Larger) %>% 
  #  unique() %>% 
  #  ## making larger interval equals 84
  #  #dplyr::mutate(Larger = ifelse(Larger==85, 84, Larger)) %>%
  #  arrange(Lower)
  
  age_intervals <- my_Probs %>% 
    dplyr::select(Lower, Larger) %>% 
    unique() %>% 
    dplyr::mutate(Larger = ifelse(Larger == 85, 84, Larger)) %>%  # Cap at 84
    arrange(Lower)
  
  # Create a vector of the breaks for the intervals
  breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
  
  # Create labels for the intervals
  labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
  
  ################# TESTING #################################################
  # Define the age range you want to keep
  age_range <- 10:84
  
  # Left join sim_no_trt[[1]]$TR with sim_no_trt[[1]]$new_CC_Death by age
  df <- sim_no_trt[[1]]$TR %>%
    left_join(sim_no_trt[[1]]$CC_Death_by_diff %>%
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
  ################# TESTING #################################################
  
  #return(df)
  sim_stalked_result[[1]]$CC_by_diff_mean_mortality <- df
  return(sim_stalked_result)
}
################################################################################


### TESTING ###
mean_CC_mortality_by_diff_result <- mean_CC_mortality_result
mean_CC_mortality_by_diff_result <-
  mean_CC_mortality_by_diff_func(sim_stalked_result =
                                   mean_CC_mortality_by_diff_result, my_Probs = 
                                   my_Probs)  

################################################################################
# B. Cancer-unrelated Mortality
other_mean_mortality_func <- function(sim_stalked_result, my_Probs) {
  
  # Extract unique age intervals
  #age_intervals <- my_Probs %>% 
  #  select(Lower, Larger) %>% 
  #  unique() %>% 
  #  ## making larger interval equals 84
  #  #dplyr::mutate(Larger = ifelse(Larger==85, 84, Larger)) %>%
  #  arrange(Lower)
  
  age_intervals <- my_Probs %>% 
    dplyr::select(Lower, Larger) %>% 
    unique() %>% 
    dplyr::mutate(Larger = ifelse(Larger == 85, 84, Larger)) %>%  # Cap at 84
    arrange(Lower)
  
  # Create a vector of the breaks for the intervals
  breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
  
  # Create labels for the intervals
  labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
  
  # Compute prevalence and average it by age intervals
  df <- sim_stalked_result[[1]]$TR %>% 
    #dplyr::select(sim, cycle, age, H, HR.HPV.infection) %>% 
    dplyr::select(everything()) %>% 
    dplyr::mutate(total_alive = H + HR.HPV.infection + CIN1 + CIN2 + CIN3 +
                    FIGO.I + FIGO.II + FIGO.III + FIGO.IV + Survival) %>%
    dplyr::mutate(other_mortality = (Other.Death / total_alive) * 10^5) %>% 
    dplyr::mutate(age_interval = cut(age, breaks = breaks, labels = labels, right = FALSE)) %>% 
    dplyr::group_by(age_interval) %>% 
    dplyr::summarise(other_mean_mortality = mean(other_mortality, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  
  #return(df)
  sim_stalked_result[[1]]$other_mean_mortality <- df
  return(sim_stalked_result)
}
################################################################################

# Initialize the result with the original structure
other_mean_mortality_result <- mean_CC_mortality_by_diff_result
# Concatenate the prevalence to the sim result 
other_mean_mortality_result <-
  other_mean_mortality_func(sim_stalked_result = 
                              other_mean_mortality_result, my_Probs = my_Probs)  

## save the results
#saveRDS(object = other_mean_mortality_result, file = "./data/stacked_sims_10x10E6x75_20241016.rds")


### ----convert .Rmd to .R
#library(knitr)
## purl("your_script.Rmd", output = "your_script.R")
## example:
#purl("Cervix_MicroSim_RMarkdown_v.072_B.Rmd", output = "cervix_microSim_stacked_list.R")
#purl("Cervix_MicroSim_RMarkdown_v.072_B.Rmd", output = "cervix_microSim_stacked_list_B.R")


## ----cost-efectiveness
####################### Cost-effectiveness analysis #############################
## store the mean costs (and MCSE) of each strategy in a new variable C (vector costs)
#v_C  <- c(sim_no_trt$tc_hat_disc, sim_trt$tc_hat_disc) 
#sd_C <- c(sd(sim_no_trt$tc_disc), sd(sim_trt$tc_disc)) / sqrt(n_i)
## store the mean QALYs (and MCSE) of each strategy in a new variable E (vector effects)
#v_E  <- c(sim_no_trt$te_hat_disc, sim_trt$te_hat_disc)
#sd_E <- c(sd(sim_no_trt$te_disc), sd(sim_trt$te_disc)) / sqrt(n_i)
#
#delta_C <- v_C[2] - v_C[1]                   # calculate incremental costs
#delta_E <- v_E[2] - v_E[1]                   # calculate incremental QALYs
## Monte Carlo Squared Error (MCSE) of incremental costs:
#sd_delta_E <- sd(sim_trt$te - sim_no_trt$te) / sqrt(n_i) 
## Monte Carlo Squared Error (MCSE) of incremental QALYs:
#sd_delta_C <- sd(sim_trt$tc_disc - sim_no_trt$tc_disc) / sqrt(n_i) 
#ICER    <- delta_C / delta_E                 # calculate the ICER
#results <- c(delta_C, delta_E, ICER)         # store the values in a new variable
#
## Create full incremental cost-effectiveness analysis table
#table_micro <- data.frame(
#  c(round(v_C, 0),  ""),           # costs per arm
#  c(round(sd_C, 0), ""),           # MCSE for costs
#  c(round(v_E, 3),  ""),           # health outcomes per arm
#  c(round(sd_E, 3), ""),           # MCSE for health outcomes
#  c("", round(delta_C, 0),   ""),  # incremental costs
#  c("", round(sd_delta_C, 0),""),  # MCSE for incremental costs
#  c("", round(delta_E, 3),   ""),  # incremental QALYs 
#  c("", round(sd_delta_E, 3),""),  # MCSE for health outcomes (QALYs) gained
#  c("", round(ICER, 0),      "")   # ICER
#)
## name the rows:
#rownames(table_micro) <- c(v_Trt, "* are MCSE values")  
## name the columns:
#colnames(table_micro) <-  
#  c("Costs", "*",  "QALYs", "*", "Incremental Costs",
#    "*", "QALYs Gained", "*", "ICER")
#table_micro  # print the table 


## ----plot curves
## This R chunk is a plot routine (not part of the main program):
library(RColorBrewer)
#ensure_library("RColorBrewer")
# Convert matrix to data frame
micro_sim_df <- sim_no_trt[[1]]$TR
micro_sim_df <- other_mean_mortality_result[[1]]$TR

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Average the specified columns by age across all simulations
averaged_micro_sim_df <- micro_sim_df %>%
  group_by(age) %>%
  summarise(across(c(CIN1, CIN2, CIN3, 
                     FIGO.I, FIGO.II, FIGO.III, FIGO.IV,
                     Survival, CC_Death, Other.Death), mean)) %>%
  ungroup()

# Reshape the data to long format for easier plotting with ggplot2
long_micro_sim_df <- averaged_micro_sim_df %>%
  pivot_longer(cols = c(CIN1, CIN2, CIN3, FIGO.I, FIGO.II, FIGO.III, FIGO.IV),
               #pivot_longer(cols = c( FIGO.IV, CC_Death),
               #Survival, CC_Death, Other.Death), 
               names_to = "Stage", 
               values_to = "Average")

# Plot the data using ggplot2
ggplot(long_micro_sim_df, aes(x = age, y = Average, color = Stage)) +
  geom_line() +
  labs(title = "Averaged CIN and FIGO Stages by Age Across All Simulations",
       x = "Age",
       y = "Average Count",
       color = "Stage") +
  theme_minimal()


## ----loading markov result
if (!require("readxl")) install.packages("readxl")
library(readxl)
# This R chunk is a plot routine (not part of the main program):


if (!require("readxl")) install.packages("readxl")
library(readxl)
#markov <-
#  readxl::read_excel("Q:/my_Q_docs/Cervix_MicroSim/CervixMicroSim_Carlos/carlos__Krijkamp_ver/data/Sortida_NoIntervencio.xlsx", sheet = "NH")
markov_df <- readxl::read_excel("./data/Sortida_NoIntervencio.xlsx")

markov_df <- markov_df %>% mutate(age = Step + 10)

# Reshape the data into long format
markov_df_long_data <- markov_df %>%
  pivot_longer(cols = c(HR.HPV.infection, CIN1, CIN2, CIN3, FIGO.I, FIGO.II,
                        #FIGO.III, FIGO.IV, Survival, CC_Death, Other.Death),
                        FIGO.III, FIGO.IV, Survival, CC_Death),
               names_to = "Health state",
               values_to = "value")

# Plot the data
ggplot(markov_df_long_data, aes(x = age, y = value, color = `Health state`)) +
  geom_line(linewidth=1, alpha=0.7) +
  labs(x = "Age", y = "Value", color = "Health state") +
  ggtitle(expression(paste("Markov cohort simulation for ", 10^6, " individuals"))) + 
  theme_minimal()  # Optional: customize the theme

################################################################################
# Comparing with the microsimulation:
# column bind `micro_df` and `markov` by the `age`column:
merged_df <- left_join(averaged_micro_sim_df, markov_df, by = 'age')

# get rid of NAs
merged_df <- merged_df %>% na.omit()
################################################################################

################################################################################
# compare them:
# Reshape the data into long format:
long_merged_data <- merged_df %>%
  pivot_longer(cols = c(#HR.HPV.infection.x, HR.HPV.infection.y, 
    CIN1.x, CIN1.y, 
    CIN2.x, CIN2.y, 
    CIN3.x, CIN3.y, 
    FIGO.I.x,FIGO.I.y,  
    FIGO.II.x, FIGO.II.y,
    #FIGO.III, FIGO.IV, Survival, CC_Death, Other.Death),
    FIGO.III.x, FIGO.III.y, FIGO.IV.x, FIGO.IV.y,
    #Survival.x, Survival.y, 
    CC_Death.y, CC_Death.y),
    names_to = "Health state",
    values_to = "value")

# Plot the data
ggplot(long_merged_data, aes(x = age, y = value, color = `Health state`)) +
  geom_line(linewidth=1, alpha=0.7) +
  labs(x = "Age", y = "Value", color = "Health state") +
  ggtitle(expression(paste("Markov cohort vs  Microsimulation for ", 10^6, " individuals"))) + 
  theme_minimal()  # Optional: customize the theme

################################################################################


## ----incidences, prevalences, and mortalities-------------------------------------------------------------------------------------------------------------------------------------------
# Markov:
markov_CN1_incidences <- c(0.00000, 204.73492, 981.96179, 1368.24200, 3006.85782, 33.48096, 1362.96678, 459.48051, 697.84223, 794.33833, 223.00222, 246.23082, 176.02167, 126.22963, 53.70939)
markov_CN2_incidences <- c(0.000000, 6.165629, 54.767952, 140.309815, 216.568392, 1476.306267, 1579.728160, 1298.914564, 466.596151, 637.661611, 442.298632, 304.784447, 250.953880, 165.628020, 116.925192)
markov_CN3_incidences <- c(0.000000, 2.090325, 9.597415, 44.467676, 148.972191, 0.000000, 3.550684, 91.881726, 12.505042, 68.377446, 25.802481, 7.952667, 1.174088, 1.177840, 2.638642)
markov_CC_incidences  <- c(0.000000, 0.000000, 0.000000, 5.520938, 8.360544, 13.282380, 22.906871, 20.825560, 15.867891, 32.483846, 8.962389, 17.681771, 11.737615, 17.354646, 14.582775)
markov_HPV_prevalences <- c(0.000000000, 0.343480414, 0.377634762, 0.087223460, 0.307341403, 0.030196332, 0.050562845, 0.050151668, 0.082952596, 0.046644059, 0.018532077, 0.034193076, 0.016407832, 0.015039027, 0.003217326)
markov_CC_mortality <- c(0.000000e+00, 0.000000e+00, 0.000000e+00, 2.977975e-06, 1.574920e-05, 2.715056e-05, 5.489929e-05, 7.284815e-05, 1.057494e-04, 5.076268e-05, 7.517773e-05, 4.960943e-05, 4.802468e-05, 4.210457e-05, 4.837655e-05) * 10^5

# Micro im:
microSim_CN1_incidences          <- other_mean_mortality_result[[1]]$mean_incidence_CIN1_per_age_interval
microSim_CN2_incidences          <- other_mean_mortality_result[[1]]$mean_incidence_CIN2_per_age_interval
microSim_CN3_incidences          <- other_mean_mortality_result[[1]]$mean_incidence_CIN3_per_age_interval
microSim_CC_incidences           <- other_mean_mortality_result[[1]]$mean_CC_incidence
microSim_HPV_prevalences         <- other_mean_mortality_result[[1]]$mean_HPV_prevalence_per_age_interval
microSim_CC_mortality            <- other_mean_mortality_result[[1]]$CC_mean_mortality
microSim_CC_by_diff_mortality    <- other_mean_mortality_result[[1]]$CC_by_diff_mean_mortality



## ----ploting incidences and prevalences
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Define age groups
age_groups <- factor(c("10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", 
                       "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", 
                       "75-79", "80-84"), 
                     levels = c("10-14", "15-19", "20-24", "25-29", "30-34", "35-39", 
                                "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", 
                                "70-74", "75-79", "80-84"))

# Create data frames from your vectors and the MicroSim data
# Replace `microSim_...` with the actual data from your `other_mean_mortality_result`

markov_data <- data.frame(
  age = age_groups,
  markov_CN1_incidences = markov_CN1_incidences,
  markov_CN2_incidences = markov_CN2_incidences,
  markov_CN3_incidences = markov_CN3_incidences,
  markov_CC_incidences = markov_CC_incidences,
  markov_HPV_prevalences = markov_HPV_prevalences,
  markov_CC_mortality = markov_CC_mortality,
  # Assign the same values from markov_CC_mortality to markov_CC_by_diff_mortality
  markov_CC_by_diff_mortality <- markov_CC_mortality
)

# Ensure all columns in markov_data are numeric
markov_data[] <- lapply(markov_data, function(x) {
  if (is.factor(x)) {
    as.character(x)
  } else {
    as.numeric(x)
  }
})


# For the MicroSim data, ensure columns are numeric if needed
# You might need to extract these from the list manually and convert them

# Example conversion if you have microSim data as tibbles
microSim_data <- data.frame(
  age = age_groups,
  microSim_CN1_incidences = as.numeric(other_mean_mortality_result[[1]]$mean_incidence_CIN1_per_age_interval$mean_incidence_CIN1),
  microSim_CN2_incidences = as.numeric(other_mean_mortality_result[[1]]$mean_incidence_CIN2_per_age_interval$mean_incidence_CIN2),
  microSim_CN3_incidences = as.numeric(other_mean_mortality_result[[1]]$mean_incidence_CIN3_per_age_interval$mean_incidence_CIN3),
  microSim_CC_incidences = as.numeric(other_mean_mortality_result[[1]]$mean_CC_incidence$CC_mean_incidence),
  microSim_HPV_prevalences = as.numeric(other_mean_mortality_result[[1]]$mean_HPV_prevalence_per_age_interval$prevalence),
  microSim_CC_mortality = as.numeric(other_mean_mortality_result[[1]]$CC_mean_mortality$CC_mean_mortality),
  microSim_CC_by_diff_mortality = as.numeric(other_mean_mortality_result[[1]]$CC_by_diff_mean_mortality$CC_by_diff_mean_mortality)
)

# Ensure all columns in microSim_data are numeric
microSim_data[] <- lapply(microSim_data, function(x) {
  if (is.factor(x)) {
    as.character(x)
  } else {
    as.numeric(x)
  }
})

# Reshape data to long format
markov_long <- markov_data %>% 
  pivot_longer(-age, names_to = "measure", values_to = "value") %>% 
  mutate(model = "Markov")

microSim_long <- microSim_data %>% 
  pivot_longer(-age, names_to = "measure", values_to = "value") %>% 
  mutate(model = "MicroSim")

# Combine data
combined_data <- bind_rows(markov_long, microSim_long)

# Plotting function
plot_comparison <- function(data, measure_name) {
  ggplot(data %>% filter(grepl(measure_name, measure)), 
         aes(x = age, y = value, fill = model)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste(measure_name, "Comparison"),
         x = "Age Group",
         y = measure_name) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create plots for each measure
plot_CN1_incidences <- plot_comparison(combined_data, "CN1_incidences")
plot_CN2_incidences <- plot_comparison(combined_data, "CN2_incidences")
plot_CN3_incidences <- plot_comparison(combined_data, "CN3_incidences")
plot_CC_incidences <- plot_comparison(combined_data, "CC_incidences")
plot_HPV_prevalences <- plot_comparison(combined_data, "HPV_prevalences")
plot_CC_mortality <- plot_comparison(combined_data, "CC_mortality")
plot_CC_by_diff_mortality <- plot_comparison(combined_data, "CC_by_diff_mortality")

# Display plots
print(plot_CN1_incidences)
print(plot_CN2_incidences)
print(plot_CN3_incidences)
print(plot_CC_incidences)
print(plot_HPV_prevalences)
print(plot_CC_mortality)
print(plot_CC_by_diff_mortality)


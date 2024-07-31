## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----preamble, echo=FALSE, results='hide', message=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------
################################################################################
# This code is a modified version of the original code from:
# [https://github.com/DARTH-git/Microsimulation-tutorial] (Krijkamp et al 2018 
# Sick-Sicker model).
# Modifications by: Carlos Dommar D'Lima - carlos.dommar@gmail.com
# This code extends the "sick-sicker" model of the original authors to a
# multi-state cervix cancer model
################################################################################
rm(list = ls())
library(tidyverse)

# Sources:
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
  #dplyr::filter(Age.group == "25-29") %>% # choose one for test
  as_tibble() # I need a tibble to use 'rename' function down there:

# tidying up a bit the transition matrix:
my_Probs <- my_Probs %>% dplyr::rename("H" = "Well")
# Rename the 'old_name' column to 'new_name'

my_Probs <- my_Probs %>% as.data.frame() # convert back to data.frame (no needed?)

###############################################################
# Function to extract and convert numbers from factor levels
extract_numbers <- function(range_factor) {
  range_string <- as.character(range_factor)
  numbers <- as.numeric(unlist(strsplit(range_string, "-")))
  return(numbers)
}
###############################################################

# Before apply the function, convert Age.group from factor to character
# my_Probs$Age.group <- as.character(my_Probs$Age.group)
# Apply the function to the Range column and create new columns
my_Probs$Lower  <- sapply(my_Probs$Age.group, function(x) extract_numbers(x)[1])
my_Probs$Larger <- sapply(my_Probs$Age.group, function(x) extract_numbers(x)[2])
# For the last cycle/iteration we need to adjust the last transition matrix:
my_Probs$Larger <- 
  ifelse(my_Probs$Larger == max(my_Probs$Larger), my_Probs$Larger + 1, my_Probs$Larger) 


## ----my_Probs, echo=TRUE, results='markup', message=FALSE, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------
library("printr")
my_Probs %>%
  head()
knitr::kable(my_Probs)


## ----model parameters-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
n_i <- 10^6                 # number of simulated individuals
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



## ----functions, include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ----sampling function----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# efficient implementation of the rMultinom() function of the Hmisc package #### 
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
                                            
  
  ##############################################################################
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


## ----probability function-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(tidy = TRUE, out.width = 60)
########################## Probability function #################################
## The Probs function that updates the transition probabilities of every cycle:
Probs <- function(M_it, my_Probs) {
  n_s <- length(v_n)
  n_i <- length(M_it)
  m_P_it <- matrix(NA, n_s, n_i) 
  rownames(m_P_it) <- v_n
  
  for (i in 1:length(v_n)) {
    state_mask <- !is.na(M_it) & M_it == v_n[i]
    
    if (sum(state_mask) > 0) {
      m_P_it[, state_mask] <- lapply(X = v_n, function(x) trans_prb(P = my_Probs, state1 = v_n[i], state2 = x)) %>%
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


## ----costs function, tidy=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Costs function
# The `Costs_per_Cancer_Diag` function estimates the costs of a diagnose 
# individual due to cancer symptoms (FIGO.I-IV) at every cycle. 
# This cost is only charged once in the patient's lifetime.
# NOTE: need to decide if the cost is applied on current time `t` or `t+1` as it is now.
Costs_per_Cancer_Diag <- function (M_it, cost_Vec, symptomatics, time_iteration, Trt = FALSE) {
  c_it <- rep(0, length(M_it))
  ci_t <- 0
  if(nrow(symptomatics) > 0 ) {
    c_it[symptomatics %>% filter(DiagnosedState == "FIGO.I" & TimeStep == time_iteration) %>% 
           select(ID) %>% as.list() %>% unlist()] <- cost_Vec[which(v_n %in% "FIGO.I")]
    
    c_it[symptomatics %>% filter(DiagnosedState == "FIGO.II" & TimeStep == time_iteration) %>% 
           select(ID) %>% as.list() %>% unlist()] <- cost_Vec[which(v_n %in% "FIGO.II")]
    
    c_it[symptomatics %>% filter(DiagnosedState == "FIGO.III" & TimeStep == time_iteration) %>% 
           select(ID) %>% as.list() %>% unlist()] <- cost_Vec[which(v_n %in% "FIGO.III")]
    
    c_it[symptomatics %>% filter(DiagnosedState == "FIGO.IV" & TimeStep == time_iteration) %>% 
           select(ID) %>% as.list() %>% unlist()] <- cost_Vec[which(v_n %in% "FIGO.IV")]
  }
  return(c_it)              		                           # return the costs
}


## ----qalys function, tidy=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Health outcome function 
# new version
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


## ----time period related functions, echo=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------
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

######### WORK IN PROGRESS #################
############################################
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


## ----symptoms, echo=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ----new cases, echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to compute (total) new cases of a state across the cycles:
new_cases <- function(state) {
  # It receives a string with the desired health-related state
  # it gives back a tibble with age and new cases of that state at that age.
  # Note that new cases are computed as all transitions going to that state
  # and coming from a different stage
  
  # Columns that containst `state`:
  col_set <- sim_no_trt$Tot_Trans_per_t %>% 
    as_tibble() %>% 
    select(contains(state)) %>%
    colnames() #%>%
  
  # Columns with transitions to `state`: `xx->state`with `xx!=state`
  new_col_set <- NULL
  for (i in 1:length(col_set)) {
    if ((stringr::str_split(string = col_set[i], pattern = "->") %>% 
         unlist() %>% .[2]) == state &&
        (stringr::str_split(string = col_set[i], pattern = "->") %>% 
         unlist() %>% .[1]) != state)
    {
      new_col_set <- append(new_col_set, col_set[i])
    } 
  }
  
  # Select columns with new_cases of `state`:
  new_cases <- sim_no_trt$Tot_Trans_per_t %>%
    as_tibble() %>%
    dplyr::select(one_of(new_col_set)) %>%
    dplyr::mutate(total_new_cases = rowSums(.)) %>%
    dplyr::select(-everything(), total_new_cases) %>%
    #dplyr::mutate(age = row_number() + 9)
    dplyr::mutate(age = row_number() + 10)
  
  return(new_cases)
}
####
new_cases_2 <- function(state1, state2, Tot_Trans_per_t){
  # This function receives two strings with the names of a state as named in the 
  # vector state `v_n` and it gives back a df with two columns: age and number of
  # new transitions in that age.
  transition_cases <- 
    Tot_Trans_per_t %>% 
    as_tibble() %>% 
    dplyr::select(paste0(state1,"->",state2)) %>% 
    dplyr::mutate(age = row_number() +10, 
                  cycle = as.numeric(row_number()))
  return(transition_cases)
}


## ----MicroSim function, tidy=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Mod: incorporate loop over simulations:
# This version stacks solution of simulations but produces a list with stacked elements
MicroSim <- function(strategy="natural_history", numb_of_sims = 1,
                     v_M_1, n_i, n_t, v_n, d_c, d_e, TR_out = TRUE, 
                     TS_out = TRUE, Trt = FALSE,  seed = 1, Pmatrix) 
{
  
  seeds <- sample(1:10000, numb_of_sims, replace = FALSE)  # Generate random seeds
  {
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
    
    for(sim in 1:numb_of_sims) {
      
      cat("Running simulation", sim, "with seed", seeds[sim], "\n")
      
      symptomatics <-
        data.frame(ID = integer(), TimeStep = integer(), 
                   DiagnosedState = character(), 
                   RecoveredFromState = logical(), stringsAsFactors = FALSE)
      
      v_dwc <- 1 / (1 + d_c) ^ (0:(n_t-1))   # calculate the cost discount weight based
                                             # on the discount rate d_c 
      print(length(v_dwc))
      
      v_dwe <- 1 / (1 + d_e) ^ (0:(n_t-1))   # calculate the QALY discount weight based 
                                             # on the discount rate d.e
      
      # Create the matrix capturing the state name/costs/health outcomes 
      # for all individuals at each time point:
      m_M <- m_C <- m_E <-  matrix(nrow = n_i, ncol = (n_t + 1), 
                                   dimnames = list( 1:n_i, 
                                                    paste0("cycle_", 1:(n_t + 1), sep = "")))  
      
      m_M[, 1] <- v_M_1  # indicate the initial health state   
      
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
      
      # new code:
      stored_list <- list()
      ######################## run over all the cycles ############################# 
      for (t in 1:(n_t)) {
        ############################################################################
        # Select the transition matrix based on the cycle `n_t`:
        # Since our age intervals start at 10 years old,
        # `age_in_loop` is calculated as `t + 10 * age_factor(cycle_period)`.
        age_in_loop <- t + 10
        # Choose corresponding transition matrix according current age:
        ## DONT USE it for cycle_period = "1yr"
        #my_age_prob_matrix <- my_Probs %>% 
        #  dplyr::filter(Lower <= 
        #                  (age_in_loop / age_factor(cycle_period)) &
        #                  Larger >= (age_in_loop / age_factor(cycle_period)) %>%
        #                  floor()) 
        my_age_prob_matrix <- my_Probs %>% 
          dplyr::filter(Lower <= age_in_loop  &
                          Larger >= age_in_loop) 
        
        # Add colnames and update `v_n`:
        rownames(my_age_prob_matrix) <- v_n <<- 
          my_age_prob_matrix %>%
          dplyr::select(-c(Age.group, Lower, Larger)) %>% 
          colnames()
        ############################################################################
        
        # update/correct n_s (<<- let change variable from inside a function):
        n_s  <<- length(v_n)  
        
        ############################################################################ 
        #new code:
        new_entries <- diagnose_column(m_M[, t], t)
        
        if (!is.null(new_entries)) {
          stored_list[[t]] <- new_entries
        }
        
        if (nrow(new_entries) > 0) {
          symptomatics <- bind_rows(symptomatics, new_entries)
        }
        ############################################################################ 
        
        
        ############################################################################    
        # Extract the transition probabilities of each individuals at cycle t
        # given the individual current state and the corresponding 
        # transition probability matrix that depends on age:
        ## Debugging:
        #cat("The time step is", t, "\n")
        m_P <- Probs(M_it =  m_M[, t], my_Probs = my_age_prob_matrix)
        
        m_M[, t + 1] <- samplev(probs = m_P, m = 1)  # sample the next health state 
                                                     # and store that state in  
                                                     # matrix m_M 
        ############################################################################    
        
        
        # m_M[, t + 1] <- update_column(m_M[, t], new_entries)
        next_col <- m_M[, t + 1]
        next_col <- update_column(m_M[, t], new_entries, next_col)
        
        # Ensure next_col updates are preserved after sampling
        m_M[, t + 1] <- ifelse(next_col == "Survival", "Survival", m_M[, t + 1])
        
        
        ############################################################################    
        ## Costs per CC diagnose at time t + 1.
        
        # Estimate costs per individual during cycle t + 1 conditional on treatment:
        m_C[, t] <-                              
          Costs_per_Cancer_Diag(M_it = m_M[, t],  
                                symptomatics = symptomatics,
                                time_iteration = t,
                                cost_Vec = cost_Vec,    
                                Trt)            
        
        
        m_E[, t + 1] <- # estimate QALYs per individual during cycle t + 1
          Effs( m_M[, t + 1], Trt, 
                utilityCoefs = utilityCoefs)                   
        ############################################################################    
        
        # Conditional on treatment
        cat('\r', paste(round(t/n_t * 100),          # display the 
                        "% done", sep = " "))        # progress of  the simulation                    
      }  
      ######################## close loop for cycles ############################### 
      
      # Combine stored entries into a single data frame
      symptomatics <- bind_rows(stored_list)
      tc_disc <- m_C[,1:n_t] %*% v_dwc       # total (discounted) cost per individual
      te_disc <- m_E[,1:n_t] %*% v_dwe       # total (discounted) QALYs per individual 
      
      tc_undisc <- m_C[,1:n_t] %*% rep(1, n_t)       # total (discounted) cost per individual
      te_undisc <- m_E[,1:n_t] %*% rep(1, n_t)       # total (discounted) QALYs per individual 
      
      tc_hat_disc <- mean(tc_disc)        # average (discounted) cost 
      te_hat_disc <- mean(te_disc)        # average (discounted) QALYs
      tc_hat_undisc <- mean(tc_undisc)        # average (discounted) cost 
      te_hat_undisc <- mean(te_undisc)        # average (discounted) QALYs
      
      # Create a matrix of transitions across states transitions from one state to the other:
      if (TS_out == TRUE) {  
        TS <- paste(m_M, cbind(m_M[, -1], NA), sep = "->")    
        
        TS <- matrix(TS, nrow = n_i)
        rownames(TS) <- paste("Ind",   1:n_i, sep = " ")   # name the rows 
        colnames(TS) <- paste0("cycle_", 1:(n_t + 1), sep = "")   # name the columns 
      } else {
        TS <- NULL
      }
      
      if (TR_out == TRUE) {
        TR <- t(apply(m_M, 2, 
                      function(x) table(factor(x, levels = v_n, ordered = TRUE))))
        #TR <- TR / n_i                                   # create a distribution 
        # trace
        
        rownames(TR) <- paste("cycle", 1:(n_t + 1), sep = "_") # name the rows 
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
        Tot_Trans_per_t <- t(apply(TS, 2, 
                                   function(x) 
                                     table(factor(x, levels 
                                                  = transitions, 
                                                  ordered = TRUE))))
        # trace
        rownames(Tot_Trans_per_t) <- paste0("cycle_", 1:(n_t + 1), sep = "") # name the rows 
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
      
      
      # Before sending back, some cleaning regarding cycle `n_t+1` which is 
      # computed but no needed as a result:
      m_M <-m_M[, 1:n_t]
      m_C <-m_C[, 1:n_t]
      m_E <-m_E[, 1:n_t]
      new_CIN1 <- new_CIN1 %>% dplyr::slice(c(1:n_t))
      new_CIN2 <- new_CIN2 %>% dplyr::slice(c(1:n_t))
      new_CIN3 <- new_CIN3 %>% dplyr::slice(c(1:n_t))
      new_Cancer <- new_Cancer %>% dplyr::slice(c(1:n_t))
      
      # Removing no needed extra row from TR:
      row_to_remove <- n_t + 1
      TR <- TR[-row_to_remove, ]
      
      rm(row_to_remove)
      
      # Removing extra column no needed in TS
      TS <- TS[, -(n_t + 1)]
      
      ### add age to TR:
      TR <- as.data.frame(TR)
      TR <- TR %>% mutate(age = row_number() + 10)
      
      ## HERE I SHOULD IMPLEMENT THE COMPUTATION OF PREVALENCE AND INCIDENCE (?):
      #prevalence_func <- function(cases, )
      
      #Remove large objects: 
      #rm(m_M, m_C, m_E)
      
      # Store the results from the simulation in a list
      results <- list(strategy = strategy,
                      seed = seeds[sim],
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
                      new_Cancer = new_Cancer)  
      
      
    #results$seed <- seeds[sim]
    simulation_results[sim] <- list(results)
    } # end of `n_t` loop
    
    ## Debugging:
    #print(results)
    
  }  # end of `numb_of_sims` loop
  
  #return(simulation_results)
  
  # stack results
  #source("./R/Sumarize_results_by_Strategy_Func.R")
  stacked_results <- 
    summarize_results_by_Strategy(results_list = simulation_results, 
                                 numb_of_sims = numb_of_sims)
  
  # Debugging:
  
  cat("The class of the stacked results is: ", stacked_results %>% class() %>% print(), "\n")
  
  return(stacked_results)
  
} # end of MicroSim function

## Or should the prevalence/incidence coputed outside the MicroSim function?
## HERE I SHOULD IMPLEMENT THE COMPUTATION OF PREVALENCE AND INCIDENCE (?):
#prevalence_func <- function(cases, )


## ----perform simulation, tidy=TRUE, echo=FALSE, include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------
########################## Run the simulation ##################################
## START SIMULATION
p = Sys.time()
# run for no treatment
#sim_no_trt  <- MicroSim(v_M_1, n_i, n_t, v_n, d_c, d_e, Trt = FALSE)
sim_no_trt  <- MicroSim(strategy = "natural_history",numb_of_sims = 100, 
                        v_M_1 = v_M_1, n_i = n_i, n_t = n_t, v_n = v_n, 
                        d_c = d_c, d_e = d_e, TR_out = TRUE, TS_out = TRUE, 
                        Trt = FALSE, seed = 1, Pmatrix = Pmatrix)

#prevalence_func <- function(sim_stalked_result) {
#  df <- sim_stalked_result[[1]]$TR %>% 
#    dplyr::select(sim, cycle, age, H, HR.HPV.infection) %>% 
#    mutate(prevalence = HR.HPV.infection / H) %>% 
#    group_by(age) %>% 
#    summarise(prevalence = mean(prevalence, na.rm = TRUE)) %>% ungroup()
#  prevalence_per_age <- df
#  #return(c(sim_stalked_result, prevalence_per_age))
#  sim_stalked_result[[1]]$mean_prevalence_per_age <- df
#  return(sim_stalked_result)
#}

################################################################################
# Prevalence is defined as number of infected divided by total alive individuals
# for that cycle/time step
prevalence_func <- function(sim_stalked_result, my_Probs) {
  # Extract unique age intervals
  age_intervals <- my_Probs %>% 
    select(Lower, Larger) %>% 
    unique() %>% 
    arrange(Lower)
  
  # Create a vector of the breaks for the intervals
  breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
  
  # Create labels for the intervals
  labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
  
  # Compute prevalence and average it by age intervals
  df <- sim_stalked_result[[1]]$TR %>% 
    #dplyr::select(sim, cycle, age, H, HR.HPV.infection) %>% 
    dplyr::select(everything()) %>% 
    mutate(total_alive = H + CIN1 + CIN2 +CIN3 +
             FIGO.I + FIGO.II + FIGO.III + FIGO.IV + Survival) %>%
    #mutate(prevalence = HR.HPV.infection / H) %>% 
    mutate(prevalence = HR.HPV.infection / total_alive) %>% 
    mutate(age_interval = cut(age, breaks = breaks, labels = labels, right = FALSE)) %>% 
    group_by(age_interval) %>% 
    summarise(prevalence = mean(prevalence, na.rm = TRUE)) %>% 
    ungroup()
  
  #return(df)
  sim_stalked_result[[1]]$mean_prevalence_per_age_interval <- df
  return(sim_stalked_result)
}
################################################################################

# Concatenate the prevalence to the sim result 
prevalence_result <-
  prevalence_func(sim_stalked_result = sim_no_trt, my_Probs = my_Probs)  

# save the results
saveRDS(object = prevalence_result, file = "./data/stacked_sims_100x10E6x75.rds")

################################################################################
## Incidence is defined by the NEW number of individuals in the state of interest
## in the time t divided by all the individuals in  the epidemiological "precedent" 
## states at time t-1 (in the previous cycle). This is how is defined in the  the 
## Markov model. We follow that definition to compare the micro sim and the Markov. 
#incidence_func <- function(sim_stalked_result, my_Probs) { # work in progress
#  
#  # Extract unique age intervals
#  age_intervals <- my_Probs %>% 
#    select(Lower, Larger) %>% 
#    unique() %>% 
#    arrange(Lower)
#  
#  # Create a vector of the breaks for the intervals
#  breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
#  
#  # Create labels for the intervals
#  labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
#  
#  # Compute prevalence and average it by age intervals
#  df <- sim_stalked_result[[1]]$TR %>% 
#    #dplyr::select(sim, cycle, age, H, HR.HPV.infection) %>% 
#    dplyr::select(everything()) %>% 
#    mutate(total_alive = H + CIN1 + CIN2 +CIN3 +
#             FIGO.I + FIGO.II + FIGO.III + FIGO.IV + Survival) %>%
#    #mutate(prevalence = HR.HPV.infection / H) %>% 
#    mutate(prevalence = HR.HPV.infection / total_alive) %>% 
#    mutate(age_interval = cut(age, breaks = breaks, labels = labels, right = FALSE)) %>% 
#    group_by(age_interval) %>% 
#    summarise(prevalence = mean(prevalence, na.rm = TRUE)) %>% 
#    ungroup()
#  
#  #return(df)
#  sim_stalked_result[[1]]$mean_prevalence_per_age_interval <- df
#  return(sim_stalked_result)
#}

# run for treatment
#sim_trt     <- MicroSim(v_M_1, n_i, n_t, v_n, d_c, d_e, Trt = TRUE)  
comp.time = Sys.time() - p
comp.time %>% print()


## ----convert .Rmd to .R---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#library(knitr)
## purl("your_script.Rmd", output = "your_script.R")
## example:
#purl("Cervix_MicroSim_RMarkdown_v.072.Rmd", output = "cervix_microSim_stacked_list.R")


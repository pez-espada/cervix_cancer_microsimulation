## ----Preamble
################################################################################
# This code is a modified version of the original code from:
# [https://github.com/DARTH-git/Microsimulation-tutorial] (Krijkamp et al 2018 
# Sick-Sicker model).
# programmed by Carlos Dommar D'Lima - carlos.dommar@gmail.com
# This code extends the "sick-sicker" model of the original authors to a
# multi-state cervix cancer model
################################################################################
rm(list = ls())
library(tidyverse)

## to prevent conflicts in the parallel environment:
#setwd(dir = "/home/07075107P/microSim/cervix_cancer_microsimulation")

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

my_Probs <- readRDS(file = "./data/probs.rds") # natural history transition matrix
my_Probs2 <- readRDS(file = "./data/probs2.rds") # vaccination transition matrix
my_Probs4 <- readRDS(file = "./data/probs3.rds") # vaccination transition matrix
my_Probs9 <- readRDS(file = "./data/probs3.rds") # vaccination transition matrix



################################################################################
# Function to extract and convert numbers from factor levels
extract_numbers <- function(range_factor) {
range_string <- as.character(range_factor)
numbers <- as.numeric(unlist(strsplit(range_string, "-")))
return(numbers)
}
################################################################################

################################################################################
my_Probs_cleaning_Func <- function(Probs_matrix) {
  # Tidying up a bit the transition matrix:
  Probs_matrix  <- Probs_matrix %>% dplyr::rename("H" = "Well")
  Probs_matrix  <- Probs_matrix %>% as.data.frame() #convert back to data.frame (no needed?)
  Probs_matrix$Lower  <- sapply(Probs_matrix$Age.group, function(x) extract_numbers(x)[1])
  Probs_matrix$Larger <- sapply(Probs_matrix$Age.group, function(x) extract_numbers(x)[2])
  return(Probs_matrix) 
}
################################################################################


# Tidying up a bit the transition matrix:
my_Probs <- my_Probs %>% as.data.frame() #convert back to data.frame (no needed?)
my_Probs <- my_Probs_cleaning_Func(Probs_matrix = my_Probs)

my_Probs2 <- my_Probs2 %>% as.data.frame() %>% #convert back to data.frame (no needed?)
  dplyr::mutate(Age.group = ifelse(Age.group == "11-14", "10-14", Age.group))
my_Probs2 <- my_Probs_cleaning_Func(Probs_matrix = my_Probs2)

my_Probs4 <- my_Probs4 %>% as.data.frame() %>% #convert back to data.frame (no needed?)
  dplyr::mutate(Age.group = ifelse(Age.group == "11-14", "10-14", Age.group))
my_Probs4 <- my_Probs_cleaning_Func(Probs_matrix = my_Probs4)

my_Probs9 <- my_Probs9 %>% as.data.frame() %>% #convert back to data.frame (no needed?)
  dplyr::mutate(Age.group = ifelse(Age.group == "11-14", "10-14", Age.group))
my_Probs9 <- my_Probs_cleaning_Func(Probs_matrix = my_Probs9)
################################################################################
## ----Model Parameters
n_i <- (2)*10^5         # number of simulated individuals
#n_i <- (5)*10^5            # number of simulated individuals
#n_i <- 10^7            # number of simulated individuals
n_i <- 10^5               # number of simulated individuals
#n_i <- 10^6               # number of simulated individuals
n_t <- 75                  # time horizon, 75 cycles (it starts from 1)
################################################################################
 

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

################################################################################
# Cost and utility inputs 
# From our Markov cervix model (CC's natural history?):
cost_Vec = c(0, 39.54, 288.91, 1552.27, 1552.27, 
           5759.81, 12903.63, 23032.41, 35323.14, 0, 0, 0)
utilityCoefs = c(1, 1, 0.987, 0.87, 0.87, 0.76, 0.67, 0.67, 0.67, 0.938, 0, 0)
################################################################################


################################################################################
## ---- FUNCTIONS -----                                                       ##  
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
transition_prob<-P[state1,state2]
return(transition_prob)
}
################################################################################
# 

################################################################################
## ---- Probability Function ----                                             ##
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
## ---- Probability Function ----                                             ##
## The Probs_2 function that updates the transition probabilities of every cycle:
## taking into account other probs than natura history
## depending on the vaccination startegies
Probs_2 <- function(M_it, my_Probs, my_Probs2, my_Probs4, my_Probs9, vacc_lbl) {
  # M_it: matrix of health states of all individuals at time t
  # my_Probs: list of distinct transition matrices for each vaccination strategy
  # vacc_vector: vector of vaccination strategies for each individual
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
## ----Sampling function
# Efficient implementation of the rMultinom() function of the Hmisc package #### 
# This function samples the next health state of each individual based on the
# transition probabilities of the current health state of each individual.
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
################################################################################

################################################################################
## ---- Costs Function ----                                                   ##
### Costs Function
# The `Costs_per_Cancer_Diag` function estimates the costs of a diagnose 
# individual due to cancer symptoms (FIGO.I-IV) at every cycle. 
# This cost is only charged once in the patient's lifetime.
# NOTE: need to decide if the cost is applied on current time `t` or `t+1` as it is now.
Costs_per_Cancer_Diag <- function (M_it, cost_Vec, symptomatics, time_iteration, Trt = FALSE) {
  c_it <- rep(0, length(M_it))
  #ci_t <- 0
  if(nrow(symptomatics) > 0 ) {
    c_it[symptomatics %>% 
           dplyr::filter(DiagnosedState == "FIGO.I" & TimeStep == time_iteration) %>% 
           select(ID) %>% as.list() %>% 
           unlist()] <- cost_Vec[which(v_n %in% "FIGO.I")]
    c_it[symptomatics %>% 
           dplyr::filter(DiagnosedState == "FIGO.II" & TimeStep == time_iteration) %>% 
           select(ID) %>% as.list() %>% 
           unlist()] <- cost_Vec[which(v_n %in% "FIGO.II")]
    c_it[symptomatics %>%
           dplyr::filter(DiagnosedState == "FIGO.III" & TimeStep == time_iteration) %>% 
           select(ID) %>% as.list() %>% 
           unlist()] <- cost_Vec[which(v_n %in% "FIGO.III")]
    c_it[symptomatics %>% 
           dplyr::filter(DiagnosedState == "FIGO.IV" & TimeStep == time_iteration) %>% 
           select(ID) %>% as.list() %>% 
           unlist()] <- cost_Vec[which(v_n %in% "FIGO.IV")]
  }
  return(c_it) # return the costs
}
################################################################################


################################################################################
## ---- Qalys Function ----                                                   ##
### Health outcome function 
# The `Effs` function estimates the QALYs of a diagnose individual due to cancer
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
  # If the TryCatch gives proble, just overrate it:
  #for (i in 1:length(utilityCoefs)) {
  #  u_it[M_it == v_n[i]] <- utilityCoefs[i]   # update the utility if healthy
  #}
  return(u_it)
}
################################################################################


#################################################################################
### ----Time period related functions
############ WORK IN PROGRESS #########################
#age_factor <- function(my_period) {
## it receives a string with the period of the cycle, and it can be:
##  - "1mth"
##  - "3mth"
##  - "4mth"
##  - "6mth"
##  - "1yr" # i.e. 12 months
## and it gives back an age factor for scaling cycle period.
#if (my_period == "1yr") {
#  my_factor <- 1
#} else if (my_period == "6mth") {
#  my_factor <- 2
#} else if (my_period == "4mth") {
#  my_factor <- 3
#} else if (my_period == "3mth") {
#  my_factor <- 4
#} else if (my_period == "1mth") {
#  my_factor <- 12
#} else {print("Cycle period can only be: '1yr', '6mth','4mth', '3mth' and '12mth'")}
#return(my_factor)
#}
########## WORK IN PROGRESS #################
#################################################################################


#################################################################################
##### ! NOT USED ! ############################
#convert_matrix_to_proper_transition <- 
#function(my_age_prob_matrix, cycle_period) {
#  my_age_prob_matrix %>% head(3)
#  ensure_library(c("expm", "pracma", "ctmcd"))
#  trans_matrix <- my_age_prob_matrix %>% 
#    select(-c("Age.group", "Lower", "Larger")) %>% 
#    as.matrix()
#  # Referenece: https://rpubs.com/crossxwill/transition_matrix
#  ## method 1: (not working atm)
#  #ensure_library(expm)
#  #TM.exp  <- expm::expm((1 / age_factor(cycle_period))) * log(trans_matrix) 
#  
#  #method 2 ;
#  #ensure_library("pracma")
#  TM_pracma <- 
#    pracma::rootm(trans_matrix, p=age_factor(cycle_period), 
#                  kmax = 20, tol = 1e-10)
#  round(TM_pracma$B, 5)
#  # Regularization with the `ctmcd` package, The code below uses the 
#  # quasi-optimization of the generator (QOG) approach from 
#  # Kreinin and Sidelnikova (2001).:
#  ensure_library("ctmcd")
#  TM_qo <- ctmcd::gm(TM_pracma$B, te=1, method = "QO") 
#}
##### ! NOT USED ! ############################
#################################################################################


## ---- Symptomatic Individuals ----                                                         ##
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
#global_diagnosed <- integer()


################################################################################
# --- Function receives a column with current state of `n_i` individuals and gives
# a dataframe with `ID, TimeStep`, `state`, and `RecoveredFromState` columns.
# The function also updates the global vector `global_diagnosed` with the IDs of
# individuals who have been diagnosed.
diagnose_column <- function(col, time_step) {
  new_entries <- data.frame(ID = integer(), 
                            TimeStep = integer(),
                            DiagnosedState = character(),
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


################################################################################
# ---- Function to update the next column based on the new entries ----       ##
# This function updates the next column based on the new entries of diagnosed
# individuals. It also updates the state of individuals who have recovered.
# The function returns the updated next column.
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
# ---- Function to add new cases to the transition matrix ----                ##
# This function adds new rows to the transition matrix for individuals who have
# been diagnosed with cancer. It also updates the age and cycle columns.
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
      
      # Modify the dataframe: Add a new row with age = 10 and
      # transition column = 0, and delete the last row (age = 85)
      transition_cases <- transition_cases %>%
        # Add row at the beginning
        add_row(!!transition_column := 0, age = 10, cycle = 1, .before = 1) %>%  
        slice(-n()) %>%  
        # Remove the last row
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
      !!!setNames(lapply(missing_cols, 
                         function(x) rep(0, nrow(Tot_Trans_per_t_tbl))), 
                  missing_cols)
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
      mutate(age = 10:(10 + n() - 1), cycle = age - 9) #%>%
    #dplyr::select(-matches("sim\\.x$")) %>%
    #dplyr::select(-sim.1)
    
    return(transition_cases)
  }
}
################################################################################


################################################################################
# Function to update the transition matrix with new cases
my_age_prob_matrix_func <- function(my_Prob_matrix, my_age_in_loop) {
  my_age_prob_matrix <- my_Prob_matrix %>% 
    dplyr::filter(Lower <= my_age_in_loop  &
                    Larger >= my_age_in_loop) 
}
################################################################################

 
################################################################################
## Function to update the transition matrix with new cases
## This function will take into account probabilities depending on 
## vaccination strategies
#my_age_prob_matrix_func_2 <- function(my_Prob_matrix, my_age_in_loop) {
#  my_age_prob_matrix <- my_Prob_matrix %>% 
#    dplyr::filter(Lower <= my_age_in_loop  &
#                    Larger >= my_age_in_loop) 
#}
################################################################################


################################################################################
## THE MICROSIMULATION MAIN FUNCTION
# This version stacks solution of simulations but produces a list with stacked elements
# check the `MicroSim` for any improvements or issues.
MicroSim <- function(strategy="natural_history", numb_of_sims = 20,
                     v_M_1, n_i, n_t, v_n, d_c, d_e, TR_out = TRUE, 
                     TS_out = TRUE, Trt = FALSE,  seed = 1, Pmatrix, vaccination = FALSE) 
{
  # Generate random seeds
  #seeds <- sample(1:10000, numb_of_sims, replace = FALSE)  
  seeds <- sample(1:100000, numb_of_sims, replace = FALSE)  
  ## fix the seeds for reproducibility::
  #seeds <- c(38222, 52130, 92742, 73352, 41494, 43929, 94560, 72382, 13846, 94537) %>% 
  #  as.integer()
  #seeds <- c(20422, 63139, 3575,  9449,  4055,  
  #           6931, 92384, 24048, 25109,  7757,
  #           25889, 32227, 57572, 36484, 38944,  
  #           4074, 45156, 93585, 48543, 57217) %>%
  #  as.integer()
  
  simulation_results <- list() 
  
  #calculate the cost discount weight based on the discount rate d_c 
  v_dwc <- 1 / (1 + d_c) ^ (0:(n_t-1))   
  # calculate the QALY discount weight based on the discount rate d_e                                             
  v_dwe <- 1 / (1 + d_e) ^ (0:(n_t-1))   
  
  # Parallel processing using foreach
  simulation_results <- 
    foreach(sim = 1:numb_of_sims, .packages = c("dplyr", "tidyr", "purrr") ) %dopar% { 
      ## clean memory:
      #if (step %% 10 == 0) gc()
      cat("Running simulation", sim, "with seed", seeds[sim], "\n")
      # Initialize a global vector to store all diagnosed individuals
      global_diagnosed <<- integer()
      symptomatics <-
        data.frame(ID = integer(), TimeStep = integer(), 
                   DiagnosedState = character(), 
                   RecoveredFromState = logical(), stringsAsFactors = FALSE)
      
      ## NOTA: PONER FUERA DEL LOOP (??)
      ##calculate the cost discount weight based on the discount rate d_c 
      #v_dwc <- 1 / (1 + d_c) ^ (0:(n_t-1))   
      ## calculate the QALY discount weight based on the discount rate d_e                                             
      #v_dwe <- 1 / (1 + d_e) ^ (0:(n_t-1))   
      
      # Create the matrix capturing the state name/costs/health outcomes 
      # for all individuals at each time point:
      #m_M <- m_C <- m_E <-  matrix(nrow = n_i, ncol = (n_t + 1), 
      m_M <- m_C <- m_E <- 
        matrix(nrow = n_i, ncol = (n_t), 
               dimnames = list( 1:n_i, 
                                paste0("cycle_", 1:(n_t),
                                       sep = "")))  
      
      m_M[, 1] <- v_M_1  # indicate the initial health state   
      
      seed <- seeds[sim]
      #seed <- 17
      set.seed(seed) # set the seed for every individual 
      
      # estimate costs per individual for the initial health state
      m_C[, 1] <- Costs_per_Cancer_Diag(M_it = m_M[, 1], 
                                        symptomatics = symptomatics,
                                        time_iteration = 1,
                                        cost_Vec = cost_Vec,  
                                        Trt)             
      
      # estimate QALYs per individual for the initial health state 
      m_E[, 1] <- Effs(m_M[, 1], Trt, utilityCoefs = utilityCoefs)  
      
      stored_list <- list()
       
      ###################### run over all the cycles ########################### 
      # loop runs over all the cycles of the simulation. It updates the
      # health state of each individual at each cycle, estimates the costs and
      # QALYs per individual at each cycle, and stores the transitions across
      # states for each individual at each cycle.
      for (t in 1:(n_t-1)) {
        ########################################################################
        # Select the transition matrix based on the cycle `n_t`:
        # Since our age intervals start at 10 years old,
        age_in_loop <- t + 9
        ########################################################################
        
        # update/correct n_s (<<- let change variable from inside a function):
        # q: there is a bad practice to use <<- in a function?
        # a: Yes, it is a bad practice to use <<- in a function.
        # q: how can I avoid it in this case?
        # a: You can avoid it by passing the variable as an argument to the function.
        #n_s  <<- length(v_n)  
        
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
        
        ## Here I need to modify the following function to extract the the right
        ## transition matrix based on the age of the individual at each cycle, and
        ## the correponding transition matrix that depends on vaccination strategies
        my_age_prob_matrix <- 
          my_age_prob_matrix_func(my_Prob_matrix = my_Probs, 
                                  my_age_in_loop = (age_in_loop + 1))
        # Add colnames and update `v_n`:
        rownames(my_age_prob_matrix) <- v_n <<- 
          my_age_prob_matrix %>%
          dplyr::select(-c(Age.group, Lower, Larger)) %>% 
          colnames()
       
         
        my_age_prob_matrix_2 <- 
          my_age_prob_matrix_func(my_Prob_matrix = my_Probs2, 
                                  my_age_in_loop = (age_in_loop + 1))
        # Add colnames and update `v_n`:
        rownames(my_age_prob_matrix_2) <- v_n <<- 
          my_age_prob_matrix_2 %>%
          dplyr::select(-c(Age.group, Lower, Larger)) %>% 
          colnames()
        
        
        my_age_prob_matrix_4 <- 
          my_age_prob_matrix_func(my_Prob_matrix = my_Probs4, 
                                  my_age_in_loop = (age_in_loop + 1))
        # Add colnames and update `v_n`:
        rownames(my_age_prob_matrix_4) <- v_n <<- 
          my_age_prob_matrix_4 %>%
          dplyr::select(-c(Age.group, Lower, Larger)) %>% 
          colnames()
        
        
        my_age_prob_matrix_9 <- 
          my_age_prob_matrix_func(my_Prob_matrix = my_Probs9, 
                                  my_age_in_loop = (age_in_loop + 1))
        # Add colnames and update `v_n`:
        rownames(my_age_prob_matrix_9) <- v_n <<- 
          my_age_prob_matrix_9 %>%
          dplyr::select(-c(Age.group, Lower, Larger)) %>% 
          colnames()
        
        
          # Extract the transition probabilities of each individuals at cycle t
        # given the individual current state and the corresponding 
        # transition probability matrix that depends on age:
        # Next time (t+1) transition
        # m_P is a (n_i x n_s) matrix with the probabilities of transitioning
        m_P <- Probs(M_it =  m_M[, t], my_Probs = my_age_prob_matrix)
        
        # for vaccination I'll need a new Probs function: 
        #m_P <- Probs_2(M_it = m_M[, t], my_Probs = c(), vacc_vector = vaccination)
        
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
        m_C[, t + 1] <-                              
          Costs_per_Cancer_Diag(M_it = m_M[, t + 1],  
                                symptomatics = symptomatics,
                                time_iteration = t,
                                cost_Vec = cost_Vec,    
                                Trt) %>% round(., 4)            
        #browser()
        
        m_E[, t + 1] <- # estimate QALYs per individual during cycle t + 1
          Effs( m_M[, t + 1], Trt, 
                utilityCoefs = utilityCoefs)                   
        ########################################################################    
        cat('\r', paste(round(t/n_t * 100),          # display the 
                        "% done\n", sep = " "))        # progress of  the simulation                    
        
      }  
      #################### close loop for cycles ############################### 
      
      # Combine stored entries in a single data frame
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
        colnames(TS) <- paste0("cycle_", 1:(n_t), sep = "")   # name the columns 
      } else {
        TS <- NULL
      }
      
      if (TR_out == TRUE) {
        TR <- t(apply(m_M, 2, 
                      function(x) table(factor(x, levels = v_n, ordered = TRUE))))
        #TR <- TR / n_i                                   # create a distribution 
        # trace
        
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
      
      new_CC_Death <- new_cases_2(state1 = c("H", "HR.HPV.infection", "CIN1", 
                                             "CIN2","CIN3","FIGO.I", 
                                             "FIGO.II", "FIGO.III", "FIGO.IV", 
                                             "Survival"),
                                  state2 = "CC_Death", 
                                  Tot_Trans_per_t = Tot_Trans_per_t)
      
      
      # Before sending back, some cleaning regarding cycle `n_t+1` which is 
      # computed but no needed as a result:
      m_M <- m_M[ , 1:n_t]
      m_C <- m_C[ , 1:n_t]
      m_E <- m_E[ , 1:n_t]
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
      TR <- TR %>% mutate(age = row_number() + 9)
      TR$sim <- sim
      
      #Remove large objects: 
      #rm(m_M, m_C, m_E)
      #rm(m_M, m_C, m_E, TS,tc_disc,tc_undisc,te_disc,te_undisc)
      
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
      
      cat("At sim number:", sim,  " reported strategy is ", strategy, "\n")
      
      # Store the results from the simulation in a list
     results <- list(#strategy = strategy,
                      #seed = seeds[sim],
                      seed = seed,
                      #sim_numb = sim, 
                      m_M = m_M, 
                      #m_C = m_C, 
                      #m_E = m_E, 
                      #tc_disc = tc_disc, 
                      #tc_undisc = tc_undisc,
                      #te_disc = te_disc,
                      #te_undisc = te_undisc,
                      tc_hat_disc = tc_hat_disc,
                      tc_hat_undisc = tc_hat_undisc,
                      te_hat_disc = te_hat_disc, 
                      te_hat_undisc = te_hat_undisc, 
                      #TS = TS,
                      TR = TR, 
                      #Tot_Trans_per_t = Tot_Trans_per_t, 
                      symptomatics = symptomatics,
                      new_CIN1 = new_CIN1,
                      new_CIN2 = new_CIN2,
                      new_CIN3 = new_CIN3,
                      new_Cancer = new_Cancer,
                      new_CC_Death = new_CC_Death,
                      CC_Death_by_diff = CC_Death_by_diff)  
      results$seed <- seeds[sim]
      #simulation_results[sim] <- list(results)
      #simulation_results[sim] <- results
      cat("At sim number:", sim,  " tc_hat_undisc is ", tc_hat_undisc, "\n")
      rm(symptomatics)
      #rm(TS) 
      return(results)
      #gc() #Force memory cleanup after each sim/batch 
      
    } # end of `foreach/dopar` loop
  
  #return(simulation_results)
  # stack results
  #source("./R/Sumarize_results_by_Strategy_Func.R")
  #source("/home/07075107P/microSim/cervix_cancer_microsimulation/R/Sumarize_results_by_Strategy_Func.R")
  #stacked_results <- 
  #  summarize_results_by_Strategy(results_list = simulation_results, 
  #                                numb_of_sims = numb_of_sims)
  
  #stopCluster(cl)  # Stop the cluster when done
  #return(stacked_results)
  return(simulation_results)
} # end of MicroSim function
################################################################################


################################################################################
                  ###################################
                  ## Pre-simulation Computations: ##
                  ###################################
################################################################################
### Prepare Parallelize code ###
library(parallel)
ensure_library("doParallel")

################################################################################
# Function to detect if running on SLURM -NOT WORKING AS INTENDED"-
is_slurm <- function() {
  slurm_id <- Sys.getenv("SLURM_JOB_ID")
  return(nzchar(slurm_id))  # Returns TRUE only if SLURM_JOB_ID is a non-empty string
}
################################################################################
 
################################################################################
# Determine number of cores
if (is_slurm()) {
  # In Slurm, use the cores requested by the job
  n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
  cat("I'm in slurm!\n")
} else {
  cat("I'm NOT in slurm!\n")
  # On local machine, use all available cores (or limit if needed)
  #n_cores <- parallel::detectCores() - 1  # Use one less than total to avoid overloading
  ## Register fewer cores (adjust based on server resources)
  n_cores <- min(detectCores() - 1, 20)  # Try using 8 or fewer cores
  #n_cores <- detectCores()  # Try using 8 or fewer cores
  #n_cores <- min(detectCores())  # Try using 8 or fewer cores
  #n_cores <- 6  # Try using 8 or fewer cores
}
# for 250000 individuals x 75 cycles x 20 sims in a Lenovo 16GB Laptop use
# five cores. It takes ca 3.5-3.7 minutes to run. Using 7 cores can run the same set
# in 3.3-3.4 minutes but the system becomes unstable and leading to crash often.
# in the office desktop with 3 cores it takes 12.1434, that's roughly 3.6 times slower
#n_cores <- 5 # for personal Lenovo .
cat("Number of cores: ", n_cores, "\n")
################################################################################
 

################################################################################
# 6-hours timeout to prevent socket drop issues
cl <- makeCluster(n_cores, timeout = 6*60*60) 
clusterExport(cl, c("Costs_per_Cancer_Diag", "Effs", "trans_prb", "Probs",
                    "my_Probs", "utilityCoefs", "v_n", "samplev", 
                    "my_age_prob_matrix_func","diagnose_column", 
                    "update_column", "states_to_check", "symptom_prob_vec",
                    "survival_prob_vec", #"global_diagnosed", 
                    "cost_Vec", "new_cases_2"))
#registerDoParallel(cl) # for parallel
registerDoSEQ()        # for sequential
################################################################################
################################################################################


################################################################################
################################################################################
## Vaccination strategies:
## 1. No vaccination
#vaccination <- TRUE
#vaccination <- FALSE
vacc2 <- FALSE
vacc4 <- FALSE
vacc9 <- FALSE

# paramters:
vacc_coverage <- c(0.3, 0.0, 0.0) # vaccination coverage for vacc 2, 4 and 9


generate_vaccine_labels <- function(n_i, vacc_coverage) {
  # Ensure the sum of coverage is valid
  if (sum(vacc_coverage) > 1) {
    stop("The sum of vacc_coverage cannot exceed 1.")
  }
  
  # Calculate the number of individuals for each vaccine
  n_vacc_2 <- round(vacc_coverage[1] * n_i)
  n_vacc_4 <- round(vacc_coverage[2] * n_i)
  n_vacc_9 <- round(vacc_coverage[3] * n_i)
  
  # Remaining individuals are "no_vacc"
  n_no_vacc <- n_i - (n_vacc_2 + n_vacc_4 + n_vacc_9)
  
  if (n_no_vacc < 0) {
    stop("The specified coverage results in more vaccinated individuals than n_i.")
  }
  
  # Create the label vector
  vacc_lbl <- c(
    rep("vacc_2", n_vacc_2),
    rep("vacc_4", n_vacc_4),
    rep("vacc_9", n_vacc_9),
    rep("no_vacc", n_no_vacc)
  )
  
  # Shuffle the vector randomly
  vacc_lbl <- sample(vacc_lbl, size = n_i, replace = FALSE)
  
  return(vacc_lbl)
}
## Example usage
#set.seed(123) # For reproducibility
#n_i <- 1000
#vacc_coverage <- c(0.6, 0.2, 0.2)
#vacc_coverage <- c(0.0, 0.0, 0.0)
#vacc_coverage <- c(0.0, 0.7, 0.0)
#vacc_coverage <- c(0.3, 0.7, 0.1)
vacc_lbl <- generate_vaccine_labels(n_i, vacc_coverage)
#
## Check the results
#table(vacc_lbl) / n_i
################################################################################



################################################################################



################################################################################
########################## Run the simulation ##################################
## START SIMULATION
Sys.setenv(OMP_NUM_THREADS = "1") # to prevent conflicts between OpenMP and R parallel
p = Sys.time()
# run for no treatment
numb_of_sims = 3
strategy <- "natural_history"
sim_no_trt  <- MicroSim(strategy = strategy, numb_of_sims = numb_of_sims, 
                        v_M_1 = v_M_1, n_i = n_i, n_t = n_t, v_n = v_n, 
                        d_c = d_c, d_e = d_e, TR_out = TRUE, TS_out = TRUE, 
                        Trt = FALSE, seed = 2, Pmatrix = Pmatrix)
# Stop the cluster when done
stopCluster(cl)  

# For stacking outside the function, we need to comment the stacking function
# inside  de the MicroSim function, and return the results as a list by commenting
# 'return(stacked_results)' and uncomment 'return(simulation_results)'. And then,
# uncomment the following lines:
#source("./R/sumarize_results_by_Strategy_Func_revised.R")
source("./R/sumarize_results_by_Strategy_Func.R")
stacked_results <- 
  summarize_results_by_Strategy(results_list = sim_no_trt, 
                                numb_of_sims = numb_of_sims)
sim_no_trt <- stacked_results

# Load computed simulation if needed here:
#sim_no_trt <- readRDS(file = "./data/stacked_sims_100x10E6x75.rds")

comp.time = Sys.time() - p
comp.time %>% print()

# adding runtime execution time:
runtime <- comp.time %>% as_tibble() %>% `colnames<-`("runtime")
sim_no_trt[[1]]$runtime <- runtime
sim_no_trt[[1]]$strategy <- strategy
sim_no_trt[[1]]$numb_of_sims   <- numb_of_sims
sim_no_trt[[1]]$numb_of_ind    <- n_i
sim_no_trt[[1]]$numb_of_cycles <- n_t
sim_no_trt[[1]]$seed <- sim_no_trt[[1]]$seed %>% 
  dplyr::select(-c("seed", "row_names")) %>% 
  dplyr::rename("seed" = "sim[[i]][[name_level_of_sim]]")

sim_no_trt[[1]]$tc_hat_undisc <- sim_no_trt[[1]]$tc_hat_undisc %>%
  dplyr::select(-c(tc_hat_undisc)) %>% 
  dplyr::rename("tc_hat_undisc" = "sim[[i]][[name_level_of_sim]]")

sim_no_trt[[1]]$tc_hat_disc <- sim_no_trt[[1]]$tc_hat_disc %>%
  dplyr::select(-c(tc_hat_disc)) %>% 
  dplyr::rename("tc_hat_disc" = "sim[[i]][[name_level_of_sim]]")

sim_no_trt[[1]]$te_hat_undisc <- sim_no_trt[[1]]$te_hat_undisc %>%
  dplyr::select(-c(te_hat_undisc)) %>% 
  dplyr::rename("te_hat_undisc" = "sim[[i]][[name_level_of_sim]]")

sim_no_trt[[1]]$te_hat_disc <- sim_no_trt[[1]]$te_hat_disc %>%
  dplyr::select(-c(te_hat_disc)) %>% 
  dplyr::rename("te_hat_undisc" = "sim[[i]][[name_level_of_sim]]")

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
################################################################################


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
  
  age_intervals <- my_Probs %>% 
    dplyr::select(Lower, Larger) %>% 
    unique() %>% 
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
  
  #return(df)
  sim_stalked_result[[1]]$CC_by_diff_mean_mortality <- df
  return(sim_stalked_result)
}
################################################################################


### TESTING ###
mean_CC_mortality_by_diff_result <- mean_CC_mortality_result
mean_CC_mortality_by_diff_result <-
  mean_CC_mortality_by_diff_func(sim_stalked_result =
                                   mean_CC_mortality_by_diff_result,
                                 my_Probs = my_Probs)  
 
################################################################################
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



################################################################################
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
  df <- sim_stalked_result[[1]]$TR %>% 
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
  
  # Storing the results in the simulation object
  sim_stalked_result[[1]]$mean_FIGO_prevalence <- df
  return(sim_stalked_result) 
}
################################################################################


# Initialize the result with the original structure
sim_result <-  other_mean_mortality_result 
## Concatenate the prevalence to the sim result 
#other_mean_mortality_result <-
#  other_mean_mortality_func(sim_stalked_result = 
#                              other_mean_mortality_result, my_Probs = my_Probs)  

# Concatenate the prevalence to the sim result 
sim_result <-
mean_FIGO_prevalence_Func(sim_stalked_result = 
                              sim_result, my_Probs = my_Probs)  
################################################################################



################################################################################
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
  df <- sim_stalked_result[[1]]$TR %>% 
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
  
  # Storing the results in the simulation object
  sim_stalked_result[[1]]$mean_FIGO <- df
  return(sim_stalked_result) 
}
################################################################################

# Concatenate the prevalence to the sim result 
sim_result <-
  mean_Figo_Func(sim_stalked_result = 
                   sim_result, my_Probs = my_Probs)  
################################################################################



################################################################################
# Mean diagnosed of Cancer averaged by age intervals (FIGO.I-.IV) and by sims
mean_Diagnosed_Func  <- function (sim_stacked_result, my_Probs) {
  age_intervals <- my_Probs %>% 
    dplyr::select(Lower, Larger) %>% 
    unique() %>% 
    arrange(Lower)
  
  # Create a vector of the breaks for the intervals
  breaks <- c(age_intervals$Lower, max(age_intervals$Larger) + 1)
  
  # Create labels for the intervals
  labels <- paste(age_intervals$Lower, age_intervals$Larger, sep = "-")
  
  # extracting diagnosed:
  sympt <- sim_stacked_result[[1]]$symptomatics
  
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
  
  # Storing the results in the simulation object
  sim_stacked_result[[1]]$mean_Diagnosed <- df 
  return(sim_stacked_result) 
}
################################################################################

# Concatenate the prevalence to the sim result 
sim_result <-
  mean_Diagnosed_Func(sim_stacked_result = sim_result, my_Probs = my_Probs)  
################################################################################



################################################################################
################################################################################
## cleaning
source("./R/Remove_columnS_from_list_Func.R")
sim_result <- 
  remove_columns_from_list(complex_list = sim_result, 
                           ... = "sim.1", "row_names")
################################################################################
################################################################################


cat("Hey, I'm done, and about to write out the results\n")

# Save the results to a file
# Get SLURM job ID from the environment variable
slurm_job_id <- Sys.getenv("SLURM_JOB_ID", unset = NA)
if (is.na(slurm_job_id)) {
  slurm_job_id <- format(Sys.time(), "%Y%m%d%H%M%S")  # Fallback to timestamp if not running in SLURM
}
cat("SLURM job ID:", slurm_job_id, "\n")

# Save simulation result:
## Use job ID in file name
#output_file <-
#  paste0("data/testing_stability/stacked_sims_20x10E6x75_20250116_madeinPADO_PARA_from_script_stackedOutside_RND_CORRECTED", slurm_job_id, ".rds")
#saveRDS(object = sim_result, file = output_file)

cat("I have written out the results\n")

### ----Convert .Rmd to .R
#library(knitr)
## purl("your_script.Rmd", output = "your_script.R")
## example:
#purl("Cervix_MicroSim_RMarkdown_v.072_B.Rmd", output = "cervix_microSim_stacked_list.R")
#purl("Cervix_MicroSim_RMarkdown_v.072_B.Rmd", output = "cervix_microSim_stacked_list_B.R")


## ----Cost-Efectivenes
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


## ----Plot curves
## This R chunk is a plot routine (not part of the main program):
library(RColorBrewer)
#ensure_library("RColorBrewer")
# Convert matrix to data frame
#micro_sim_df <- sim_no_trt[[1]]$TR
#micro_sim_df <- other_mean_mortality_result[[1]]$TR
micro_sim_df <- sim_result[[1]]$TR

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


## ----Loading Markov result
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

## Plot the data
#ggplot(long_merged_data, aes(x = age, y = value, color = `Health state`)) +
#  geom_line(linewidth=1, alpha=0.7) +
#  labs(x = "Age", y = "Value", color = "Health state") +
#  ggtitle(expression(paste("Markov cohort vs  Microsimulation for ", 10^6, " individuals"))) + 
#  theme_minimal()  # Optional: customize the theme

################################################################################


## ----Incidences, Prevalences, and Mortalities
# Markov:
markov_CN1_incidences <- c(0.00000, 204.73492, 981.96179, 1368.24200, 3006.85782, 33.48096, 1362.96678, 459.48051, 697.84223, 794.33833, 223.00222, 246.23082, 176.02167, 126.22963, 53.70939)
markov_CN2_incidences <- c(0.000000, 6.165629, 54.767952, 140.309815, 216.568392, 1476.306267, 1579.728160, 1298.914564, 466.596151, 637.661611, 442.298632, 304.784447, 250.953880, 165.628020, 116.925192)
markov_CN3_incidences <- c(0.000000, 2.090325, 9.597415, 44.467676, 148.972191, 0.000000, 3.550684, 91.881726, 12.505042, 68.377446, 25.802481, 7.952667, 1.174088, 1.177840, 2.638642)
markov_CC_incidences  <- c(0.000000, 0.000000, 0.000000, 5.520938, 8.360544, 13.282380, 22.906871, 20.825560, 15.867891, 32.483846, 8.962389, 17.681771, 11.737615, 17.354646, 14.582775)
markov_HPV_prevalences <- c(0.000000000, 0.343480414, 0.377634762, 0.087223460, 0.307341403, 0.030196332, 0.050562845, 0.050151668, 0.082952596, 0.046644059, 0.018532077, 0.034193076, 0.016407832, 0.015039027, 0.003217326)
markov_CC_mortality <- c(0.000000e+00, 0.000000e+00, 0.000000e+00, 2.977975e-06, 1.574920e-05, 2.715056e-05, 5.489929e-05, 7.284815e-05, 1.057494e-04, 5.076268e-05, 7.517773e-05, 4.960943e-05, 4.802468e-05, 4.210457e-05, 4.837655e-05) * 10^5

## MicroSim:
#microSim_CN1_incidences          <- other_mean_mortality_result[[1]]$mean_incidence_CIN1_per_age_interval
#microSim_CN2_incidences          <- other_mean_mortality_result[[1]]$mean_incidence_CIN2_per_age_interval
#microSim_CN3_incidences          <- other_mean_mortality_result[[1]]$mean_incidence_CIN3_per_age_interval
#microSim_CC_incidences           <- other_mean_mortality_result[[1]]$mean_CC_incidence
#microSim_HPV_prevalences         <- other_mean_mortality_result[[1]]$mean_HPV_prevalence_per_age_interval
#microSim_CC_mortality            <- other_mean_mortality_result[[1]]$CC_mean_mortality
#microSim_CC_by_diff_mortality    <- other_mean_mortality_result[[1]]$CC_by_diff_mean_mortality

microSim_CN1_incidences          <-sim_result[[1]]$mean_incidence_CIN1_per_age_interval
microSim_CN2_incidences          <-sim_result[[1]]$mean_incidence_CIN2_per_age_interval
microSim_CN3_incidences          <-sim_result[[1]]$mean_incidence_CIN3_per_age_interval
microSim_CC_incidences           <-sim_result[[1]]$mean_CC_incidence
microSim_HPV_prevalences         <-sim_result[[1]]$mean_HPV_prevalence_per_age_interval
microSim_CC_mortality            <-sim_result[[1]]$CC_mean_mortality
microSim_CC_by_diff_mortality    <-sim_result[[1]]$CC_by_diff_mean_mortality


## ----Ploting incidences and prevalences
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
  #microSim_CN1_incidences = as.numeric(other_mean_mortality_result[[1]]$mean_incidence_CIN1_per_age_interval$mean_incidence_CIN1),
  #microSim_CN2_incidences = as.numeric(other_mean_mortality_result[[1]]$mean_incidence_CIN2_per_age_interval$mean_incidence_CIN2),
  #microSim_CN3_incidences = as.numeric(other_mean_mortality_result[[1]]$mean_incidence_CIN3_per_age_interval$mean_incidence_CIN3),
  #microSim_CC_incidences = as.numeric(other_mean_mortality_result[[1]]$mean_CC_incidence$CC_mean_incidence),
  #microSim_HPV_prevalences = as.numeric(other_mean_mortality_result[[1]]$mean_HPV_prevalence_per_age_interval$prevalence),
  #microSim_CC_mortality = as.numeric(other_mean_mortality_result[[1]]$CC_mean_mortality$CC_mean_mortality),
  #microSim_CC_by_diff_mortality = as.numeric(other_mean_mortality_result[[1]]$CC_by_diff_mean_mortality$CC_by_diff_mean_mortality)
  
  microSim_CN1_incidences  = as.numeric(sim_result[[1]]$mean_incidence_CIN1_per_age_interval$mean_incidence_CIN1),
  microSim_CN2_incidences  = as.numeric(sim_result[[1]]$mean_incidence_CIN2_per_age_interval$mean_incidence_CIN2),
  microSim_CN3_incidences  = as.numeric(sim_result[[1]]$mean_incidence_CIN3_per_age_interval$mean_incidence_CIN3),
  microSim_CC_incidences   = as.numeric(sim_result[[1]]$mean_CC_incidence$CC_mean_incidence),
  microSim_HPV_prevalences = as.numeric(sim_result[[1]]$mean_HPV_prevalence_per_age_interval$prevalence),
  microSim_CC_mortality    = as.numeric(sim_result[[1]]$CC_mean_mortality$CC_mean_mortality),
  microSim_CC_by_diff_mortality = as.numeric(sim_result[[1]]$CC_by_diff_mean_mortality$CC_by_diff_mean_mortality)
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

plot_comparison <- function(data, measure_name) {
  ggplot(data %>% dplyr::filter(grepl(measure_name, measure)), 
         aes(x = age, y = value, fill = model)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = paste(
        measure_name, "Comparison\n",
        "N =", sim_result[[1]]$numb_of_ind, 
        ";  cycles=", sim_result[[1]]$numb_of_cycles, 
        "Parallelized\n", 
        "Numb. of avg. sims =", numb_of_sims
      ),
      x = "Age Group",
      y = measure_name
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 15, hjust = 0.5),  # Adjust size and center the title
      plot.margin = margin(15, 5, 5, 5)                   # Add extra margin
    )
}


################################################################################
## Plotting FIGO prevalences
figo_data_prevalence <- sim_result[["No Intervention"]]$mean_FIGO_prevalence

# Reshape the data into a long format
data_long <- tidyr::pivot_longer(
  figo_data_prevalence,
  cols = starts_with("mean_FIGO"),
  names_to = "FIGO_stage",
  values_to = "prevalence"
)

# Update the FIGO_stage names for better readability
data_long$FIGO_stage <- gsub("mean_FIGO_", "FIGO ", data_long$FIGO_stage)

# Create the plot
plot_FIGO_prevalence <- 
  ggplot(data_long, aes(x = age_interval, y = prevalence, color = FIGO_stage, group = FIGO_stage)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Mean FIGO Prevalence by Age Interval",
    x = "Age Interval",
    y = "Mean Prevalence (%)",
    color = "FIGO Stage"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################################################################################
## Plotting mean FIGOs:
mean_FIGO <- sim_result[["No Intervention"]]$mean_FIGO

# Reshape the data into a long format
data_long <- tidyr::pivot_longer(
  mean_FIGO,
  cols = starts_with("mean_FIGO"),
  names_to = "FIGO_stage",
  values_to = "diagnosed"
)

# Update the FIGO_stage names for better readability
data_long$FIGO_stage <- gsub("diag_FIGO_", "FIGO ", data_long$FIGO_stage)

# Create the plot
plot_mean_FIGO <- 
  ggplot(data_long, aes(x = age_interval, y = diagnosed, 
                        color = FIGO_stage, group = FIGO_stage)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Mean FIGO by Age Interval",
    x = "Age Interval",
    y = "Mean",
    color = "FIGO Stage"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
################################################################################



################################################################################
## Plotting mean Diagnosed:
mean_Diagnosed <- sim_result[["No Intervention"]]$mean_Diagnosed

# Reshape the data into a long format
data_long <- tidyr::pivot_longer(
  mean_Diagnosed,
  cols = starts_with("mean_Diagnosed_FIGO"),
  names_to = "FIGO_stage",
  values_to = "diagnosed"
)

# Update the FIGO_stage names for better readability
data_long$FIGO_stage <- gsub("diag_FIGO_", "FIGO ", data_long$FIGO_stage)

# Create the plot
plot_mean_Diagnosed_FIGO <- 
  ggplot(data_long, aes(x = age_interval, y = diagnosed, 
                        color = FIGO_stage, group = FIGO_stage)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Mean Diagnosed FIGO by Age Interval",
    x = "Age Interval",
    y = "Mean",
    color = "FIGO Stage"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
################################################################################



# Create plots for each measure
plot_CN1_incidences <- plot_comparison(combined_data, "CN1_incidences")
plot_CN2_incidences <- plot_comparison(combined_data, "CN2_incidences")
plot_CN3_incidences <- plot_comparison(combined_data, "CN3_incidences")
plot_CC_incidences <- plot_comparison(combined_data, "CC_incidences")
plot_HPV_prevalences <- plot_comparison(combined_data, "HPV_prevalences")
plot_CC_mortality <- plot_comparison(combined_data, "CC_mortality")
#plot_CC_by_diff_mortality <- plot_comparison(combined_data, "CC_by_diff_mortality")

# Display plots
print(plot_CN1_incidences)
print(plot_CN2_incidences)
print(plot_CN3_incidences)
print(plot_CC_incidences)
print(plot_HPV_prevalences)
print(plot_CC_mortality)
#print(plot_FIGO_prevalence)
#print(plot_mean_FIGO)
print(plot_mean_Diagnosed_FIGO)
#print(plot_CC_by_diff_mortality)


if (numb_of_sims >=60) {
  ################################################################################
  # For number of simulations of 60 we can analize the cost results to check
  # whether there is a numerical artifact or logic code problem producing
  # a tendency of decreas tc_hat_undisc along simulations:
  #average_cost <-
  #  other_mean_mortality_result[["No Intervention"]]$tc_hat_undisc$`sim[[i]][[name_level_of_sim]]`
  #average_cost <-
  #  sim_result[["No Intervention"]]$tc_hat_undisc$`sim[[i]][[name_level_of_sim]]`
  
  average_cost <-
    sim_result[["No Intervention"]]$tc_hat_undisc$tc_hat_undisc
  
  # Calculate confidence intervals for groups of 10 simulations
  grouped_means <- tapply(average_cost, (seq_along(average_cost) - 1) %/% 10, mean)
  grouped_sd <- tapply(average_cost, (seq_along(average_cost) - 1) %/% 10, sd)
  group_size <- 10
  z_value <- 1.96 # for 95% confidence
  
  # Calculate CI for each group
  CI <- grouped_means + c(-1, 1) * (z_value * (grouped_sd / sqrt(group_size)))
  
  # Plot the moving average
  library(zoo)
  moving_avg <- rollmean(average_cost, 10, align = "center")
  
  plot(average_cost, type = "l", main = "Average Cost over Simulations")
  lines(moving_avg, col = "red")
}

df <- sim_result[["No Intervention"]]$TR %>% select(FIGO.I, FIGO.II, FIGO.III, FIGO.IV) 
# select(FIGO.I, FIGO.II, FIGO.III, FIGO.IV) and summarize by columns
df <- df %>% summarise(across(everything(), sum, na.rm = TRUE))

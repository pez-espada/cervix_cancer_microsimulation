## ---- Preamble 
################################################################################
# This code is a modified version of the original code from:
# [https://github.com/DARTH-git/Microsimulation-tutorial] (Krijkamp et al 2018 
# Sick-Sicker model). 
# programmed by Carlos Dommar D'Lima - carlos.dommar@gmail.com
# This code extends the "sick-sicker" model of the original authors to a
# multi-state cervix cancer model
# CORRECTED TRANSITIONS ADDED ON 2025/04/14
# v.02 Screening by Cytology added (Work In Progress)
################################################################################
rm(list = ls())
library(tidyverse)
#library(future)

## For debuging purposes:
#options(error = recover)
## after debugging, you can set the error option back to the default:
#options(error = NULL)

## To prevent conflicts in the parallel environment:
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
# OLD TRANSITIONS:
my_Probs_old <- readRDS(file = "./data/probs.rds") # natural history transition matrix
my_Probs2_old <- readRDS(file = "./data/probs2.rds") # vaccination transition matrix
################################################################################

# CORRECTED TRANSITIONS (since 2025/04/14):
library(readxl)
my_Probs <- read_excel("data/corrected_transitions_20250414/Probs_20250414.xls")
my_Probs <- my_Probs %>% dplyr::rename(Age.group = `Age group`)
#my_Probs2 <- read_excel("data/corrected_transitions_20250414/Probs2_20250414.xlsx")
#my_Probs2<- my_Probs2 %>% dplyr::rename(Age.group = `Age group`)


################################################################################
adjust_infection_probs <- function(my_Probs, infection_reduction = 0.7) {
  # Convert to data frame
  my_Probs <- as.data.frame(my_Probs)
  my_Probs_adjusted <- my_Probs
  
  # Add 'state' column from column names (excluding the first column)
  my_Probs_adjusted$state <- names(my_Probs_adjusted)[2:ncol(my_Probs_adjusted)]
  
  # Adjust infection probability for "Well" state
  my_Probs_adjusted[my_Probs_adjusted$state == "Well", "HR.HPV.infection"] <- 
    my_Probs_adjusted[my_Probs_adjusted$state == "Well", "HR.HPV.infection"] * (1 - infection_reduction)
  
  # Recalculate probability of staying "Well"
  my_Probs_adjusted[my_Probs_adjusted$state == "Well", "Well"] <- 
    1 - (my_Probs_adjusted[my_Probs_adjusted$state == "Well", "HR.HPV.infection"] + 
           my_Probs_adjusted[my_Probs_adjusted$state == "Well", "Other.Death"])
  
  # Drop 'state' column
  my_Probs_adjusted$state <- NULL
  
  return(my_Probs_adjusted)
}
################################################################################

my_Probs2 <- adjust_infection_probs(my_Probs, infection_reduction = 0.7)


## Obtaining 'my_Probs2' from 'my_Probs' programatically (Sandra's code):
#infection_reduction <- 0.7 # due to vaccination
#my_Probs <- my_Probs %>% as.data.frame()
#my_Probs2 <- my_Probs
#my_Probs2$state <- names(my_Probs2[2:length(my_Probs2)])
#my_Probs2[my_Probs2$state == "Well", "HR.HPV.infection"  ] <- 
#  my_Probs2[my_Probs2$state=="Well", "HR.HPV.infection"  ]*(1 - infection_reduction)
#my_Probs2[my_Probs2$state == "Well", "Well" ] <- 
#  1-(my_Probs2[my_Probs2$state == "Well", "HR.HPV.infection"] + my_Probs2[my_Probs2$state == "Well", "Other.Death"])
#my_Probs2$state <- NULL
##Test 'adjust_infection_probs()' function:
#my_Probs_adjusted <- adjust_infection_probs(my_Probs, infection_reduction = 0.7)
#identical(my_Probs_adjusted, my_Probs2) # if TRUE then they're identical

# vaccination 2 associated immunity transition matrix
my_Probs2_nat_immunity <- readRDS(file = "./data/probs3.rds") 
my_Probs4 <- readRDS(file = "./data/probs3.rds") # vaccination transition matrix
my_Probs9 <- readRDS(file = "./data/probs3.rds") # vaccination transition matrix

## arbitrary correction of an obvious error (for old matrices);
#my_Probs2$Other.Death[my_Probs2$Other.Death == 8.150000e+08] <- 8.150000e-08
#my_Probs2 <- my_Probs2 %>%
#  dplyr::rename("CC_Death" = "Death")

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
my_Probs <- my_Probs_cleaning_Func(Probs_matrix = my_Probs)
my_Probs <- my_Probs %>% as.data.frame() #convert back to data.frame (no needed?)

my_Probs2 <- my_Probs_cleaning_Func(Probs_matrix = my_Probs2)
my_Probs2 <- my_Probs2 %>% as.data.frame() #convert back to data.frame (no needed?)

my_Probs2_nat_immunity <- my_Probs_cleaning_Func(Probs_matrix = my_Probs2_nat_immunity)
my_Probs2_nat_immunity <- my_Probs2_nat_immunity %>% as.data.frame() #convert back to data.frame (no needed?)

my_Probs4 <- my_Probs_cleaning_Func(Probs_matrix = my_Probs4)
my_Probs4 <- my_Probs4 %>% as.data.frame() #convert back to data.frame (no needed?)

my_Probs9 <- my_Probs_cleaning_Func(Probs_matrix = my_Probs9)
my_Probs9 <- my_Probs9 %>% as.data.frame() #convert back to data.frame (no needed?)
################################################################################
## ----Model Parameters
n_i <- 10^4               # number of simulated individuals
#n_t <- 3                  # time horizon, 3 cycles (it starts from 1)
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
v_n <- colnames(my_Probs)
v_n <- v_n[-c(1,14,15)]
n_s   <- length(v_n)                # the number of health states
v_M_1 <- rep("H", n_i)              # everyone begins in the healthy state 
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

n_dose_vacc2 <- 2
cost_vacc2 <- 34.8 * n_dose_vacc2 # cost per dose
#cost_vacc2 <- 69.6
cost_vacc4 <- 0
cost_vacc9 <- 0
################################################################################

################################################################################
## ---- FUNCTIONS -----                                                       ##  
### For extracting the probabilities of transitions given the transition matrix:
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

################################################################################
### ---- Probability Function ---- (original)                                 ##
### The Probs function that updates the transition probabilities of every cycle:
#Probs <- function(M_it, my_Probs) {
#  n_s <- length(v_n)
#  n_i <- length(M_it)
#  m_P_it <- matrix(NA, n_s, n_i) 
#  rownames(m_P_it) <- v_n
#  for (i in 1:length(v_n)) {
#    state_mask <- !is.na(M_it) & M_it == v_n[i]
#    
#    if (sum(state_mask) > 0) {
#      m_P_it[, state_mask] <- 
#        lapply(X = v_n, function(x) trans_prb(P = my_Probs, state1 =
#                                                v_n[i], state2 = x)) %>%
#        unlist()
#    } else {
#      ## Debugging:
#      #cat("State", v_n[i], "is not present in M_it at this time step\n")
#    }
#  }
#  if (any(is.na(m_P_it))) {
#    # Diagnostic message
#    #cat("Transition probabilities contain NA values\n")
#  }
#  
#  ifelse(colSums(m_P_it, na.rm = TRUE) >= .991, 
#         return(t(m_P_it)), 
#         stop("Probabilities do not sum to 1"))
#}
################################################################################


# New Probs fnct:
################################################################################
Probs <- function(M_it, my_Probs) {
  n_s <- length(v_n)
  n_i <- dim(M_it)[1]
  m_P_it <- matrix(NA, n_s, n_i) 
  ID <- M_it$ID 
  #ID <- as.integer(ID) # fix for parallelization
  M_it<-M_it$health_state 
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
  
  #if(colSums(m_P_it, na.rm = TRUE) >= .991){
  #if(any(colSums(m_P_it, na.rm = TRUE) > 1)) {
  #  stop("Probabilities do not sum to 1")  # HERE IT BREAKS AT CYCLE 5!!
  #}else{
    t_m_P_it<-t(m_P_it)
    t_m_P_it<-cbind(ID,t_m_P_it)
    return(t_m_P_it)
  #}
}
################################################################################

################################################################################
Probs_3_optimized <- function(M_it, v_n, n_i, seed, prob_matrix, prob_matrix_2, 
                              prob_matrix_2_nat_immunity, prob_matrix_4, 
                              prob_matrix_9, vacc_lbl, age) {
  M_it <- tibble(ID = as.integer(1:length(M_it)), health_state = M_it)
  
  # merge data with vacc_lbl:
  M_it <- M_it %>% dplyr::left_join(vacc_lbl, by = "ID")
  
  # Create a matrix to store the probabilities
  P_combined <- matrix(NA, nrow = n_i, ncol = length(v_n) + 1) # additional column for ID
  P_combined[, 1] <- 1:n_i # ID column
  
  current_row <- 1
  # Use data.table joins and vectorized selection
  for (vacc_status in unique(M_it$vacc_state)) {
    # subset all individuals with a specific vacc_state
    #state_subset <- M_it[vacc_state == vacc_status]
    state_subset <- M_it %>% dplyr::filter(vacc_state == vacc_status)
    
    # determine the prob matrix corresponding to the vacc_state of the subset
    prob_mat <- switch(vacc_status,
                       "no_vacc" = prob_matrix,
                       "vacc_2" = if (any(state_subset$immuned)) {
                         prob_matrix_2_nat_immunity
                       } else {
                         prob_matrix_2
                       },
                       "vacc_4" = prob_matrix_4,
                       "vacc_9" = prob_matrix_9,
                       stop("Unknown vaccination status detected")
    )
    # Drop unnecessary columns:
    #prob_mat <- prob_mat[, c("Age.group", "Lower", "Larger") := NULL]
    prob_mat <-  prob_mat %>% 
      dplyr::select(-c(Age.group, Lower, Larger))
    
    #############################################################################
    ## let's try something different: use the tested Probs() function
    ## for each state_subset
    ## Before calling the Probs() func lets prepare its arguments:
    M_it_2 <- state_subset  %>% dplyr::select(ID, health_state)
    
    P_combined[current_row : (dim(state_subset)[1] + current_row - 1), ] <-  
      Probs(M_it = M_it_2, my_Probs = prob_mat)
    #Test_Prob <-   Probs(M_it = M_it_2, my_Probs = prob_mat)
    
    current_row <- current_row + dim(state_subset)[1]
    #############################################################################
  } # for vacc_status
  
  # Sort by ID:
  P_combined <- P_combined[order(P_combined[,1]), ]
    
  P_combined <- P_combined[, -1]  # remove first ID column
  colnames(P_combined) <- v_n
  return(P_combined)
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
  
  ## Handle NA in probs
  #if (any(is.na(probs))) {
  #  warning("NA detected in transition probabilities, replacing with uniform distribution")
  #  probs[is.na(probs)] <- 1 / k
  #}
  
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
  #U[k, ] <- 1  # Force last row to be exactly 1
  if (any((U[k, ] - 1) > 1e-04))
    stop("error in multinom: probabilities do not sum to 1")
  ##############################################################################
  ##############################################################################
  
  ### Random sampling, binning, and moving states: 
  for (j in 1:m) {
    # RANDOM GEN CODE HERE:
    un <- rep(runif(n), rep(k, n)) # repeat `runif(n)` `rep(k,n)`times
    # this create a numeric of `n_i x n_s` that 
    # sample  an uniformed distributed number 
    # between 0 and 1. The generated random number
    # repeats itself `n_s` times and then another 
    # rand unif number is drawn. This process is 
    # carried out `n_i` times. NOTE: every time
    # runif() is run it produces a new random sample
    # i.e. it does not seem dependent on the seed
    
    ## Here's where we choose the individuals' next states:
    #ran[, j] <- lev[1 + colSums(un > U)]
    ran[, j] <- lev[1 + pmin(colSums(un > U), k - 1)]
  }
  #cat("\n")
  #cat("Unique states computed by samplev() is/are:\n")
  #ran %>% unique() %>% print()
  #cat("\n")
  
  return(ran)
}
################################################################################

################################################################################
## ---- Costs Function ----                                                   ##
### Costs Function
# The `Costs_per_Cancer_Diag` function estimates the costs of a diagnose 
# individual due to cancer symptoms (FIGO.I-IV) at every cycle. 
# This cost is only charged once in the patient's lifetime.
# NOTE: need to decide if the cost is applied on current time `t` or `t+1` as it is now.
Costs_per_Cancer_Diag <- function (M_it, cost_Vec, symptomatics, 
                                   time_iteration, Trt = FALSE) {
  #ci_t <- 0
  c_it <- rep(0, length(M_it))
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
Effs <- function (M_it, trt = FALSE, cl = 1, utilityCoefs) {
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
# If the TryCatch gives problems, just overrate it:
#for (i in 1:length(utilityCoefs)) {
#  u_it[M_it == v_n[i]] <- utilityCoefs[i]   # update the utility if healthy
#}
return(u_it)
}
################################################################################


ensure_library("dplyr", "tidyverse", "purrr", "data.table")
################################################################################
# --- Function receives a column with current state of `n_i` individuals and gives
# a dataframe with `ID, TimeStep`, `state`, and `RecoveredFromState` columns.
# The function also updates the global vector `global_diagnosed` with the IDs of
# individuals who have been diagnosed.
# RANDOM FUNCTION:
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
        # (RANDOM GEN CODE HERE):
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
    #new_entries$ID <- as.integer(new_entries$ID) # fix for parallelization
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
# Altenative function 4:
new_cases_2 <- function(state1, state2, Tot_Trans_per_t) {
  Tot_Trans_per_t_tbl <- as_tibble(Tot_Trans_per_t)
  
  if (length(state1) == 1) {
    transition_column <- paste0(state1, "->", state2)
    
    if (transition_column %in% colnames(Tot_Trans_per_t_tbl)) {
      transition_cases <- Tot_Trans_per_t_tbl %>%
        select(all_of(transition_column)) %>%  
        mutate(age = row_number() + 10, cycle = age - 9) %>% 
        add_row(!!rlang::sym(transition_column) := 0, age = 10, cycle = 1, .before = 1) %>% 
        slice(-n()) %>%
        mutate(age = 10:(10 + n() - 1), cycle = age - 9)
    } else {
      warning(paste0("Transition '", transition_column, "' not found! Using a column of zeros."))
      transition_cases <- tibble(
        !!rlang::sym(transition_column) := rep(0, nrow(Tot_Trans_per_t_tbl)),
        #age = row_number() + 10, cycle = age - 9
        age = seq_len(nrow(Tot_Trans_per_t_tbl)) + 10,  
        cycle = seq_len(nrow(Tot_Trans_per_t_tbl)) + 1,
        cat("cycle, ", cycle, "\n")
      ) %>%
        add_row(!!rlang::sym(transition_column) := 0, age = 10, cycle = 1, .before = 1) %>% 
        slice(-n()) %>%
        mutate(age = 10:(10 + n() - 1), cycle = age - 9)
    }
    return(transition_cases)
  } else if (length(state1) > 1) {
    transition_columns <- paste0(state1, "->", state2)
    
    existing_cols <- intersect(transition_columns, colnames(Tot_Trans_per_t_tbl))
    missing_cols <- setdiff(transition_columns, colnames(Tot_Trans_per_t_tbl))
    
    if (length(missing_cols) > 0) {
      warning(paste0("Some transitions not found: ", paste(missing_cols, collapse = ", "), ". Using columns of zeros for these."))
    }
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
      select(all_of(existing_cols)) %>%
      bind_cols(missing_df) %>%
      rowwise() %>%
      mutate(!!rlang::sym(paste0(state2, "_per_t")) := sum(c_across(everything()), na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(age = row_number() + 10, cycle = age - 9) %>%
      add_row(!!rlang::sym(paste0(state2, "_per_t")) := 0, age = 10, cycle = 1, .before = 1) %>%
      slice(-n()) %>%
      mutate(age = 10:(10 + n() - 1), cycle = age - 9)
    
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
extract_screening_days <- function(strategy_string) {
  # Extract the age range and period using regular expressions
  matches <-
    regmatches(strategy_string, regexec("(\\d+)-(\\d+).*?(\\d+)", 
                                        strategy_string))[[1]]
  start_age <- as.integer(matches[2])
  end_age <- as.integer(matches[3])
  period <- as.integer(matches[4])
  
  # Generate screening days
  days <- seq(start_age, end_age, by = period)
  return(days)
}


## additonal initialization:
CIN1_followup_IDs <- character(0)
log_transitions <- FALSE  # toggle to TRUE for debugging
state_transition_log <- data.table()  # only used if log_transitions == TRUE

## Example usage
#extract_screening_days("25-29 cito 3 anys")
## Output: [1] 25 28
################################################################################





 

################################################################################
################################################################################
################################################################################
## THE MICROSIMULATION MAIN FUNCTION
# This version stacks solution of simulations but produces a list with stacked elements
MicroSim <- function(strat=strat, 
                     numb_of_sims = 20,
                     v_M_1,
                     n_i, 
                     n_t, 
                     v_n, 
                     d_c, 
                     d_e, 
                     TR_out = TRUE, 
                     TS_out = TRUE, 
                     Trt = FALSE,  
                     Pmatrix,
                     use_parallel = FALSE, 
                     reproducible = TRUE, 
                     master_seed = 123,
                     cost_vacc2, 
                     cost_vacc4, 
                     cost_vacc9,
                     #screening_strategies, 
                     screening_coverage ,
                     vacc_coverage,
                     ScreenPrice,
                     costCoeff_md,
                     citoSpecif
)
{
  cat("I HAVE ENTERED THE SIMULATOR \n")
  cl <- NULL  # Ensure cl exists in all cases
  
  if (use_parallel) {
    n_cores <- min(detectCores() - 1, 20)  # Try using 8 or fewer cores
    cat("Number of cores: ", n_cores, "\n")
    # 6-hours timeout to prevent socket drop issues
    cl <- makeCluster(n_cores, timeout = 6*60*60) 
    
    clusterExport(cl, c("Costs_per_Cancer_Diag", "Effs", "trans_prb", "Probs",
                        "vacc_lbl", "Probs_3_optimized",
                        "my_Probs", "my_Probs2","my_Probs4", "my_Probs9", 
                        "my_Probs2_nat_immunity", "utilityCoefs", "v_n", "samplev", 
                        "my_age_prob_matrix_func","diagnose_column", 
                        "update_column", "states_to_check", 
                        "symptom_prob_vec", "survival_prob_vec", #"global_diagnosed", 
                        "cost_Vec", "new_cases_2",
                        "extract_screening_days", "IDs", "screenSensi",
                        "screenProbs"))
    
    registerDoParallel(cl)
  } else {
    registerDoSEQ()  # Runs sequentially for debugging
  }
  
  # Generate independent seeds for each simulation run
  if (reproducible && !is.null(master_seed)) {
    set.seed(master_seed)
    seeds <- sample.int(1e6, numb_of_sims)  # Generate unique seeds
    cat("THE RANDOM SEEDS ARE:", seeds, "\n")
  } else {
    seeds <- NULL  # No reproducibility
  }
  
  # Some Initiliazations: 
  simulation_results <- list() 
  
  # calculate the cost discount weight based on the discount rate d_c 
  v_dwc <- 1 / (1 + d_c) ^ (0:(n_t-1))   
  # calculate the QALY discount weight based on the discount rate d_e                                             
  v_dwe <- 1 / (1 + d_e) ^ (0:(n_t-1))   
  
  # If vaccination, apply vaccination cost to those vaccinated individuals
  # ONLY ONCE per sim batch:
  vacc_cost <- rep(0, n_i)
  if (any(vacc_coverage != 0)) { 
    # vacc_covverage pos1 is bivalent, pos2 is 4-valent and pos3 is 9-valent
    if (vacc_coverage[1] != 0) {
      vaccinated_id <- which(vacc_lbl$vacc_state == "vacc_2")
      vacc_cost[vaccinated_id] <- cost_vacc2
      #cat("we have vaccinated here!\n")
    }
    if (vacc_coverage[2] != 0) {
      vaccinated_id <- which(vacc_lbl$vacc_state == "vacc_4")
      vacc_cost[vaccinated_id] <- cost_vacc4
    }
    if (vacc_coverage[2] != 0) {
      vaccinated_id <- which(vacc_lbl$vacc_state == "vacc_9")
      vacc_cost[vaccinated_id] <- cost_vacc9
    }
  }
  
  # source("./R/sumarize_results_by_Strategy_Func.R")
  #source("./R/sumarize_results_by_Strategy_Func_v2.R")
  
  stacked_results <- NULL
  # initialize joined_batches_per_strategy
  #joined_batches_per_strategy <- list()
  
  # Parallel processing using foreach - batch-level loop
  simulation_results <- 
    foreach(sim = 1:numb_of_sims, .packages = c("dplyr", 
                                                "tidyr", "purrr", 
                                                "data.table") ) %dopar%
    { 
      
      current_sim <- sim  # Create an isolated snapshot
      cat(">>> STARTING SIM:", sim, "\n")
      
      
      seed <- seeds[sim]
      
      cat("\n")
      cat("\n")
      cat("\n")
      cat("----------------------------------------------------\n")
      cat("----------------------------------------------------\n")
      cat("Running simulation", sim, "with seed", seeds[sim], "\n")
      cat("----------------------------------------------------\n")
      
      # Initialize a global vector to store all diagnosed individuals
      global_diagnosed <<- integer()
      symptomatics <-
        data.frame(ID = integer(), TimeStep = integer(), 
                   DiagnosedState = character(), 
                   RecoveredFromState = logical(), stringsAsFactors = FALSE)
      
      
      # Create the matrix capturing the state name/costs/health outcomes 
      # for all individuals at each time point:
      #m_M <- m_C <- m_E <-  matrix(nrow = n_i, ncol = (n_t + 1), 
      m_M <- m_C <- m_E <- 
        matrix(nrow = n_i, ncol = (n_t), 
               dimnames = list( 1:n_i, 
                                paste0("cycle_", 1:(n_t),
                                       sep = "")))  
      
      m_M[, 1] <- v_M_1  # Indicate the initial health state   
      
      
      # estimate costs per individual for the initial health state
      m_C[, 1] <- Costs_per_Cancer_Diag(M_it = m_M[, 1], 
                                        symptomatics = symptomatics,
                                        time_iteration = 1,
                                        cost_Vec = cost_Vec,  
                                        Trt)             
      # account for vaccination cost:
      m_C[, 1] <- m_C[, 1] + vacc_cost
      
      # estimate QALYs per individual for the initial health state 
      m_E[, 1] <- Effs(m_M[, 1], Trt, utilityCoefs = utilityCoefs)  
      
      stored_list <- list()
      
      ### TEST 20250509
      # Before the t loop
      current_age_group <- NA
      my_age_prob_matrix <- NULL
      my_age_prob_matrix_2 <- NULL
      my_age_prob_matrix_2_nat_immunity <- NULL
      my_age_prob_matrix_4 <- NULL
      my_age_prob_matrix_9 <- NULL
      
      
      # --- Initialize screening logic state per simulation ---
      cyto_screening_days <- extract_screening_days(strat)
      
      detected_IDs <- character()              # CIN2+ detections (screened or from CIN1 follow-up)
      CIN1_followup_IDs <- character()         # Under CIN1 follow-up
      screened_registry <- data.table(sim = integer(), age = integer(), ID = character())  # Per (sim, age) uniqueness
      cost_log <- data.table()
      
      
      
      ########################################################################
      #################### run over all the cycles ########################### 
      # Loop runs over all the cycles of the simulation. It updates the
      # health state of each individual at each cycle, estimates the costs and
      # QALYs per individual at each cycle, and stores the transitions across
      # states for each individual at each cycle.
      for(t in 1:(n_t-1)) {
        ######################################################################
        # Select the transition matrix based on the cycle `n_t`:
        # Since our age intervals start at 10 years old,
        age_in_loop <- t + 9
        
        #cat("Simulation:", sim, "Cycle:", t, ", ", "Age:", age_in_loop, ", ", "seed:", seed, "\n")
        
        ########################################################################
        
        ######################################################################## 
        # Computation of Symptomatics:
        # RANDOM FUNCTION:
        new_entries <- diagnose_column(m_M[, t], t) 
        
        if (!is.null(new_entries)) {
          stored_list[[t]] <- new_entries
        }
        if (nrow(new_entries) > 0) {
          symptomatics <- bind_rows(symptomatics, new_entries)
        }
        ######################################################################## 
        
        ### TEST 20250509:
        # Inside the loop
        # floor()  es la parte entera de la división
        # Note: we use 'age in loop + 1' because we ask for the transitions to 
        # move states ahead in the future t + 1. This is because of model design.
        age_group <- floor((age_in_loop + 1) / 5)
        
        if (is.na(current_age_group) || age_group != current_age_group) {
          current_age_group <- age_group
          
          my_age_prob_matrix <- my_age_prob_matrix_func(my_Probs, age_in_loop + 1)
          rownames(my_age_prob_matrix) <- v_n
          
          my_age_prob_matrix_2 <- my_age_prob_matrix_func(my_Probs2, age_in_loop + 1) %>%
            dplyr::mutate(
              Age.group = ifelse(Age.group == "11-14", "10-14", Age.group),
              Lower = ifelse(Lower == "11", "10", Lower)
            )
          rownames(my_age_prob_matrix_2) <- v_n
          
          my_age_prob_matrix_2_nat_immunity <- my_age_prob_matrix_func(my_Probs2_nat_immunity, age_in_loop + 1) %>%
            dplyr::mutate(
              Age.group = ifelse(Age.group == "11-14", "10-14", Age.group),
              Lower = ifelse(Lower == "11", "10", Lower)
            )
          rownames(my_age_prob_matrix_2_nat_immunity) <- v_n
          
          my_age_prob_matrix_4 <- my_age_prob_matrix_func(my_Probs4, age_in_loop + 1) %>%
            dplyr::mutate(
              Age.group = ifelse(Age.group == "11-14", "10-14", Age.group),
              Lower = ifelse(Lower == "11", "10", Lower)
            )
          rownames(my_age_prob_matrix_4) <- v_n
          
          my_age_prob_matrix_9 <- my_age_prob_matrix_func(my_Probs9, age_in_loop + 1) %>%
            dplyr::mutate(
              Age.group = ifelse(Age.group == "11-14", "10-14", Age.group),
              Lower = ifelse(Lower == "11", "10", Lower)
            )
          rownames(my_age_prob_matrix_9) <- v_n
        } #endif
        
        
        # Extract the transition probabilities of each individuals at cycle t
        # given the individual current state and the corresponding 
        # transition probability matrix that depends on age:
        # Next time (t+1) transition:
        ## m_P is a (n_i x n_s) matrix with the probabilities of transitioning
        #m_P <- Probs(M_it =  m_M[, t], my_Probs = my_age_prob_matrix)
        
        # Function to obtain individual transition probabilities based on 
        # current state and the prob of transition one cycle/t ahead.
        m_P <- Probs_3_optimized(M_it = m_M[, t], v_n = v_n, n_i = n_i, 
                                 prob_matrix = my_age_prob_matrix, 
                                 prob_matrix_2 = my_age_prob_matrix_2,
                                 prob_matrix_2_nat_immunity = my_age_prob_matrix_2_nat_immunity,
                                 prob_matrix_4 = my_age_prob_matrix_4, 
                                 prob_matrix_9 = my_age_prob_matrix_9, 
                                 vacc_lbl = vacc_lbl,
                                 age = age_in_loop)
        #cat("Dimension of m_P is (outside the function): ",dim(m_P),"\n")
        
        # Make actual transition by  random sampling (RANDOM FUNCTION): 
        m_M[, t + 1] <- samplev(probs = m_P, m = 1)  # sample the next health state 
        # and store that state in  
        ## matrix m_M 
        #cat("Dimension of m_M is ",dim(m_M),"\n")
        ########################################################################    
        
        # m_M[, t + 1] <- update_column(m_M[, t], new_entries)
        next_col <- m_M[, t + 1]
        next_col <- update_column(m_M[, t], new_entries, next_col)
        
        # Ensure next_col updates for `Survival` are preserved after sampling
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
        
        m_E[, t + 1] <- # estimate QALYs per individual during cycle t + 1
          Effs( m_M[, t + 1], Trt, 
                utilityCoefs = utilityCoefs)                   
        ########################################################################    
        ########################################################################    
        #cat('\r', paste(round(t/n_t * 100),          # display the 
        #                "% done\n", sep = " "))      # progress of  the simulation
        
        ### # ########################################################################    
        ### # ## ------------------- Cytology Screening Block ----------------------
        ### # ## Version 2-F (Fixed Follow-up Timing)
        ### # ## Version 2-G (Now they have similar costs, check with commented code
        ### # at the end of CIN1 follow-up )
        ### # ## Version 2-I 
        ### Add newly diagnosed CIN1 cases from previous cycle to follow-up list
        ##if (exists("newly_diagnosed_CIN1_IDs") && length(newly_diagnosed_CIN1_IDs) > 0) {
        ##  CIN1_followup_IDs <- unique(c(CIN1_followup_IDs, newly_diagnosed_CIN1_IDs))
        ##  newly_diagnosed_CIN1_IDs <- NULL  # reset for this cycle
        ##} else {
        ##  newly_diagnosed_CIN1_IDs <- NULL
        ##}
        ##
        ### ------------------- Cytology Screening Block ----------------------
        ##if (age_in_loop %in% cyto_screening_days) {
        ##  cat(" Performing cytology screening at age", age_in_loop, "for sim", current_sim, "\n")
        ##  
        ##  not_detected <- !IDs %in% detected_IDs
        ##  already_screened <- screened_registry[sim == current_sim & age == age_in_loop, ID]
        ##  eligible_ids <- setdiff(IDs[not_detected], already_screened)
        ##  
        ##  if (length(eligible_ids) > 0) {
        ##    eligible_screened <- runif(length(eligible_ids)) < screening_coverage
        ##    screened_ids <- eligible_ids[eligible_screened]
        ##    
        ##    screened_registry <- rbind(screened_registry, data.table(
        ##      sim = current_sim,
        ##      age = age_in_loop,
        ##      ID = screened_ids
        ##    ))
        ##    
        ##    cost_log <- rbindlist(list(cost_log, data.table(
        ##      sim = current_sim,
        ##      age = age_in_loop,
        ##      ID = screened_ids,
        ##      cost_type = strat,
        ##      cost = ScreenPrice
        ##    )), use.names = TRUE)
        ##    
        ##    screened_states <- m_M[match(screened_ids, IDs), t]
        ##    state_indices <- match(screened_states, v_n)
        ##    diagnose_probs <- screenSensi[state_indices]
        ##    diagnosed <- runif(length(diagnose_probs)) < diagnose_probs
        ##    
        ##    diagnosed_ids <- screened_ids[diagnosed]
        ##    diagnosed_states <- screened_states[diagnosed]
        ##    diagnosed_indices <- state_indices[diagnosed]
        ##    
        ##    if (length(diagnosed_ids) > 0) {
        ##      # Log diagnosis cost at current cycle
        ##      followup_costs <- costCoeff_md[diagnosed_indices]
        ##      cost_log <- rbindlist(list(cost_log, data.table(
        ##        sim = current_sim,
        ##        age = age_in_loop,
        ##        ID = diagnosed_ids,
        ##        cost_type = paste0("diagnosed_by_cyto_", diagnosed_states),
        ##        cost = followup_costs
        ##      )), use.names = TRUE)
        ##      
        ##      # Track CIN1 diagnosed IDs but DO NOT add to follow-up yet — defer to next cycle
        ##      CIN1_diagnosed <- diagnosed_states == "CIN1"
        ##      newly_diagnosed_CIN1_IDs <- diagnosed_ids[CIN1_diagnosed]
        ##      
        ##      # CIN2+ detected — update detected list immediately
        ##      CIN2plus_mask <- diagnosed_states %in% c("CIN2", "CIN3", "FIGO.I", "FIGO.II", "FIGO.III", "FIGO.IV")
        ##      CIN2plus_new <- diagnosed_ids[CIN2plus_mask & !(diagnosed_ids %in% detected_IDs)]
        ##      detected_IDs <- unique(c(detected_IDs, CIN2plus_new))
        ##      
        ##      # Recovery logic unchanged
        ##      recovery_probs <- screenProbs[diagnosed_indices]
        ##      recovery_mask <- runif(length(diagnosed_ids)) < recovery_probs
        ##      
        ##      if (any(recovery_mask)) {
        ##        recovered_ids <- diagnosed_ids[recovery_mask]
        ##        recovered_states <- diagnosed_states[recovery_mask]
        ##        recovered_rows <- match(recovered_ids, IDs)
        ##        
        ##        to_survival <- recovered_states %in% c("FIGO.I", "FIGO.II", "FIGO.III", "FIGO.IV")
        ##        to_H        <- recovered_states %in% c("CIN1", "CIN2", "CIN3")
        ##        
        ##        #m_M[recovered_rows[to_survival], t]     <- "Survival"
        ##        m_M[recovered_rows[to_survival], t + 1] <- "Survival"
        ##        #m_M[recovered_rows[to_H],        t]     <- "H"
        ##        m_M[recovered_rows[to_H],        t + 1] <- "H"
        ##        
        ##        cost_log <- rbindlist(list(cost_log, data.table(
        ##          sim = current_sim,
        ##          #age = age_in_loop,
        ##          age = age_in_loop + 1,  # Recovery occurs after state transition
        ##          ID = recovered_ids,
        ##          cost_type = paste0("recovery_from_", recovered_states),
        ##          cost = 0
        ##        )), use.names = TRUE)
        ##      }
        ##    }
        ##  }
        ##}
        ##
        ### ----------------- CIN1 Follow-Up Block 1.3 ------------------
        ##if (length(CIN1_followup_IDs) > 0) {
        ##  followup_rows <- match(CIN1_followup_IDs, IDs)
        ##  
        ##  # Step 1: Log CIN1 follow-up if current state is still CIN1 and not yet detected
        ##  current_states <- m_M[followup_rows, t]
        ##  still_CIN1 <- current_states == "CIN1"
        ##  still_CIN1_IDs <- CIN1_followup_IDs[still_CIN1]
        ##  
        ##  # Only log follow-up cost if not already detected
        ##  still_CIN1_IDs <- setdiff(still_CIN1_IDs, detected_IDs)
        ##  
        ##  if (length(still_CIN1_IDs) > 0) {
        ##    cost_log <- rbindlist(list(cost_log, data.table(
        ##      sim = current_sim,
        ##      age = age_in_loop,
        ##      ID = still_CIN1_IDs,
        ##      cost_type = "CIN1_followup",
        ##      cost = costCoeff_md[match("CIN1", v_n)]
        ##    )), use.names = TRUE)
        ##  }
        ##  
        ##  # Step 2: Check progression to CIN2+ at next state
        ##  next_states <- m_M[followup_rows, t + 1]
        ##  progressed <- next_states %in% c("CIN2", "CIN3", "FIGO.I", "FIGO.II", "FIGO.III", "FIGO.IV")
        ##  
        ##  if (any(progressed)) {
        ##    progressed_IDs <- CIN1_followup_IDs[progressed]
        ##    # Only new progressions (not already treated)
        ##    progressed_IDs <- setdiff(progressed_IDs, detected_IDs)
        ##    
        ##    if (length(progressed_IDs) > 0) {
        ##      progressed_states <- next_states[match(progressed_IDs, CIN1_followup_IDs)]
        ##      progressed_indices <- match(progressed_states, v_n)
        ##      
        ##      cost_log <- rbindlist(list(cost_log, data.table(
        ##        sim = current_sim,
        ##        age = age_in_loop + 1,
        ##        ID = progressed_IDs,
        ##        cost_type = paste0("progressed_from_CIN1_", progressed_states),
        ##        cost = costCoeff_md[progressed_indices]
        ##      )), use.names = TRUE)
        ##      
        ##      detected_IDs <- unique(c(detected_IDs, progressed_IDs))
        ##      CIN1_followup_IDs <- setdiff(CIN1_followup_IDs, progressed_IDs)
        ##    }
        ##  }
        ##  
        ##  # Step 3: Regressed → stop follow-up
        ##  regressed <- next_states %in% c("H", "HPV.infection")
        ##  if (any(regressed)) {
        ##    regressed_IDs <- CIN1_followup_IDs[regressed]
        ##    CIN1_followup_IDs <- setdiff(CIN1_followup_IDs, regressed_IDs)
        ##  }
        ##}
        
       
        # --------------------------------------------------------------------
        # ------------------- Cytology Screening Block -----------------------
        # Version 2-I (Fixed Follow-up Timing + Recovery Handling)
        # --------------------------------------------------------------------
        
        # Defer CIN1 diagnosed IDs from previous cycle
        # Si en el ciclo anterior se detectaron nuevos CIN1, se incorporan al seguimiento
        if (exists("newly_diagnosed_CIN1_IDs") && length(newly_diagnosed_CIN1_IDs) > 0) {
          CIN1_followup_IDs <- unique(c(CIN1_followup_IDs, newly_diagnosed_CIN1_IDs))
          newly_diagnosed_CIN1_IDs <- NULL
        } else {
          newly_diagnosed_CIN1_IDs <- NULL
        }
        
        # Si la edad actual corresponde a un día de cribado:
        if (age_in_loop %in% cyto_screening_days) {
          cat("🧪 Performing cytology screening at age", age_in_loop, "for sim", current_sim, "\n")
          
          # Excluimos IDs ya detectados o ya cribados en esta edad y simulación:
          not_detected <- !IDs %in% detected_IDs
          already_screened <- screened_registry[sim == current_sim & age == age_in_loop, ID]
          eligible_ids <- setdiff(IDs[not_detected], already_screened)
          
          if (length(eligible_ids) > 0) {
            eligible_screened <- runif(length(eligible_ids)) < screening_coverage
            screened_ids <- eligible_ids[eligible_screened]
            
            # Log screening
            # Registramos que han sido cribados en esta edad/simulación:
            screened_registry <- rbind(screened_registry, data.table(
              sim = current_sim,
              age = age_in_loop,
              ID = screened_ids
            ))
            
            # Registramos el coste del cribado:
            cost_log <- rbindlist(list(cost_log, data.table(
              sim = current_sim,
              age = age_in_loop,
              ID = screened_ids,
              cost_type = strat,
              cost = ScreenPrice
            )), use.names = TRUE)
            
            # Screening results
            # Obtenemos el estado de salud actual y 
            # lo diagnosticamos con cierta sensibilidad/probabilidad:
            screened_states <- m_M[match(screened_ids, IDs), t]
            state_indices <- match(screened_states, v_n)
            diagnose_probs <- screenSensi[state_indices]
            diagnosed <- runif(length(diagnose_probs)) < diagnose_probs
            
            diagnosed_ids <- screened_ids[diagnosed]
            diagnosed_states <- screened_states[diagnosed]
            diagnosed_indices <- state_indices[diagnosed]
            
            if (length(diagnosed_ids) > 0) {
              # log diagnosis cost at current cycle
              # Registramos el coste del diagnóstico por citología
              followup_costs <- costCoeff_md[diagnosed_indices]
              cost_log <- rbindlist(list(cost_log, data.table(
                sim = current_sim,
                age = age_in_loop,
                ID = diagnosed_ids,
                cost_type = paste0("diagnosed_by_cyto_", diagnosed_states),
                cost = followup_costs
              )), use.names = TRUE)
              
              # Defer CIN1 diagnosed to next cycle
              # Track CIN1 diagnosed IDs but DO NOT add to follow-up yet — defer to next cycle
              # Guardamos los CIN1 diagnosticados para iniciar seguimiento en el próximo ciclo
              # (esto se hace para evitar que se inicien seguimientos en el mismo ciclo)
              CIN1_diagnosed <- diagnosed_states == "CIN1"
              newly_diagnosed_CIN1_IDs <- diagnosed_ids[CIN1_diagnosed]
              
              # CIN2+ detected — track immediately
              # Para individuos con lesiones CIN2+ o cáncer FIGO, registramos detección inmediata
              CIN2plus_mask <- diagnosed_states %in% c("CIN2", "CIN3", "FIGO.I", "FIGO.II", "FIGO.III", "FIGO.IV")
              CIN2plus_new <- diagnosed_ids[CIN2plus_mask & !(diagnosed_ids %in% detected_IDs)]
              detected_IDs <- unique(c(detected_IDs, CIN2plus_new))
              
              # ---------------- Recovery Block ----------------
              recovery_probs <- screenProbs[diagnosed_indices]
              recovery_mask <- runif(length(diagnosed_ids)) < recovery_probs
              
              if (any(recovery_mask)) {
                recovered_ids <- diagnosed_ids[recovery_mask]
                recovered_states <- diagnosed_states[recovery_mask]
                recovered_rows <- match(recovered_ids, IDs)
                
                # Apply recovery in next state
                # Los cánceres FIGO recuperan a estado de "Survival"
                to_survival <- recovered_states %in% c("FIGO.I", "FIGO.II", "FIGO.III", "FIGO.IV")
                # CIN1,2,3 regresan a estado "H"
                to_H        <- recovered_states %in% c("CIN1", "CIN2", "CIN3")
                
                # Para los individuos que se han recuperado tras diagnóstico y cuya patología era FIGO.X,
                # asignamos el estado "Survival" en el siguiente ciclo (t + 1)
                # NOTA: No recuperan inmediatamente; la recuperación se aplica al ciclo siguiente.
                if (any(to_survival)) {
                  m_M[recovered_rows[to_survival], t + 1] <- "Survival"
                }
                
                # Para los individuos recuperados que tenían CIN1, CIN2 o CIN3, se les asigna el estado "H" (sano)
                # también en el siguiente ciclo. Esto simula la recuperación natural tras el diagnóstico.
                if (any(to_H)) {
                  m_M[recovered_rows[to_H], t + 1] <- "H"
                }
                
                # Registramos la recuperación (sin coste) 
                # en el log para análisis posterior
                # Note: age is incremented by 1 because recovery occurs after state transition
                cost_log <- rbindlist(list(cost_log, data.table(
                  sim = current_sim,
                  age = age_in_loop + 1,
                  ID = recovered_ids,
                  cost_type = paste0("recovery_from_", recovered_states),
                  cost = 0
                )), use.names = TRUE)
                
                # Remove recovered from follow-up if present
                CIN1_followup_IDs <- setdiff(CIN1_followup_IDs, recovered_ids)
              }
            }
          }
        }
        
        # --------------------------------------------------------------------
        # -------------------- CIN1 Follow-Up Block 1.4 -----------------------
        # --------------------------------------------------------------------
        if (length(CIN1_followup_IDs) > 0) {
          followup_rows <- match(CIN1_followup_IDs, IDs)
          
          # Step 1: Follow-up for still CIN1 and not already detected
          # (Seguimiento si siguen siendo CIN1 y no están ya detectados)
          current_states <- m_M[followup_rows, t]
          still_CIN1 <- current_states == "CIN1"
          still_CIN1_IDs <- CIN1_followup_IDs[still_CIN1]
          still_CIN1_IDs <- setdiff(still_CIN1_IDs, detected_IDs)
          
          if (length(still_CIN1_IDs) > 0) {
            cost_log <- rbindlist(list(cost_log, data.table(
              sim = current_sim,
              age = age_in_loop,
              ID = still_CIN1_IDs,
              cost_type = "CIN1_followup",
              cost = costCoeff_md[match("CIN1", v_n)]
            )), use.names = TRUE)
          }
          
          # Step 2: Progression to CIN2+ at t + 1
          # Check if any CIN1 cases progressed to CIN2+ or cancer
          # Detectamos progresión a CIN2+ en el estado del próximo ciclo
          next_states <- m_M[followup_rows, t + 1]
          progressed <- next_states %in% c("CIN2", "CIN3", "FIGO.I", "FIGO.II", "FIGO.III", "FIGO.IV")
          
          if (any(progressed)) {
            progressed_IDs <- CIN1_followup_IDs[progressed]
            # Only new progressions (not already treated)
            progressed_IDs <- setdiff(progressed_IDs, detected_IDs) # Evitar duplicados
            
            if (length(progressed_IDs) > 0) {
              progressed_states <- next_states[match(progressed_IDs, CIN1_followup_IDs)]
              progressed_indices <- match(progressed_states, v_n)
              
              # Coste de tratamiento tras progresión desde CIN1
              cost_log <- rbindlist(list(cost_log, data.table(
                sim = current_sim,
                age = age_in_loop + 1, # tratamiento se registra en el ciclo siguiente
                ID = progressed_IDs,
                cost_type = paste0("progressed_from_CIN1_", progressed_states),
                cost = costCoeff_md[progressed_indices]
              )), use.names = TRUE)
              
              detected_IDs <- unique(c(detected_IDs, progressed_IDs))
              CIN1_followup_IDs <- setdiff(CIN1_followup_IDs, progressed_IDs)
            }
          }
          
          # Step 3: Regressed -> Stop follow-up
          # Regresión natural — si vuelven a "H" o "HPV.infection", se detiene el seguimiento
          regressed <- next_states %in% c("H", "HPV.infection")
          if (any(regressed)) {
            regressed_IDs <- CIN1_followup_IDs[regressed]
            CIN1_followup_IDs <- setdiff(CIN1_followup_IDs, regressed_IDs)
          }
        }
        
        
         
        
        ## Comparison table (KEEP THIS COMMENTED WHEN RUNNING)
        #cost_comparison <- map_dfr(names(sim_result), function(strat) {
        #  micro_costs <- sim_result[[strat]]$tc_hat_undisc$tc_hat_undisc
        #  markov_cost <- sim_result[[strat]]$markov_cost_undi
        #  
        #  tibble(
        #    strategy = strat,
        #    mean_microsim_cost = mean(micro_costs),
        #    sd_microsim_cost = sd(micro_costs),
        #    markov_cost = markov_cost,
        #    difference = mean(micro_costs) - markov_cost,
        #    percent_diff = 100 * (mean(micro_costs) - markov_cost) / markov_cost
        #  )
        #})
        #
        #print(cost_comparison %>% arrange(desc(abs(percent_diff))))
        
        
      }
      #################### close loop for cycles ############################# 
      ########################################################################
      
     
      cat("----- DEBUG: cost_log for ID == 1246 -----\n")
      print(cost_log[ID == "1246"])
      cat("----- END DEBUG -----\n")
      
      
       
      # ---- UPDATE m_C MATRIX WITH CURRENT COST_LOG ENTRIES ----
      ########################################################################
      # data.table version (pick only one):
      # Versión usando data.table para eficiencia.
      # Este bloque se ejecuta después de terminar el loop (n_t),
      # y agrega a la matriz m_C los costes acumulados durante el ciclo.
      if (nrow(cost_log) > 0) {
        # Asignamos los índices de fila y columna para m_C:
        # - row_i: fila correspondiente al ID del individuo (convertido a entero)
        # - col_t: columna correspondiente al ciclo temporal (edad - 9)
        #          asumiendo que la edad mínima es 9 y corresponde a la columna 1
        temp_cost_log <- copy(cost_log)
        temp_cost_log[, row_i := as.integer(ID)]
        temp_cost_log[, col_t := age - 9]
        
        # Filter invalid indices
        # Filtramos los índices que están fuera de los límites de m_C
        temp_cost_log <- temp_cost_log[row_i >= 1 & row_i <= nrow(m_C) &
                                         col_t >= 1 & col_t <= ncol(m_C)]
        
        # Aggregate costs by row_i and col_t
        agg_costs <- temp_cost_log[, .(total_cost = sum(cost)), by = .(row_i, col_t)]
        
        # Update m_C with aggregated costs:
        for (i in seq_len(nrow(agg_costs))) {
          m_C[agg_costs$row_i[i], agg_costs$col_t[i]] <- 
            m_C[agg_costs$row_i[i], agg_costs$col_t[i]] + agg_costs$total_cost[i]
        }
      }
      ########################################################################
      #########################################################################
      ## dplyr version (pick only one):
      #if (nrow(cost_log) > 0) {
      # cost_log_df <- as.data.frame(cost_log)
      # 
      # row_i <- as.integer(cost_log_df$ID)
      # col_t <- cost_log_df$age - 9
      # 
      # valid_idx <- row_i >= 1 & row_i <= nrow(m_C) &
      #   col_t >= 1 & col_t <= ncol(m_C)
      # 
      # agg_df <- tibble(row = row_i[valid_idx],
      #                  col = col_t[valid_idx],
      #                  cost = cost_log_df$cost[valid_idx]) %>%
      #   group_by(row, col) %>%
      #   summarise(cost = sum(cost), .groups = "drop")
      # 
      # for (i in seq_len(nrow(agg_df))) {
      #   m_C[agg_df$row[i], agg_df$col[i]] <-
      #     m_C[agg_df$row[i], agg_df$col[i]] + agg_df$cost[i]
      # }
      #}
      #########################################################################
      
      # reset cyto_screen_days
      cyto_screening_days <- NULL
      
      
      # Combine stored entries in a single data frame
      symptomatics <- bind_rows(stored_list)
      cat("=================================================\n")
      cat("symprotamics dimensions: ", dim(symptomatics),"\n")
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
      
      Tot_Trans_per_t <- Tot_Trans_per_t %>% as_tibble()
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
      
      #NEW CODE 26.06.25:
      cost_log <- cost_log[!duplicated(cost_log[, .(sim, age, ID, cost_type)]), ]
      
      #cost_log <- unique(cost_log, by = c("sim", "age", "ID", "cost_type"))
      
      #cat("At sim number:", sim,  " reported strategy is ", strategy, "\n")
      
      # Store the results from the simulation in a list
      results <- list(#strategy = strategy,
        seed = seeds[sim],
        #seed = seed,
        #sim_numb = sim, 
        m_M = m_M, 
        m_C = m_C, 
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
        CC_Death_by_diff = CC_Death_by_diff, 
        screening_cost = cost_log) 
      
      results$seed <- seeds[sim]
      #results$seed <- seed
      #simulation_results[sim] <- list(results)
      #simulation_results[sim] <- results
      cat("At sim number:", sim,  " tc_hat_undisc is ", tc_hat_undisc, "\n")
      rm(symptomatics)
      #rm(TS) 
      
      ## Write to a log file to track worker outputs
      #cat(sprintf("Simulation %d, Length: %d\n", sim, length(output)), 
      #    file = "debug_log.txt", append = TRUE)
      return(results)
      #gc() #Force memory cleanup after each sim/batch 
      
    } # end of `foreach/dopar` loop
  
  return(simulation_results)
  #cat("Lenght of simulation_results = ", length(simulation_results), "\n")
  #
  #stacked_results <- 
  #  summarize_results_by_Strategy_v2(
  #    #strategy = screening_strategies[[strat]]$sim.name,
  #    strategy = strat,
  #    results_list = simulation_results, 
  #    numb_of_sims = numb_of_sims)
  
  #joined_batches_per_strategy[[strat]] <-  stacked_results
  #return(simulation_results)
  
  return(stacked_results)
  #return(joined_batches_per_strategy)
  
} # end of MicroSim function
################################################################################
################################################################################


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
## Vaccination strategies:
# Paramters:
# vaccination coverage for vacc 2, 4 and 9:

vacc_coverage <- c(0.0, 0.0, 0.0) 

# natural immunity associated with vacc 2, 4, and 9:
nat_immunity_linked_to_vacc <- c(0.0, 0.0, 0.0)

################################################################################
# RANDOM FUNCTION
generate_vaccine_labels <- function(n_i, vacc_coverage, nat_immunity, seed) {
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
 
  # RANDOM GEN line: 
  # Shuffle the vector randomly
  seed_2 <- 123
  set.seed(seed_2)
  vacc_lbl <- sample(vacc_lbl, size = n_i, replace = FALSE)
  
  # Initialize the 'immuned' vector with FALSE for everyone
  immuned <- rep(FALSE, n_i)
  
  # For vaccinated individuals, check if they overcome their immunity probability:
  for (i in 1:n_i) {
    if (vacc_lbl[i] == "vacc_2") {
      # Check if individual overcomes immunity probability for vacc_2
      immuned[i] <- runif(1) < nat_immunity[1]
    } else if (vacc_lbl[i] == "vacc_4") {
      # Check if individual overcomes immunity probability for vacc_4
      immuned[i] <- runif(1) < nat_immunity[2]
    } else if (vacc_lbl[i] == "vacc_9") {
      # Check if individual overcomes immunity probability for vacc_9
      immuned[i] <- runif(1) < nat_immunity[3]
    }
    # Individuals with "no_vacc" remain FALSE for immunity
  }
  
  # Create the final data frame with ID, vacc_state, and immuned status
  result_df <- data.frame(
    ID = seq_len(n_i),
    vacc_state = vacc_lbl,
    immuned = immuned
  )
  
  return(result_df)
}
################################################################################
# RANDOM FUNCTION:
vacc_lbl <-
  generate_vaccine_labels(n_i, vacc_coverage, nat_immunity_linked_to_vacc, seed)
################################################################################


################################################################################
################################################################################
########################## Run the simulation ##################################
##  START SIMULATION
Sys.setenv(OMP_NUM_THREADS = "1") # to prevent conflicts between OpenMP and R parallel
p = Sys.time()
numb_of_sims = 3
#numb_of_sims =  4
#numb_of_sims = 1
#numb_of_sims = 20

# Initialize individual IDs
IDs <- 1:n_i
#strategy <- "natural_history"
#strategy <- "vacc_2_coverage_0.0"

# Screening Strategies:
source(file = "R/params_only_cyto_AMontoliu.R") # load Parameters_strategies()
screening_coverage = 0.8; vacc_coverage = 0
ScreenPrice.md = ScreenPrice = 27.86
# Direct medical costs of monitoring and treatment in each state:
costCoeff_md <- c(0, 39.54, 288.91, 1552.27, 1552.27, 5759.81,
                   12903.63, 23032.41, 35323.14, 0, 0, 0)

# Direct non-medical costs of monitoring and treatment in each state:
costCoeff_nmd <-  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# Indirect costs of monitoring and treatment in each state:
costCoeff_i <-  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# Sensitivity of cytology as a primary test for each state:
screenSensi <- c(0, 0, 0.177, 0.5, 0.523, 1, 1, 1, 1, 0, 0, 0)

# Cytology Specificity:
citoSpecif <- 0

screening_strategies <- Parameters_strategy(Coverage = screening_coverage, 
                                            cobertura_vacuna = vacc_coverage)
##Test:
#screening_strategy_2 <- screening_strategies[2]
#screening_strategies <- screening_strategy_2

figoSymProb <- c(0.11, 0.23, 0.66, 0.9) 

# prob of recovery during cytology screening:
screenProbs <- c(0, 0, 1, 1, 1, 0.9688, 0.9066, 0.7064, 0.3986, 0, 0, 0)
### FOR TESTING 
#screenProbs <- c(0, 0, .5, .5, .5, 0.9688, 0.9066, 0.7064, 0.3986, 0, 0, 0)
#screenProbs <- c(0, 0, 0, 0, 0, 0.9688, 0.9066, 0.7064, 0.3986, 0, 0, 0)
symptom_prob_vec <- figoSymProb
survival_prob_vec <- screenProbs[6:9]
states_to_check <- c("FIGO.I", "FIGO.II", "FIGO.III", "FIGO.IV")

stored_list <- vector("list", n_t)



# Init Storage for Costs Output
cost_log <- 
  data.table(sim = integer(), 
             age = integer(), 
             ID = integer(), 
             cost_type = character(),
             cost = numeric())

numb_screening_strat <- screening_strategies %>% length()

sim_result <- list()
source("./R/sumarize_results_by_Strategy_Func.R")

for (n_strat in 1:length(screening_strategies)) {
  strat <- screening_strategies[[n_strat]]$sim.name
  
  cat("#######################################################################\n")
  cat ("The Strategy is ", strat, "\n")
  cat ("n_strat is ", n_strat, "\n")
  cat("#######################################################################\n")
  
  #sim_result[[strat]] <- MicroSim(strat = strat, 
  sim_raw_result <- MicroSim(strat = strat, 
                             numb_of_sims = numb_of_sims, 
                             v_M_1 = v_M_1,
                             n_i = n_i, 
                             n_t = n_t, 
                             v_n = v_n, 
                             d_c = d_c, 
                             d_e = d_e, 
                             TR_out = TRUE, 
                             TS_out = TRUE, 
                             Trt = FALSE, 
                             Pmatrix = Pmatrix,
                             master_seed = 123,
                             reproducible = TRUE, 
                             use_parallel = FALSE,
                             #use_parallel = TRUE,
                             cost_vacc2, 
                             cost_vacc4, 
                             cost_vacc9,
                             screening_coverage = screening_coverage,
                             vacc_coverage = vacc_coverage,
                             ScreenPrice = ScreenPrice,
                             costCoeff_md = costCoeff_md,
                             citoSpecif = citoSpecif 
  )
  
  # Check MicroSim output
  if (is.null(sim_raw_result)) {
    warning(paste("MicroSim returned NULL for strategy:", strat))
    next
  }
  
  # Summarize results
  stacked_results <- summarize_results_by_Strategy(strategy = strat,
                                                   results_list = sim_raw_result, 
                                                   numb_of_sims = numb_of_sims)
  
  # Save
  sim_result[[strat]] <- stacked_results 
  #sim_result[[strat]] <- stacked_results[[strat]] 
}

comp.time = Sys.time() - p
comp.time %>% print()





################################################################################
# POST-PROCESSING:
source(file = "R/post_process_strategies.R", local = environment())
################################################################################


################################################################################
# ADD MARKOV RESULTS
source(file = "R/add_markov_results.R")

################################################################################

cat("Markov results added, stop here for the moment\n")
stop()


#################################################################################  
### ----Incidences, Prevalences, and Mortalities
## ADDING MARKOV RESULTS (corrected):
##load(file = "data/markov_results/markov_vacc_incidences_vectors.RData")
## Loading the CORRECTED-TRANSITIONs results, the 'markov_sim_vacc_incidences' object:
##load(file = "data/markov_results/markov_vacc_CORRECTED_incidences_vectors.RData")
#load(file = "data/markov_results/markov_vacc_CORRECTED_incidences_vectors_20250425.RData")
## Markov:
##markov_CN1_incidences  <- c(0.00000, 204.73492, 981.96179, 1368.24200, 3006.85782, 33.48096, 1362.96678, 459.48051, 697.84223, 794.33833, 223.00222, 246.23082, 176.02167, 126.22963, 53.70939)
#sim_result[[1]]$markov_CN1_incidences  <- markov_sim_vacc_incidences[[paste0("Markov_CIN1_Incidence_", vacc_coverage[1]*10^2)]] 
#
##markov_CN2_incidences  <- c(0.000000, 6.165629, 54.767952, 140.309815, 216.568392, 1476.306267, 1579.728160, 1298.914564, 466.596151, 637.661611, 442.298632, 304.784447, 250.953880, 165.628020, 116.925192)
#sim_result[[1]]$markov_CN2_incidences  <- markov_sim_vacc_incidences[[paste0("Markov_CIN2_Incidence_", vacc_coverage[1]*10^2)]]
#
##markov_CN3_incidences  <- c(0.000000, 2.090325, 9.597415, 44.467676, 148.972191, 0.000000, 3.550684, 91.881726, 12.505042, 68.377446, 25.802481, 7.952667, 1.174088, 1.177840, 2.638642)
#sim_result[[1]]$markov_CN3_incidences  <- markov_sim_vacc_incidences[[paste0("Markov_CIN3_Incidence_", vacc_coverage[1]*10^2)]]
#
##markov_CC_incidences   <- c(0.000000, 0.000000, 0.000000, 5.520938, 8.360544, 13.282380, 22.906871, 20.825560, 15.867891, 32.483846, 8.962389, 17.681771, 11.737615, 17.354646, 14.582775)
#sim_result[[1]]$markov_CC_incidences   <- markov_sim_vacc_incidences[[paste0("Markov_CC_Incidence_", vacc_coverage[1]*10^2)]]
#
##markov_HPV_prevalences <- c(0.000000000, 0.343480414, 0.377634762, 0.087223460, 0.307341403, 0.030196332, 0.050562845, 0.050151668, 0.082952596, 0.046644059, 0.018532077, 0.034193076, 0.016407832, 0.015039027, 0.003217326)
#sim_result[[1]]$markov_HPV_prevalences <- markov_sim_vacc_incidences[[paste0("Markov_HPVPrevalence_", vacc_coverage[1]*10^2)]]
#
## markov_CC_mortality <- c(0.000000e+00, 0.000000e+00, 0.000000e+00, 2.977975e-06, 
##                          1.574920e-05, 2.715056e-05, 5.489929e-05, 7.284815e-05,
##                          1.057494e-04, 5.076268e-05, 7.517773e-05, 4.960943e-05,
##                          4.802468e-05, 4.210457e-05, 4.837655e-05) * 10^5
#sim_result[[1]]$markov_CC_mortality <- markov_sim_vacc_incidences[[paste0("Markov_CCMortality_", vacc_coverage[1]*10^2)]]
#
#
## NOTE: change for corresponding vacc strategy 0, 60, 70, or 80:
#sim_result[[1]]$markov_new_CIN1   <- markov_sim_vacc_incidences[[paste0("Markov_n CIN1_", vacc_coverage[1]*10^2)]]
#sim_result[[1]]$markov_new_CIN2   <- markov_sim_vacc_incidences[[paste0("Markov_n CIN2_", vacc_coverage[1]*10^2)]]
#sim_result[[1]]$markov_new_CIN3   <- markov_sim_vacc_incidences[[paste0("Markov_n CIN3_", vacc_coverage[1]*10^2)]]
#sim_result[[1]]$markov_new_Cancer <- markov_sim_vacc_incidences[[paste0("Markov_n CC_", vacc_coverage[1]*10^2)]]
#################################################################################  
#
#################################################################################  
### Adding corresponding (to vacc strategy) Markov QUALYs and Costs to sim_result
#master_markov_vacc_results_CORRECTED <- 
#  get(load(file = "data/corrected_transitions_20250414/Markov_results_only_vaccination_Strategies_only_vac_20250414.rda"))
## Select vaccination-associated QALYs, Costs (discounted and undiscounted):
#rm(df)
#
## SELECT VACCINATION LEVEL:
##markov_vacc_lvl <- 80 # it can be 0, 60, 70 or 80
#markov_vacc_lvl <- vacc_coverage[1]*10^2 # it can be 0, 60, 70 or 80
#
#if(markov_vacc_lvl == 0) {
#  Mark_vacc_lvl <- "Vaccination coverage: 0%"
#} else if (markov_vacc_lvl == 60) {
#  Mark_vacc_lvl <- "Vaccination coverage: 60%"
#} else if (markov_vacc_lvl == 70) {
#  Mark_vacc_lvl <- "Vaccination coverage: 70%"
#} else if (markov_vacc_lvl == 80) {
#  Mark_vacc_lvl <- "Vaccination coverage: 80%"
#} else {
#  cat("Not valid Markov vacc level")
#}
#
#sim_result[[1]]$markov_qaly_undis <- master_markov_vacc_results_CORRECTED %>% 
#  dplyr::filter(sim.name == Mark_vacc_lvl) %>%
#  dplyr::select(`Per person QALYs und`) %>% 
#  as.numeric()
#
#sim_result[[1]]$markov_qaly_undis_Tot <- master_markov_vacc_results_CORRECTED %>% 
#  dplyr::filter(sim.name == Mark_vacc_lvl) %>%
#  dplyr::select(`Total QALYs und`) %>% 
#  as.numeric()
#
#sim_result[[1]]$markov_qaly_disc <- master_markov_vacc_results_CORRECTED %>% 
#  dplyr::filter(sim.name == Mark_vacc_lvl) %>%
#  dplyr::select(`Per person QALYs disc`) %>% 
#  as.numeric()
#
#sim_result[[1]]$markov_qaly_disc_Tot <- master_markov_vacc_results_CORRECTED %>% 
#  dplyr::filter(sim.name == Mark_vacc_lvl) %>%
#  dplyr::select(`Total QALYs disc`) %>% 
#  as.numeric()
#
### -- ##
#
#sim_result[[1]]$markov_cost_undis <- master_markov_vacc_results_CORRECTED %>% 
#  dplyr::filter(sim.name == Mark_vacc_lvl) %>%
#  dplyr::select(`Per person D cost und`) %>% 
#  as.numeric()
#
#sim_result[[1]]$markov_cost_undis_Tot <- master_markov_vacc_results_CORRECTED %>% 
#  dplyr::filter(sim.name == Mark_vacc_lvl) %>%
#  dplyr::select(`Total D cost und`) %>% 
#  as.numeric()
#
#sim_result[[1]]$markov_cost_disc <- master_markov_vacc_results_CORRECTED %>% 
#  dplyr::filter(sim.name == Mark_vacc_lvl) %>%
#  dplyr::select(`Per person D cost disc`) %>% 
#  as.numeric()
#
#sim_result[[1]]$markov_cost_disc_Tot <- master_markov_vacc_results_CORRECTED %>% 
#  dplyr::filter(sim.name == Mark_vacc_lvl) %>%
#  dplyr::select(`Total D cost disc`) %>% 
#  as.numeric()
#################################################################################  
#################################################################################  
  
################################################################################  
## Adding corresponding (to vacc strategy) Markov age-averaged new_cases to sim_result
library(tibble)
library(stringr)
library(purrr)

# Define your state suffixes
states <- c("CIN1", "CIN2", "CIN3", "Cancer")

# Ensure sim_result exists and has an entry at [[1]]
if (!exists("sim_result")) sim_result <- list()
if (is.null(sim_result[[1]])) sim_result[[1]] <- list()

# Helper function to convert vector to tibble
convert_markov_vector <- function(vec) {
  tibble(
    age_interval = names(vec) %>%
      str_extract("[0-9]{2}-[0-9]{2}") %>%
      factor(levels = unique(.)),
    mean_new_cases = as.numeric(vec)
  )
}

## Convert and assign each result into sim_result[[1]]
purrr::walk(states, function(state) {
  obj_name <- paste0("markov_new_", state)  # Only the field name
  new_name <- paste0("new_averaged_", state, "_Markov_per_age_interval")
  
  # access directly inside sim_result[[1]], no get()
  vec <- sim_result[[1]][[obj_name]]
  
  # assign the converted result back into sim_result[[1]]
  sim_result[[1]][[new_name]] <<- convert_markov_vector(vec)
})

################################################################################ 

################################################################################  
################################################################################  
################################################################################  
## Adding MicroSim results:
microSim_CN1_incidences          <- sim_result[[1]]$mean_incidence_CIN1_per_age_interval
microSim_CN2_incidences          <- sim_result[[1]]$mean_incidence_CIN2_per_age_interval
microSim_CN3_incidences          <- sim_result[[1]]$mean_incidence_CIN3_per_age_interval
microSim_CC_incidences           <- sim_result[[1]]$mean_CC_incidence
microSim_HPV_prevalences         <- sim_result[[1]]$mean_HPV_prevalence_per_age_interval
microSim_CC_mortality            <- sim_result[[1]]$CC_mean_mortality
microSim_CC_by_diff_mortality    <- sim_result[[1]]$CC_by_diff_mean_mortality
microSim_new_CIN1                <- sim_result[[1]]$new_averaged_CIN1_per_age_interval
microSim_new_CIN2                <- sim_result[[1]]$new_averaged_CIN2_per_age_interval
microSim_new_CIN3                <- sim_result[[1]]$new_averaged_CIN3_per_age_interval
microSim_new_Cancer              <- sim_result[[1]]$new_averaged_Cancer_per_age_interval
################################################################################  
################################################################################  


cat("Hey, I'm done, and about to write out the results\n")

args <- commandArgs(trailingOnly = TRUE)

#if (length(args) == 0) {
#  stop("Output tag argument missing. Usage: Rscript run_sim.R <unique_tag>")
#}
#
#unique_tag <- args[1]

# Use passed SLURM job ID as output tag
slurm_job_id <- if (length(args) >= 1) args[1] else NA

if (is.na(slurm_job_id) || slurm_job_id == "") {
  slurm_job_id <- format(Sys.time(), "%Y%m%d%H%M%S")
  cat("Warning: SLURM_JOB_ID not provided. Using timestamp fallback:", slurm_job_id, "\n")
} else {
  cat("SLURM job ID:", slurm_job_id, "\n")
}


## Extract first vaccine coverage value for filename
#vacc_tag <- sprintf("%.1f", vacc_coverage[1])  # Format as 0.8, 0.0, etc.

## Define directory and static filename components
##output_dir <- "data/TESTING_20250429"
#output_dir <- "data/cyto_screening/"
##base_filename <- "stacked_sims_20x10E6x75_vacc2_0.8_NEW_TRANSITIONS_PARA_20250506_sim_"
##base_filename <- paste0("stacked_sims_20x10E6x75_vacc2_", vacc_tag,
##                        "_update_WITH_select_floorswitch_NEW_TRANSITIONS_PARA_20250522_A_sim_")
#base_filename <- paste0("cyto_screening_sims_20x10E6x75_coverage_", screening_coverage,
#                        "_recovery_CIN123_", screenProbs[3])
#
### Construct full path
##output_file <- file.path(output_dir, paste0(base_filename, unique_tag, ".rds"))
#
#output_file <- file.path(output_dir, paste0(base_filename, slurm_job_id, ".rds"))
#
#
#cat("Saving simulation result to:", output_file, "\n")
#saveRDS(object = sim_result, file = output_file)


################################################################################
################################################################################
################################################################################

# LOAD SIMULATION:
##### TO LOAD PRE-RUN MICRO-SIMULATIONS:
#sim_result <- readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_vacc2_0.0_NEW_TRANSITIONS_PARA_20250425_sim_20495.rds")
### IF YOU LOAD A PRE-RUN SIMULATION AND WANT TO POST-PROCESS RUN SCRIPT FROM HERE ALL
### WAY TO THE BOTTOM AND INCLUDE THE NEEDED FOLLOWING VARIABLES:
### (IF NOT LOADING PRE-RUNNED SIMULATION KEEP THE FOLLOWING 4 LINES COMMENTED)
##rm(list = ls())
#n_i <- 10^6
#numb_of_sims <- 20
#n_t = 75
#vacc_coverage <- c(0.0,0,0) # for correct plot titles 

# Load microsim results with vacc strategies (sequentially runned):
sim_result_0 <- 
  readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_vacc2_0.0_SEQ_20250416_sim_20015.rds")

sim_result_60 <- 
  readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_vacc2_0.6_SEQ_20250414_sim_19942.rds")

sim_result_70 <- 
  readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_vacc2_0.7_SEQ_20250415_sim_19962.rds")

sim_result_80 <- 
  readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_vacc2_0.8_SEQ_20250414_sim_19941.rds")

## Load microsim results with vacc strategies (parallel runned):
#sim_result_0 <- 
#  readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_vacc2_0.0_NEW_TRANSITIONS_PARA_20250425_sim_20495.rds")
#sim_result_0 <- 
#  readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_vacc2_0.0_NEW_TRANSITIONS_PARA_20250424_sim_20379.rds")
#
#sim_result_60 <- 
#  readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_vacc2_0.6_PARA_20250416_sim_20032.rds")
#
#sim_result_70 <- 
#  readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_vacc2_0.7_PARA_20250416_sim_20031.rds")
#
#sim_result_80 <- 
#  readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_vacc2_0.8_PARA_20250416_sim_20030.rds")

# Load old natural history (no vaccination strategy implmented so it should match vacc_0.0 strategy):
#sim_natural_history <- 
#  readRDS(file = "data/last_results_20250324/stacked_sims_PARA_20x10E6x75_20250417_NATURAL_HISTORY_sim_20042.rds")

# Old transitions:
#sim_natural_history <- 
#  readRDS(file = "data/last_results_20250324/stacked_sims_PARA_20x10E6x75_20250424_NATURAL_HISTORY_REPROD_sim_20324.rds")
#sim_natural_history <- 
#  readRDS(file = "data/natural_history/stacked_sims_20x10E6x75_20250211_madeinPADO_PARA_NATURAL_HISTORY7071.rds")

# Load result for vacc = 0 and using the old transitiosn for debuging purposes:
sim_result_0_old_trans <- 
  readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_vacc2_0.0_OLD_TRANSITIONS_PARA_20250424_sim_20330.rds")

# load  Markov vaccination computation strategies:
#load(file = "data/markov_vacc_vectors.RData")
load(file = "data/markov_vacc_CORRECTED_vectors.RData")

## Use this only to pots-process some of the previous results:

#sim_result <- sim_result_0
#sim_result <- sim_natural_history
#sim_result <- sim_result_0_old_trans


cat("I have written out the results\n")

### ----Convert .Rmd to .R
#library(knitr)
## purl("your_script.Rmd", output = "your_script.R")
## example:
#purl("Cervix_MicroSim_RMarkdown_v.072_B.Rmd", output = "cervix_microSim_stacked_list.R")

## ---- COST-EFECTIVENES
####################### Cost-effectiveness analysis #############################
## store the mean costs (and MCSE) of each strategy in a new variable C (vector costs)
#v_C  <- c(sim_result$tc_hat_disc, sim_trt$tc_hat_disc) 
#sd_C <- c(sd(sim_result$tc_disc), sd(sim_trt$tc_disc)) / sqrt(n_i)
## store the mean QALYs (and MCSE) of each strategy in a new variable E (vector effects)
#v_E  <- c(sim_result$te_hat_disc, sim_trt$te_hat_disc)
#sd_E <- c(sd(sim_result$te_disc), sd(sim_trt$te_disc)) / sqrt(n_i)
#
#delta_C <- v_C[2] - v_C[1]                   # calculate incremental costs
#delta_E <- v_E[2] - v_E[1]                   # calculate incremental QALYs
## Monte Carlo Squared Error (MCSE) of incremental costs:
#sd_delta_E <- sd(sim_trt$te - sim_result$te) / sqrt(n_i) 
## Monte Carlo Squared Error (MCSE) of incremental QALYs:
#sd_delta_C <- sd(sim_trt$tc_disc - sim_result$tc_disc) / sqrt(n_i) 
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

################################################################################
################################################################################
###                      PLOTTING ROUTINES                                    ##
################################################################################
################################################################################
## ---- Plot curves
## This R chunk is a plot routine (not part of the main program):
library(RColorBrewer)
#ensure_library("RColorBrewer")
# Convert matrix to data frame
#micro_sim_df <- sim_result[[1]]$TR
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

## ---- Loading Markov result
if (!require("readxl")) install.packages("readxl")
library(readxl)
# This R chunk is a plot routine (not part of the main program):


if (!require("readxl")) install.packages("readxl")
library(readxl)
#markov <-
# readxl::read_excel("Q:/my_Q_docs/Cervix_MicroSim/CervixMicroSim_Carlos/carlos__Krijkamp_ver/data/Sortida_NoIntervencio.xlsx", sheet = "NH")
markov_df <- readxl::read_excel("./data/Sortida_NoIntervencio.xlsx")

markov_df <- markov_df %>% mutate(age = Step + 10)

# Reshape the data into long format
markov_df_long_data <- markov_df %>%
  pivot_longer(cols = c(HR.HPV.infection, CIN1, CIN2, CIN3, FIGO.I, FIGO.II,
                        #FIGO.III, FIGO.IV, Survival, CC_Death, Other.Death),
                        FIGO.III, FIGO.IV, Survival, CC_Death),
               names_to = "Health state",
               values_to = "value")

## Plot the data
#ggplot(markov_df_long_data, aes(x = age, y = value, color = `Health state`)) +
#  geom_line(linewidth=1, alpha=0.7) +
#  labs(x = "Age", y = "Value", color = "Health state") +
#  ggtitle(expression(paste("Markov cohort simulation for ", 10^6, " individuals"))) + 
#  theme_minimal()  # Optional: customize the theme

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


################################################################################  
################################################################################  
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
  markov_CN1_incidences = sim_result[[1]]$markov_CN1_incidences,
  markov_CN2_incidences = sim_result[[1]]$markov_CN2_incidences,
  markov_CN3_incidences = sim_result[[1]]$markov_CN3_incidences,
  markov_CC_incidences = sim_result[[1]]$markov_CC_incidences,
  markov_HPV_prevalences = sim_result[[1]]$markov_HPV_prevalences,
  markov_CC_mortality = sim_result[[1]]$markov_CC_mortality,
  # Assign the same values from markov_CC_mortality to markov_CC_by_diff_mortality
  markov_CC_by_diff_mortality <- sim_result[[1]]$markov_CC_mortality
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

# Conversion if you have microSim data as tibbles
microSim_data <- data.frame(
  age = age_groups,
  
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

################################################################################
plot_comparison <- function(data, measure_name) {
  ggplot(data %>% dplyr::filter(grepl(measure_name, measure)), 
         aes(x = age, y = value, fill = model)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = paste(
        measure_name, "Comparison\n",
        "N =", sim_result[[1]]$numb_of_ind, 
        ";  cycles=", sim_result[[1]]$numb_of_cycles, 
        "Para.", 
        "Avgd. sims =", numb_of_sims, "\n",
        "Vacc.=", vacc_coverage[1]
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

################################################################################
## Plotting FIGO prevalences
figo_data_prevalence <- sim_result[[1]]$mean_FIGO_prevalence

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
mean_FIGO <- sim_result[[1]]$mean_FIGO

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
mean_Diagnosed <- sim_result[[1]]$mean_Diagnosed

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
## Plot new individuals in epi classes averaged by age interval:
plot_mean_new_CIN1 <-
  #ggplot(microSim_new_CIN1, aes(x = age_interval, y = mean_new_cases)) +
  ggplot(sim_result[[1]]$new_averaged_CIN1_per_age_interval, aes(x = age_interval, y = mean_new_cases)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Average New CIN1 Cases by Age Group",
       x = "Age Group", y = "Mean New Cases") + 
  theme_minimal(base_size = 14) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_mean_new_CIN2 <-
  #ggplot(microSim_new_CIN2, aes(x = age_interval, y = mean_new_cases)) +
  ggplot(sim_result[[1]]$new_averaged_CIN2_per_age_interval, aes(x = age_interval, y = mean_new_cases)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Average New CIN2 Cases by Age Group",
       x = "Age Group", y = "Mean New Cases") + 
  theme_minimal(base_size = 14) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plot_mean_new_CIN3 <-
  #ggplot(microSim_new_CIN3, aes(x = age_interval, y = mean_new_cases)) +
  ggplot(sim_result[[1]]$new_averaged_CIN3_per_age_interval, aes(x = age_interval, y = mean_new_cases)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Average New CIN3 Cases by Age Group",
       x = "Age Group", y = "Mean New Cases") + 
  theme_minimal(base_size = 14) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_mean_new_Cancer <-
  #ggplot(microSim_new_Cancer, aes(x = age_interval, y = mean_new_cases)) +
  ggplot(sim_result[[1]]$new_averaged_Cancer_per_age_interval, aes(x = age_interval, y = mean_new_cases)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Average New Cancer Cases by Age Group",
       x = "Age Group", y = "Mean New Cases") + 
  theme_minimal(base_size = 14) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_mean_new_CC <-
  #ggplot(sim_result[[1]]$new_averaged_CC_Death_per_age_interval, aes(x = age_interval, y = mean_new_cases)) +
  ggplot(sim_result[[1]]$new_averaged_CC_Death_per_age_interval, aes(x = age_interval, y = mean_new_cases)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Average New Cancer Death Cases by Age Group",
       x = "Age Group", y = "Mean New Cases") + 
  theme_minimal(base_size = 14) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################################################################################
# Create comparison for new individuals in diverse epi clases. 
# # These are raw numbers, not normalized numbers such as incidencem and prevalences
compare_models_plot <- function(markov_vector, microsim_tbl, 
                                outcome_label = "Outcome", 
                                microsim_col = "mean_new_cases",
                                #N = "1e+06", cycles = 75, sims = 20, vacc = 0.8) {
                                N = n_i, cycles = n_t,
                                sims = numb_of_sims, 
                                vacc = vacc_coverage[1]) {
  # Clean age group labels from Markov vector
  age_labels <- sub("^n [^ ]+ ", "", names(markov_vector))
  
  # Build Markov dataframe
  markov_df <- tibble(
    age_interval = factor(age_labels, levels = unique(age_labels)),
    value = as.numeric(markov_vector),
    model = "Markov"
  )
  
  # Prepare MicroSim dataframe
  micro_df <- microsim_tbl %>%
    rename(value = all_of(microsim_col)) %>%
    mutate(model = "MicroSim",
           age_interval = factor(age_interval, levels = levels(markov_df$age_interval)))
  
  # Combine both
  combined_df <- bind_rows(markov_df, micro_df)
  
  # Plot
  ggplot(combined_df, aes(x = age_interval, y = value, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
    labs(
      title = paste("New", outcome_label, "- Comparison"),
      subtitle = paste0("N = ", N, 
                        " ; cycles= ", cycles, 
                        " ; Para. Avgd. sims = ", sims, 
                        "\nVacc.= ", vacc),
      x = "Age Group",
      y = outcome_label
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_manual(values = c("Markov" = "salmon", "MicroSim" = "turquoise3"))
}
################################################################################

plot_comparison_new_CIN1 <- compare_models_plot(
  markov_vector = sim_result[[1]]$markov_new_CIN1,
  #microsim_tbl = microSim_new_CIN1,
  microsim_tbl = sim_result[[1]]$new_averaged_CIN1_per_age_interval,
  outcome_label = "CIN1"
)

plot_comparison_new_CIN2 <- compare_models_plot(
  markov_vector = sim_result[[1]]$markov_new_CIN2,
  #microsim_tbl = microSim_new_CIN2,
  microsim_tbl = sim_result[[1]]$new_averaged_CIN2_per_age_interval,
  outcome_label = "CIN2"
)

plot_comparison_new_CIN3 <- compare_models_plot(
  markov_vector = sim_result[[1]]$markov_new_CIN3,
  #microsim_tbl = microSim_new_CIN3,
  microsim_tbl = sim_result[[1]]$new_averaged_CIN3_per_age_interval,
  outcome_label = "CIN3"
)

plot_comparison_new_Cancer <- compare_models_plot(
  markov_vector = sim_result[[1]]$markov_new_Cancer,
  #microsim_tbl = microSim_new_Cancer,
  microsim_tbl = sim_result[[1]]$new_averaged_Cancer_per_age_interval,
  outcome_label = "Cancer"
)


################################################################################
difference_plot <- function(markov_vector, microsim_tbl, 
                            outcome_label = "Outcome", 
                            microsim_col = "mean_new_cases",
                            type = c("relative", "absolute"),
                            N = n_i, cycles = n_t,
                            sims = numb_of_sims, 
                            vacc = vacc_coverage[1]) {
  
  type <- match.arg(type)
  
  # Extract age group labels
  age_labels <- sub("^n [^ ]+ ", "", names(markov_vector))
  
  # Markov data frame
  markov_df <- tibble(
    age_interval = factor(age_labels, levels = unique(age_labels)),
    markov_value = as.numeric(markov_vector)
  )
  
  # MicroSim data frame
  microsim_df <- microsim_tbl %>%
    rename(microsim_value = all_of(microsim_col)) %>%
    mutate(age_interval = factor(age_interval, levels = levels(markov_df$age_interval)))
  
  # Join and compute both differences
  diff_df <- left_join(markov_df, microsim_df, by = "age_interval") %>%
    mutate(
      absolute_difference = markov_value - microsim_value,
      relative_difference = (markov_value - microsim_value) / markov_value
    )
  
  # Choose y-axis and label
  if (type == "relative") {
    y_col <- diff_df$relative_difference
    y_label <- "Relative Difference (%)"
    y_format <- scales::percent_format(accuracy = 1)
  } else {
    y_col <- diff_df$absolute_difference
    y_label <- "Absolute Difference (Markov - MicroSim)"
    y_format <- scales::comma_format()
  }
  
  # Plot
  ggplot(diff_df, aes(x = age_interval, y = y_col, 
                      fill = factor(ifelse(y_col > 0, "Markov > MicroSim", "MicroSim > Markov"))
)) +
    geom_col() +
    scale_y_continuous(labels = y_format) +
    scale_fill_manual(
      #values = c("TRUE" = "salmon", "FALSE" = "turquoise3"),
      #labels = c("Markov > MicroSim", "MicroSim > Markov")
      values = c("Markov > MicroSim" = "salmon", "MicroSim > Markov" = "turquoise3"),
      name = "Comparison"
    ) +
    labs(
      title = paste(ifelse(type == "relative", "Relative", "Absolute"), 
                    "Difference in", outcome_label),
      subtitle = paste0(
        ifelse(type == "relative", "(Markov - MicroSim) / Markov", "Markov - MicroSim"), 
        "\nN = ", N, 
        " ; cycles = ", cycles, 
        " ; Para. Avgd. sims = ", sims, 
        "\nVacc. = ", vacc
      ),
      x = "Age Group",
      y = y_label,
      fill = "Comparison"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}
################################################################################
## Relative and absolute differences in new averaged cases

## Relative difference plot
#plot_rel_diff_CIN1 <- difference_plot(
#  markov_vector = markov_new_CIN1,
#  microsim_tbl = microSim_new_CIN1,
#  outcome_label = "CIN1",
#  type = "relative"
#)
#
## Absolute difference plot
#plot_abs_diff_CIN1 <- difference_plot(
#  markov_vector = markov_new_CIN1,
#  microsim_tbl = microSim_new_CIN1,
#  outcome_label = "CIN1",
#  type = "absolute"
#)
#
## Relative difference plot
#plot_rel_diff_CIN2 <- difference_plot(
#  markov_vector = markov_new_CIN2,
#  microsim_tbl = microSim_new_CIN2,
#  outcome_label = "CIN2",
#  type = "relative"
#)
#
## Absolute difference plot
#plot_abs_diff_CIN2 <- difference_plot(
#  markov_vector = markov_new_CIN2,
#  microsim_tbl = microSim_new_CIN2,
#  outcome_label = "CIN2",
#  type = "absolute"
#)
#
## Relative difference plot
#plot_rel_diff_CIN3 <- difference_plot(
#  markov_vector = markov_new_CIN3,
#  microsim_tbl = microSim_new_CIN3,
#  outcome_label = "CIN3",
#  type = "relative"
#)
#
## Absolute difference plot
#plot_abs_diff_CIN3 <- difference_plot(
#  markov_vector = markov_new_CIN3,
#  microsim_tbl = microSim_new_CIN3,
#  outcome_label = "CIN3",
#  type = "absolute"
#)
#
## Relative difference plot
#plot_rel_diff_Cancer <- difference_plot(
#  markov_vector = markov_new_Cancer,
#  microsim_tbl = microSim_new_Cancer,
#  outcome_label = "Cancer",
#  type = "relative"
#)
#
## Absolute difference plot
#plot_abs_diff_Cancer <- difference_plot(
#  markov_vector = markov_new_Cancer,
#  microsim_tbl = microSim_new_Cancer,
#  outcome_label = "Cancer",
#  type = "absolute"
#)


# Create plots for each measure
plot_CN1_incidences <- plot_comparison(combined_data, "CN1_incidences")
plot_CN2_incidences <- plot_comparison(combined_data, "CN2_incidences")
plot_CN3_incidences <- plot_comparison(combined_data, "CN3_incidences")
plot_CC_incidences <- plot_comparison(combined_data, "CC_incidences")
plot_HPV_prevalences <- plot_comparison(combined_data, "HPV_prevalences")
plot_CC_mortality <- plot_comparison(combined_data, "CC_mortality")
#plot_CC_by_diff_mortality <- plot_comparison(combined_data, "CC_by_diff_mortality")


## Display plots
#print(plot_CN1_incidences)
#print(plot_CN2_incidences)
#print(plot_CN3_incidences)
#print(plot_CC_incidences)
#print(plot_HPV_prevalences)
#print(plot_CC_mortality)
##print(plot_FIGO_prevalence)
##print(plot_mean_FIGO)
#print(plot_mean_Diagnosed_FIGO)
##print(plot_CC_by_diff_mortality)
#print(plot_mean_new_CIN1)
#print(plot_mean_new_CIN2)
#print(plot_mean_new_CIN3)
#print(plot_mean_new_Cancer)
#print(plot_mean_new_CC)


##### Combining plots in a single image:
library("patchwork") # disable to run as job script with sbatch:
combined_plot <- (plot_CN1_incidences | plot_CN2_incidences | plot_CN3_incidences) /
  #(plot_CC_incidences | plot_HPV_prevalences | plot_CC_mortality)
  (plot_CC_incidences | plot_HPV_prevalences | plot_mean_Diagnosed_FIGO)# plot_CC_mortality)

##combined_plot2 <- (plot_mean_new_CIN1 | plot_mean_new_CIN2 | plot_mean_new_CIN3) |
combined_plot3 <- (plot_comparison_new_CIN1 | plot_comparison_new_CIN2) /
  (plot_comparison_new_CIN3 | plot_comparison_new_Cancer) 
# View it
print(combined_plot)
print(combined_plot3)

#ggsave("figures/combined_plots_incidence_vacc_80.pdf", combined_plot, width = 15, height = 10, dpi = 300)

#################################################################################
#### Check visually for patterns in the difference of microSims and Markov sims:
#combined_plot_rel_diff_new_cases <- (plot_rel_diff_CIN1 | plot_rel_diff_CIN2) /
#                  (plot_rel_diff_CIN3 | plot_rel_diff_Cancer)  
#
#combined_plot_abs_diff_new_cases <- (plot_abs_diff_CIN1 | plot_abs_diff_CIN2) /
#                  (plot_abs_diff_CIN3 | plot_abs_diff_Cancer)  
#
#print(combined_plot_rel_diff_new_cases)
#print(combined_plot_abs_diff_new_cases)
#################################################################################


# For checking trend existence:
################################################################################
if (numb_of_sims >=60) {
  ##############################################################################
  # For number of simulations of 60 we can analize the cost results to check
  # whether there is a numerical artifact or logic code problem producing
  # a tendency of decrease tc_hat_undisc along simulations:
  average_cost <-
    sim_result[[1]]$tc_hat_undisc$tc_hat_undisc
  
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

df <- sim_result[[1]]$TR %>% select(FIGO.I, FIGO.II, FIGO.III, FIGO.IV) 
# select(FIGO.I, FIGO.II, FIGO.III, FIGO.IV) and summarize by columns
df <- df %>% summarise(across(everything(), sum, na.rm = TRUE))
################################################################################


# DEBUGGING
cat("\n")
cat("I HAVE REACHED THE END OF THE SCRIPT FINE.\n")
cat("WITH n_i = ", n_i,  " , numb_of_sims = ", numb_of_sims, "\n")

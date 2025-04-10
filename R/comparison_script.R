# Comparing vaccination startegy results  Markov vs Microsim.

rm(list = ls())

## Load a pre-run simulation:
# Natural history (code with no implementaation of vacc whatsoever)
sim_result_NH <- readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_20250401_NATURAL_HISTORY_sim_13324.rds")
# vacc 0%, 60%. 7'% amnd 80%:
sim_result_0  <- readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_20250401_NATURAL_HISTORY_sim_13324.rds")
sim_result_60 <- readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_20250401_NATURAL_HISTORY_sim_13324.rds")
sim_result_70 <- readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_20250401_NATURAL_HISTORY_sim_13324.rds")
sim_result_80 <- readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_20250401_NATURAL_HISTORY_sim_13324.rds")

# Markov sims: A=80% vaccinated, B=70% vaccinated, C=60% vaccinated:
markov_strategy_A <- load(file = "data/markov_results/results_only_vaccination_Strategy_A_paral.rda")
df_markov_strag_A <- df

markov_strategy_B <- load(file = "data/markov_results/results_only_vaccination_Strategy_B_paral.rda")
df_markov_strag_B <- df

markov_strategy_C <- load(file = "data/markov_results/results_only_vaccination_Strategy_C_paral.rda")
df_markov_strag_C <- df

# Build vectors of Markov's vaccination strategies 0%, 60%, 70%, 80% 
# of CIN1-CIN3, Cervix Cancer and HPV infection:

# Define your prefixes of interest
prefixes <- c("CIN1", "CIN2", "CIN3", "CC", "VPH")

# Define age brackets (used to construct column names)
age_from <- seq(10, 80, by = 5)
age_to   <- seq(14, 84, by = 5)
ages <- paste0(age_from, "-", age_to)

# Loop through each prefix to extract the named vectors
for (prefix in prefixes) {
  # Construct full column names for the current prefix
  cols <- paste0("n ", prefix, ages)
  
  # Extract values from the first row (vaccination strategy)
  values <- as.numeric(df_markov_strag_A[1, cols])
  
  # Assign original column names as names to the vector
  names(values) <- cols
  
  # Assign to a new variable like CIN1_Markov_80 in the global environment
  var_name <- paste0(prefix, "_Markov_80")
  assign(var_name, values)
}


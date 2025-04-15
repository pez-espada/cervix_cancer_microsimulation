# Comparing vaccination startegy results  Markov vs Microsim.

#rm(list = ls())
library(dplyr)

## Load a pre-run simulation:
# Natural history (code with no implementaation of vacc whatsoever)
sim_result_NH <- readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_20250401_NATURAL_HISTORY_sim_13324.rds")
# vacc 0%, 60%. 7'% amnd 80%:
sim_result_0  <- readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_20250401_NATURAL_HISTORY_sim_13324.rds")
sim_result_60 <- readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_20250401_NATURAL_HISTORY_sim_13324.rds")
sim_result_70 <- readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_20250401_NATURAL_HISTORY_sim_13324.rds")
sim_result_80 <- readRDS(file = "data/last_results_20250324/stacked_sims_20x10E6x75_20250401_NATURAL_HISTORY_sim_13324.rds")

# Markov sims: A=80% vaccinated, B=70% vaccinated, C=60% vaccinated:
# NOTE: these ARE NOT INCIDENCES (it appears to be new cases, but yet to be confirmed)
markov_strategy_A <- load(file = "data/markov_results/results_only_vaccination_Strategy_A_paral.rda")
df_markov_strag_A <- df

markov_strategy_B <- load(file = "data/markov_results/results_only_vaccination_Strategy_B_paral.rda")
df_markov_strag_B <- df

markov_strategy_C <- load(file = "data/markov_results/results_only_vaccination_Strategy_C_paral.rda")
df_markov_strag_C <- df
################################################################################

#################################################################################
## Preparing the Markov data for plotting
#
## List of strategies and associated row sources
#strategies <- list(
#  "60" = df_markov_strag_C,
#  "70" = df_markov_strag_B,
#  "80" = df_markov_strag_A
#)
#
## Prefixes of interest (disease stages / categories)
#prefixes <- c("CIN1", "CIN2", "CIN3", "CC", "VPH")
#
## Age brackets
#age_from <- seq(10, 80, by = 5)
#age_to   <- seq(14, 84, by = 5)
#ages <- paste0(age_from, "-", age_to)
#
## Loop over strategies and prefixes
#for (strategy in names(strategies)) {
#  df <- strategies[[strategy]]
#  
#  for (prefix in prefixes) {
#    # Construct full column names
#    cols <- paste0("n ", prefix, ages)
#    
#    # Extract values from first row (i.e., vaccination scenario)
#    values <- as.numeric(df[1, cols])
#    
#    # Add names to vector from original column names
#    names(values) <- cols
#    
#    # New name format: Markov_CIN1_60
#    var_name <- paste0("Markov_", prefix, "_", strategy)
#    
#    # Assign to global environment
#    assign(var_name, values)
#  }
#}
#################################################################################

# For the new all-included Markov vaccination strategies master df (incidences, prevalences and all)

markov_sim_vacc_w_incidences <- 
  #load(file = "data/markov_results/results_only_vaccination_Strategies_only_vac_Incidences.rda")
  # Corrected (with corrected transistions) Markov results:
  load(file = "data/corrected_transitions_20250414/Markov_results_only_vaccination_Strategies_only_vac_20250414.rda")
df_markov_all_vacc_w_incidence <- df

# Rearranging some colnames for easier data wrangling:
## First: fix the 'n CIN340-44' → 'n CIN3 40-44'
#names(df_markov_all_vacc_w_incidence) <- gsub(
#  pattern = "^(n\\s+[A-Z]+)(\\d)(\\d{2}-\\d{2})$",
#  replacement = "\\1\\2 \\3",
#  x = names(df_markov_all_vacc_w_incidence)
#)
#
## Second: capitalize patterns like 'Cin1', 'Cin2', 'Cin3' at the beginning of names
#names(df_markov_all_vacc_w_incidence) <- gsub(
#  pattern = "^Cin([123])",
#  replacement = "CIN\\1",
#  x = names(df_markov_all_vacc_w_incidence)
#)



# Pattern: any letters + digits followed by a hyphenated number group (e.g., "CIN340-44", "CC5555-1010")
names(df_markov_all_vacc_w_incidence) <- gsub(
  pattern = "(\\b[A-Z]+)(\\d+)(-\\d+\\b)",
  replacement = "\\1 \\2\\3",
  x = names(df_markov_all_vacc_w_incidence)
)

# Optional: make 'Cin1', 'Cin2', etc., uppercase as 'CIN1', 'CIN2' anywhere in the name
names(df_markov_all_vacc_w_incidence) <- gsub(
  pattern = "(\\b)Cin([123])",
  replacement = "\\1CIN\\2",
  x = names(df_markov_all_vacc_w_incidence),
  ignore.case = FALSE
)




################################################################################
# Preparing the Markov data for plotting

# List of strategies and associated row sources
strategies <- list(
  "0"  = df_markov_all_vacc_w_incidence %>% 
    dplyr::filter(sim.name == "Vaccination coverage: 0%"),
  
  #"60" = df_markov_strag_C,
  "60" = df_markov_all_vacc_w_incidence %>%
    dplyr::filter(sim.name == "Vaccination coverage: 60%"),
  
  #"70" = df_markov_strag_B,
  "70" = df_markov_all_vacc_w_incidence %>%
    dplyr::filter(sim.name == "Vaccination coverage: 70%"),
  
  #"80" = df_markov_strag_A
  "80" = df_markov_all_vacc_w_incidence %>%
    dplyr::filter(sim.name == "Vaccination coverage: 80%")
)

# Prefixes of interest (disease stages / categories)
#prefixes <- c("CIN1", "CIN2", "CIN3", "CC", "VPH")
prefixes <- c("n CIN1", "n CIN2","n CIN3", "n CC",
              "CIN1_Incidence", "CIN2_Incidence", "CIN3_Incidence", 
              "CC_Incidence", "HPVPrevalence", "CCMortality")

# Age brackets
age_from <- seq(10, 80, by = 5)
age_to   <- seq(14, 84, by = 5)
ages <- paste0(age_from, "-", age_to)

# Loop over strategies and prefixes
for (strategy in names(strategies)) {
  df <- strategies[[strategy]]
  
  for (prefix in prefixes) {
    # Construct full column names
    #cols <- paste0("n ", prefix, ages)
    cols <- paste(prefix, ages)
    
    # Extract values from first row (i.e., vaccination scenario)
    #values <- as.numeric(df[1, cols])
    #values <- df %>% dplyr::select(cols) %>% as.numeric()
    values <- df %>% dplyr::select(all_of(cols)) %>% as.numeric()
    
    # Add names to vector from original column names
    names(values) <- cols
    
    # New name format: Markov_CIN1_60
    var_name <- paste0("Markov_", prefix, "_", strategy)
    
    # Assign to global environment
    assign(var_name, values)
  }
}
################################################################################

# Create a list of all Markov_* objects in the global environment
markov_sim_vacc_incidences <- mget(ls(pattern = "^Markov_"))

# Save the list to an .RData file
#save(markov_sim_vacc, file = "data/markov_vacc_vectors.RData")
save(markov_sim_vacc, file = "data/markov_vacc_CORRECTED_vectors.RData")
#save(markov_sim_vacc_incidences, file = "data/markov_results/markov_vacc_incidences_vectors.RData")
save(markov_sim_vacc_incidences, file = "data/markov_results/markov_vacc_CORRECTED_incidences_vectors.RData")


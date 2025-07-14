## ADDING MARKOV RESULTS TO COMPARE AND VALIDATE MICROSIM:
# Load corresponding Markov results for validation:
load(file = "data/markov_results/results_only_cyto_Strategies_only_cyto_L_20250526.rda")
markov_cyto_result_df <- df; rm(df)

#markov_cyto_result_df %>% t() %>% View()

# Data Wrangling on markov_cyto_result_df:
library(dplyr)
library(stringr)

# Fix names
names(markov_cyto_result_df) <- names(markov_cyto_result_df) %>%
  # Capitalize Cyto, VPH, CIN, Cin, CC where appropriate
  str_replace_all("(?<=\\bn )Cyto", "CYTO") %>%
  str_replace_all("(?<=\\bn )VPH", "VPH") %>%
  str_replace_all("(?<=\\bn )CIN1", "CIN1") %>%
  str_replace_all("^Cin1_Incidence", "CIN1_Incidence") %>%
  str_replace_all("(?<=\\bn )CC", "CC") %>%
  # Add space before age groups like 10-14
  str_replace_all("(?<=[A-Za-z])(?=\\d{2}-\\d{2})", " ")

names(markov_cyto_result_df) <- names(markov_cyto_result_df) %>%
  # Capitalize 'Cin1', 'Cin2', 'Cin3' -> 'CIN1', etc.
  str_replace_all("\\bCin([123])", "CIN\\1") %>%
  # Ensure a space before the age range (e.g., 10-14, 15-19, ..., 80-84)
  str_replace_all("(?<! )(?=\\d{2}-\\d{2}\\b)", " ")


## Check the updated names
#names(markov_cyto_result_df)
markov_cyto_result_df <- markov_cyto_result_df %>%
  mutate(
    vaccination_coverage = 0,
    screening_coverage = 0.8,
    .before = 1  # Places the new columns at the beginning
  )

#markov_cyto_result_df %>% t() %>% View()


# For the moment we consider only the following Markov results in the loop:
for (strategy_name in names(sim_result[names(sim_result) != "runtime"])) {
  
  cat("Processing:", strategy_name, "\n")
  cat("For each strategy combination, corresponding Markov result is added")
  
  # Adding Markov incidences and prevalences:
  sim_result[[strategy_name]]$markov_CN1_incidences <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    #dplyr::select(matches("^n CIN1 "))
    dplyr::select(matches("^CIN1_Incidence "))
  
  sim_result[[strategy_name]]$markov_CN2_incidences <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    #dplyr::select(matches("^n CIN2 "))
    dplyr::select(matches("^CIN2_Incidence "))
  
  sim_result[[strategy_name]]$markov_CN3_incidences <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    #dplyr::select(matches("^n CIN3 "))
    dplyr::select(matches("^CIN3_Incidence "))
 
  sim_result[[strategy_name]]$markov_CC_incidences <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    dplyr::select(matches("^CC_Incidence "))
  
  sim_result[[strategy_name]]$markov_HPV_prevalences <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    dplyr::select(matches("^HPVPrevalence "))
  
  sim_result[[strategy_name]]$markov_CC_mortality <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    dplyr::select(matches("^CCMortality "))
  
  ## --
  
  sim_result[[strategy_name]]$markov_qaly_undisc <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    dplyr::select(`Per person QALYs und`) %>%
    as.numeric()
  
  sim_result[[strategy_name]]$markov_qaly_undis_Tot <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    dplyr::select(`Total QALYs und`) %>%
    as.numeric()
   
  sim_result[[strategy_name]]$markov_qaly_disc <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    dplyr::select(`Per person QALYs disc`) %>%
    as.numeric()
  
  sim_result[[strategy_name]]$markov_qaly_disc_Tot <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    dplyr::select(`Total QALYs disc`) %>%
    as.numeric()
  
  ## --
  
  sim_result[[strategy_name]]$markov_cost_undi <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    dplyr::select(`Per person D cost und`) %>%
    as.numeric()
  
  sim_result[[strategy_name]]$markov_cost_undis_Tot <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    dplyr::select(`Total D cost und`) %>%
    as.numeric()
  
  sim_result[[strategy_name]]$markov_cost_dis <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    dplyr::select(`Per person DM cost disc`) %>%
    as.numeric()
  
  sim_result[[strategy_name]]$markov_cost_disc_Tot <- 
    markov_cyto_result_df %>%
    dplyr::filter(sim.name == strategy_name) %>%
    dplyr::select(`Total D cost disc`) %>%
    as.numeric()
  
  
} # endfor strategy_name






















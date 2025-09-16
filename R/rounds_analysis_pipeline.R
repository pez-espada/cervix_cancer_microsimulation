#===============================================================================
# rounds_analysis_pipeline.R
# Script de análisis de "rounds" de cribado citológico
#===============================================================================

library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)

# ---------------- Helper ----------------
# Formatear screening ages para visualización
format_rounds <- function(ages_list) {
  ages <- unlist(ages_list)
  if(length(ages) <= 4) {
    paste(ages, collapse = ", ")
  } else if(length(ages) > 0) {
    paste0(min(ages), ", ..., ", max(ages))
  } else {
    NA_character_
  }
}

# ---------------- Pipeline ----------------
analyze_rounds <- function(sim_result) {
  
  # Ignorar elementos que no son estrategias
  strategy_names <- names(sim_result)[names(sim_result) != "runtime"]
  
  strategy_summaries <- map(strategy_names, function(sname) {
    
    rounds_df <- sim_result[[sname]]$rounds
    
    # Formato legible de screening ages
    rounds_df <- rounds_df %>%
      mutate(screening_fmt = sapply(screening_ages, format_rounds))
    
    # Resumen por batch / simulación
    batch_summary <- rounds_df %>%
      group_by(sim) %>%
      summarise(
        mean_rounds   = mean(rounds),
        sd_rounds     = sd(rounds),
        median_rounds = median(rounds),
        n_ind         = n(),
        prop_no_screen = mean(rounds == 0),
        .groups = "drop"
      )
    
    # Resumen agregado a nivel estrategia
    strategy_summary <- batch_summary %>%
      summarise(
        mean_rounds    = mean(mean_rounds),
        sd_rounds      = sd(mean_rounds),
        median_rounds  = median(median_rounds),
        prop_no_screen = mean(prop_no_screen)
      )
    
    # ---------------- Gráficos ----------------
    
    # 1. Distribución de rounds por estrategia
    p_rounds <- ggplot(rounds_df, aes(x = rounds)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
      labs(title = paste("Distribución de rounds -", sname),
           x = "Número de rounds",
           y = "Número de individuos")
    
    # 2. Distribución de edades de cribado (aplanar listas)
    age_df <- rounds_df %>%
      select(ID, screening_ages) %>%
      unnest_longer(screening_ages) %>%
      mutate(screening_ages = as.numeric(screening_ages))
    
    p_ages <- ggplot(age_df, aes(x = screening_ages)) +
      geom_histogram(binwidth = 1, fill = "darkgreen", color = "black") +
      labs(title = paste("Distribución de edades de cribado -", sname),
           x = "Edad cribado",
           y = "Número de individuos")
    
    list(
      rounds_individual = rounds_df,
      batch_summary     = batch_summary,
      strategy_summary  = strategy_summary,
      plot_rounds       = p_rounds,
      plot_ages         = p_ages
    )
    
  }) 
  
  names(strategy_summaries) <- strategy_names
  return(strategy_summaries)
}

# ---------------- Uso ----------------
# Guardar sim_result en un RData o RDS antes de correr este script
# load("sim_result.RData")
# sim_analysis <- analyze_rounds(sim_result)

# Luego puedes inspeccionar:
# sim_analysis$`25-34 cito 3 anys`$strategy_summary
# sim_analysis$`25-34 cito 3 anys`$plot_rounds
# sim_analysis$`25-34 cito 3 anys`$plot_ages


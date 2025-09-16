#===============================================================================
# rounds_analysis_pipeline_optimized.R
# Script de análisis de "rounds" optimizado para datasets grandes
#===============================================================================

library(data.table)
library(ggplot2)
library(furrr)
library(purrr)

# library(future)
# options(future.globals.maxSize = 5*1024^3)  # 5 GiB


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
analyze_rounds_optimized <- function(sim_result, sample_for_plot = 5000, n_cores = 4) {
  
  # Ignorar elementos que no son estrategias
  strategy_names <- names(sim_result)[names(sim_result) != "runtime"]
  
  plan(multisession, workers = n_cores)
  
  strategy_summaries <- future_map(strategy_names, function(sname) {
    
    rounds_df <- as.data.table(sim_result[[sname]]$rounds)
    
    # Formato legible de screening ages (muy rápido con data.table + map_chr)
    rounds_df[, screening_fmt := purrr::map_chr(screening_ages, format_rounds)]
    
    # Resumen por batch / simulación
    batch_summary <- rounds_df[, .(
      mean_rounds    = mean(rounds),
      sd_rounds      = sd(rounds),
      median_rounds  = median(rounds),
      n_ind          = .N,
      prop_no_screen = mean(rounds == 0)
    ), by = sim]
    
    # Resumen agregado a nivel estrategia
    strategy_summary <- batch_summary[, .(
      mean_rounds    = mean(mean_rounds),
      sd_rounds      = sd(mean_rounds),
      median_rounds  = median(median_rounds),
      prop_no_screen = mean(prop_no_screen)
    )]
    
    # ---------------- Gráficos ----------------
    # Muestreo para gráficos
    if(nrow(rounds_df) > sample_for_plot) {
      plot_sample <- rounds_df[sample(.N, sample_for_plot)]
    } else {
      plot_sample <- rounds_df
    }
    
    # 1. Distribución de rounds
    p_rounds <- ggplot(plot_sample, aes(x = rounds)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
      labs(title = paste("Distribución de rounds -", sname),
           x = "Número de rounds",
           y = "Número de individuos")
    
    # 2. Distribución de edades de cribado
    age_df <- plot_sample[, .(ID, screening_ages)]
    age_df <- age_df[, .(screening_ages = unlist(screening_ages)), by = ID]
    age_df[, screening_ages := as.numeric(screening_ages)]
    
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
    
  }, .options = furrr_options(seed = TRUE))
  
  names(strategy_summaries) <- strategy_names
  return(strategy_summaries)
}

# ---------------- Uso ----------------
# Guardar sim_result en un RData o RDS antes de correr este script
# load("sim_result.RData")
# sim_analysis <- analyze_rounds_optimized(sim_result, sample_for_plot = 5000, n_cores = 4)

# Inspeccionar resúmenes y gráficos
# sim_analysis$`25-34 cito 3 anys`$strategy_summary
# print(sim_analysis$`25-34 cito 3 anys`$plot_rounds)
# print(sim_analysis$`25-34 cito 3 anys`$plot_ages)

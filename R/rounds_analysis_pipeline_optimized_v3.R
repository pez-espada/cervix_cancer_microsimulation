#===============================================================================
# rounds_analysis_pipeline_optimized_v3.R
# Script de análisis de "rounds" optimizado para datasets grandes
# Incluye fix para as.data.table() en workers de future
#===============================================================================

library(data.table)
library(ggplot2)
library(furrr)
library(purrr)

# ---------------- Helper ----------------
# Formatear screening ages para visualización
format_rounds <- function(ages_list) {
  if (length(ages_list) == 0 || all(sapply(ages_list, is.null))) return(NA_character_)
  ages <- unlist(ages_list)
  if(length(ages) <= 4) {
    paste(ages, collapse = ", ")
  } else {
    paste0(min(ages), ", ..., ", max(ages))
  }
}

# ---------------- Pipeline ----------------
analyze_rounds_optimized_v3 <- function(sim_result, sample_for_plot = 5000, n_cores = 4) {
  
  strategy_names <- names(sim_result)[names(sim_result) != "runtime"]
  
  # Configurar paralelización
  plan(multisession, workers = n_cores)
  
  strategy_summaries <- future_map(strategy_names, function(sname) {
    
    # ---------------- Cargar paquetes dentro del worker ----------------
    data.table::setDTthreads(1) # evita conflictos de threads
    library(data.table)
    library(ggplot2)     # <--- esto asegura ggplot2 disponible en cada worker
    
    # Extraer solo el subset de rounds de esta estrategia (mucho más pequeño)
    rounds_df <- data.table::as.data.table(sim_result[[sname]]$rounds)
    
    # Aseguramos que screening_ages siempre tiene lista
    rounds_df[is.null(screening_ages), screening_ages := list(list())]
    
    # Formato legible
    rounds_df[, screening_fmt := purrr::map_chr(screening_ages, format_rounds)]
    
    # ---------------- Resumen por batch / simulación ----------------
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
    if(nrow(rounds_df) > sample_for_plot) {
      plot_sample <- rounds_df[sample(.N, sample_for_plot)]
    } else {
      plot_sample <- rounds_df
    }
    
    # Distribución de rounds
    p_rounds <- ggplot(plot_sample, aes(x = rounds)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
      labs(title = paste("Distribución de rounds -", sname),
           x = "Número de rounds",
           y = "Número de individuos")
    
    # Distribución de edades de cribado
    age_df <- plot_sample[, .(ID, screening_ages)]
    age_df <- age_df[, .(screening_ages = unlist(screening_ages)), by = ID]
    age_df[, screening_ages := as.numeric(screening_ages)]
    
    p_ages <- ggplot(age_df, aes(x = screening_ages)) +
      geom_histogram(binwidth = 1, fill = "darkgreen", color = "black") +
      labs(title = paste("Distribución de edades de cribado -", sname),
           x = "Edad cribado",
           y = "Número de individuos")
    
    # ---------------- Resultado por estrategia ----------------
    list(
      rounds_individual = rounds_df,
      batch_summary     = batch_summary,
      strategy_summary  = strategy_summary,
      plot_rounds       = p_rounds,
      plot_ages         = p_ages
    )
    
  }, .options = furrr_options(
    seed = TRUE, 
    #packages = "data.table",  # <--- asegura que data.table esté cargado en los workers
    packages = c("data.table", "ggplot2"),  # <--- ambos paquetes
    globals = c("format_rounds", "sample_for_plot")
  ))
  
  names(strategy_summaries) <- strategy_names
  return(strategy_summaries)
}

# ---------------- Uso ----------------
# load("sim_result.RData")
# sim_analysis <- analyze_rounds_optimized_v3(sim_result, sample_for_plot = 5000, n_cores = 4)

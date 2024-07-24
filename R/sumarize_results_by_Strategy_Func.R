summarize_results_by_Strategy <- function(results_list, numb_of_sims) {
  library(dplyr)
  
  result<-list()
  
  Strategy_name<-"No Intervention" #Sandra: en format funció ho passem com a variable incorporada als parametres de la estrategia.
  result[[Strategy_name]]<-list()
  #for(i in 1:numb_of_sims ){
  for(i in 1:length(results_list)){
    results_list[[i]]$tc_disc<-NULL
    results_list[[i]]$tc_undisc<-NULL
    results_list[[i]]$te_disc<-NULL
    results_list[[i]]$te_undisc<-NULL
    results_list[[i]][["Tot_Trans_per_t"]]<-NULL
  }
  
  sim<-results_list 
  names_sim<-names(sim[[1]])
  
  
  #for (name_level_of_sim in names_sim) {
  #  result[[Strategy_name]][[name_level_of_sim]] <- bind_rows(lapply(seq_along(sim), function(i) {
  #    df <- as.data.frame(sim[[i]][[name_level_of_sim]])  # Extraemos el data frame de cada lista
  #    df["sim"] <- i       # Agregamos una columna con el numero de simulacion
  #    
  #    ## Debugging
  #    #rownames(df) %>% print()  
  #    
  #    df["row_names"] <- as.numeric(rownames(df))
  #    #df <- df %>% mutate(cycle=row_number())
  #    rownames(df)<-NULL 
  #    df<-df[,c("sim","row_names",colnames(df)[1:(ncol(df)-2)])]
  #    return(df)
  #  }))
  
  for (name_level_of_sim in names_sim) {
    result[[Strategy_name]][[name_level_of_sim]] <- bind_rows(lapply(seq_along(sim), function(i) {
      df <- as.data.frame(sim[[i]][[name_level_of_sim]])  # Extraemos el data frame de cada lista
      df["sim"] <- i       # Agregamos una columna con el numero de simulacion
      if(name_level_of_sim!="TR"){
        df["row_names"] <- as.numeric(rownames(df))
        rownames(df)<-NULL 
        df<-df[,c("sim","row_names",colnames(df)[1:(ncol(df)-2)])]
      }else{
        rownames(df)<-NULL 
        df$cycle<-seq(1:nrow(df))
      }
      return(df)
    })) 
  
      
    #depuració
    if(sum(names(result[[Strategy_name]][[name_level_of_sim]])=="V1")>0 |
       sum( names(result[[Strategy_name]][[name_level_of_sim]])=="sim[[i]][[name_level_of_sim]]")>0) {
      colnames(result[[Strategy_name]][[name_level_of_sim]])[ ncol(result[[Strategy_name]][[name_level_of_sim]])]<-name_level_of_sim
    } 
    if(sum(endsWith( names(result[[Strategy_name]][[name_level_of_sim]]),"_disc"))>0 |sum(endsWith( names(result[[Strategy_name]][[name_level_of_sim]]),"_undisc"))>0  ){
      result[[Strategy_name]][[name_level_of_sim]][,"row_names"]<-NULL
      
    } 
    
  }
  
  new_result<-list()
  
  for (name_level_of_sim in names_sim) {
    new_result[[Strategy_name]][[name_level_of_sim]]<-as.data.frame( result[[Strategy_name]][[name_level_of_sim]])
    
    
  }
  
  return(new_result)
  
} # end of function

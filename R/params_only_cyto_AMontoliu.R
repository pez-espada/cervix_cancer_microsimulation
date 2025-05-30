#--------------------------------------------------
#Afegim edat final de cribat de HPV : Age_final_hpv
#--------------------------------------------------

provaCoverage<-function(grups_prova,Coverage){
  provaCoverage<-rep(0,15)
  for(i in grups_prova){
    
    provaCoverage[i]<-Coverage
  }
  
  
  return(provaCoverage)
  
}
 # edats de cribat<-




Parameters_strategy<-function(Coverage,cobertura_vacuna){
  
  ALL.PARAMETERS = list(
    
    ### Bloc 17
    
    list(sim.name="25-29 cito 3 anys",
       Screening = TRUE, ScreenType = 1,
         Age_first_cyto=25, Age_first_hpv=30, Age_final_hpv=30,
         screenPeriod = 3, 
         screenCoverage =provaCoverage(c(4),Coverage),
         vaccineCov=cobertura_vacuna
         ,additional_screening=0
         
    ),
    
    #### Bloc 18
    
    list(sim.name="25-34 cito 3 anys",
       Screening = TRUE, ScreenType = 1,
         Age_first_cyto=25, Age_first_hpv=35,Age_final_hpv=35,
         screenPeriod = 3,
         screenCoverage =provaCoverage(c(4,5),Coverage),
         vaccineCov=cobertura_vacuna
         ,additional_screening=0
    ),
    
    
    
    
    
    
    ### Bloc 19
    
    list(sim.name="30-34 cito 3 anys ",
       Screening = TRUE, ScreenType = 1,
         Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=35,
         screenPeriod = 3, 
         screenCoverage =provaCoverage(c(5),Coverage),
         vaccineCov=cobertura_vacuna
         ,additional_screening=0
         
    ),
    
    
    
    
    ### Bloc 20
    
    list(sim.name="25-29 cito 5 anys",
       Screening = TRUE, ScreenType = 1,
         Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=30,
         screenPeriod = 5, 
         screenCoverage =provaCoverage(c(4),Coverage),
         vaccineCov=cobertura_vacuna
         ,additional_screening=0
    ),
    
    
    ### Bloc 21
    
    list(sim.name="25-34 cito 5 anys",
       Screening = TRUE, ScreenType = 1,
         Age_first_cyto=25, Age_first_hpv=35, Age_final_hpv=35,
         screenPeriod = 5, 
         screenCoverage =provaCoverage(c(4,5),Coverage),
         vaccineCov=cobertura_vacuna
         ,additional_screening=0
         
    ),
    
    
    
    ### Bloc 22
    
    list(sim.name="30-34 cito 5 anys",
       Screening = TRUE, ScreenType = 1,
         Age_first_cyto=30, Age_first_hpv=35, Age_final_hpv=35,
         screenPeriod = 5, 
         screenCoverage =provaCoverage(c(5),Coverage),
         vaccineCov=cobertura_vacuna
         ,additional_screening=0
         
    )
  )
  
     return(ALL.PARAMETERS)
}

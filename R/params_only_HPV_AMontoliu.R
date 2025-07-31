

provaCoverage<-function(grups_prova,Coverage){
  provaCoverage<-rep(0,15)
  for(i in grups_prova){
    
    provaCoverage[i]<-Coverage
  }
  
  
  return(provaCoverage)
  
}
  




Parameters_strategy<-function(Coverage,cobertura_vacuna){
  
      
      ALL.PARAMETERS = list(
     
       
        ### Bloc 13
        
        list(sim.name="30-64 VPH 5 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=5,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
             
        ),
        list(sim.name="30-64 VPH 6 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=6,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="30-64 VPH 7 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=7,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="30-64 VPH 8 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=8,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="30-64 VPH 9 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=9,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="30-64 VPH 10 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=10,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        
        ### Bloc 14
        
        list(sim.name="35-64 VPH 5 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=35, Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=5,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="35-64 VPH 6 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=6,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="35-64 VPH 7 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=7,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="35-64 VPH 8 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=8,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="35-64 VPH 9 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=9,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="35-64 VPH 10 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=10,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        
        ### Bloc 15
        
        list(sim.name="30-69 VPH 5 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=5,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="30-69 VPH 6 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=6,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="30-69 VPH 7 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=7,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="30-69 VPH 8 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=8,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="30-69 VPH 9 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=30, Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=9,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="30-69 VPH 10 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=10,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        
        ### Bloc 16
        
        list(sim.name="35-69 VPH 5 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=5,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="35-69 VPH 6 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=6,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="35-69 VPH 7 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=7,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="35-69 VPH 8 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=8,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="35-69 VPH 9 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=9,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        ),
        list(sim.name="35-69 VPH 10 anys",
             Screening = TRUE, ScreenType = 1,
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=10,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             ,additional_screening=0
        )
      )
      return(ALL.PARAMETERS)
}
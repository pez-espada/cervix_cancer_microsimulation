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
  



Parameters_strategy<-function(Coverage,cobertura_vacuna){
  
      
      ALL.PARAMETERS = list(
     
        ### Bloc 1
        
        list(sim.name="25-29 cito 3 anys + 30-64 VPH 5 anys",
             Age_first_cyto=25, Age_first_hpv=30, Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=5, 
             screenCoverage =provaCoverage(c(4),Coverage),

             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
        
        ),
        list(sim.name="25-29 cito 3 anys + 30-64 VPH 6 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=6,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna 
        ),
        list(sim.name="25-29 cito 3 anys + 30-64 VPH 7 anys",
             Age_first_cyto=25, Age_first_hpv=30, Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=7,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
        ),
        list(sim.name="25-29 cito 3 anys + 30-64 VPH 8 anys",
             Age_first_cyto=25, Age_first_hpv=30, Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=8,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
        ),
        list(sim.name="25-29 cito 3 anys + 30-64 VPH 9 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=9,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
        ),
        list(sim.name="25-29 cito 3 anys + 30-64 VPH 10 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=10,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
        ),
        
        #### Bloc 2
        
        list(sim.name="25-34 cito 3 anys + 35-64 VPH 5 anys",
             Age_first_cyto=25, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=5,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
        ),
        list(sim.name="25-34 cito 3 anys + 35-64 VPH 6 anys",
             Age_first_cyto=25, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=6,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
        ),
        list(sim.name="25-34 cito 3 anys + 35-64 VPH 7 anys",
             Age_first_cyto=25, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=7,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
        ),
        list(sim.name="25-34 cito 3 anys + 35-64 VPH 8 anys",
             Age_first_cyto=25, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=8,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 3 anys + 35-64 VPH 9 anys",
             Age_first_cyto=25, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=9,
             screenCoverage =provaCoverage(c(4,5),Coverage),

             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 3 anys + 35-64 VPH 10 anys",
             Age_first_cyto=25, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=10,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),  
        
        ### Bloc 3
        
        list(sim.name="25-29 cito 3 anys + 30-69 VPH 5 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=5,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 3 anys + 30-69 VPH 6 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=6,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 3 anys + 30-69 VPH 7 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=7,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 3 anys + 30-69 VPH 8 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=8,
             screenCoverage =provaCoverage(c(4),Coverage), 
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 3 anys + 30-69 VPH 9 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=9,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 3 anys + 30-69 VPH 10 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=10,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        
        ### Bloc 4:
        
        list(sim.name="25-34 cito 3 anys + 35-69 VPH 5 anys",
             Age_first_cyto=25, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=5,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 3 anys + 35-69 VPH 6 anys",
             Age_first_cyto=25, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=6,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 3 anys + 35-69 VPH 7 anys",
             Age_first_cyto=25, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=7,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 3 anys + 35-69 VPH 8 anys",
             Age_first_cyto=25, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=8,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 3 anys + 35-69 VPH 9 anys",
             Age_first_cyto=25, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=9,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 3 anys + 35-69 VPH 10 anys",
             Age_first_cyto=25, Age_first_hpv=35, Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=10,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        
        ### Bloc 5
        
        list(sim.name="30-34 cito 3 anys + 35-64 VPH 5 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=5,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 3 anys + 35-64 VPH 6 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=6,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 3 anys + 35-64 VPH 7 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=7,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 3 anys + 35-64 VPH 8 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=8,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 3 anys + 35-64 VPH 9 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=9,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 3 anys + 35-64 VPH 10 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=10,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        
        ### Bloc 6
        
        list(sim.name="30-34 cito 3 anys + 35-69 VPH 5 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=5,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage( 6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 3 anys + 35-69 VPH 6 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=6,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage( 6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 3 anys + 35-69 VPH 7 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=7,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage( 6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 3 anys + 35-69 VPH 8 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=8,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage( 6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 3 anys + 35-69 VPH 9 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=9,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage( 6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 3 anys + 35-69 VPH 10 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=10,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage( 6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        
        ### Bloc 7
        
        list(sim.name="25-29 cito 5 anys + 30-64 VPH 5 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 5, HPVPeriod=5,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage( 5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 5 anys + 30-64 VPH 6 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 5, HPVPeriod=6,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage( 5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 5 anys + 30-64 VPH 7 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 5, HPVPeriod=7,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage( 5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 5 anys + 30-64 VPH 8 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 5, HPVPeriod=8,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage( 5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 5 anys + 30-64 VPH 9 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 5, HPVPeriod=9,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage( 5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 5 anys + 30-64 VPH 10 anys",
             Age_first_cyto=25, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 5, HPVPeriod=10,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage( 5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        
        ### Bloc 8
        
        list(sim.name="25-34 cito 5 anys + 35-64 VPH 5 anys",
             Age_first_cyto=25, Age_first_hpv=35, Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 5, HPVPeriod=5,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 5 anys + 35-64 VPH 6 anys",
             Age_first_cyto=25, Age_first_hpv=35, Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 5, HPVPeriod=6,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 5 anys + 35-64 VPH 7 anys",
             Age_first_cyto=25, Age_first_hpv=35, Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 5, HPVPeriod=7,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 5 anys + 35-64 VPH 8 anys",
             Age_first_cyto=25, Age_first_hpv=35, Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 5, HPVPeriod=8,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 5 anys + 35-64 VPH 9 anys",
             Age_first_cyto=25, Age_first_hpv=35, Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 5, HPVPeriod=9,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 5 anys + 35-64 VPH 10 anys",
             Age_first_cyto=25, Age_first_hpv=35, Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 5, HPVPeriod=10,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        
        ### Bloc 9
        
        list(sim.name="25-29 cito 5 anys + 30-69 VPH 5 anys",
             Age_first_cyto=25, Age_first_hpv=30, Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 5, HPVPeriod=5,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 5 anys + 30-69 VPH 6 anys",
             Age_first_cyto=25, Age_first_hpv=30,  Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 5, HPVPeriod=6,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 5 anys + 30-69 VPH 7 anys",
             Age_first_cyto=25, Age_first_hpv=30,  Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 5, HPVPeriod=7,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 5 anys + 30-69 VPH 8 anys",
             Age_first_cyto=25, Age_first_hpv=30,  Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 5, HPVPeriod=8,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 3 anys + 30-69 VPH 9 anys",
             Age_first_cyto=25, Age_first_hpv=30,  Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=9,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-29 cito 5 anys + 30-69 VPH 10 anys",
             Age_first_cyto=25, Age_first_hpv=30,  Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 5, HPVPeriod=10,
             screenCoverage =provaCoverage(c(4),Coverage),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        
        ### Bloc 10
        
        list(sim.name="25-34 cito 5 anys + 35-69 VPH 5 anys",
             Age_first_cyto=25, Age_first_hpv=35,  Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 5, HPVPeriod=5,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 5 anys + 35-69 VPH 6 anys",
             Age_first_cyto=25, Age_first_hpv=35, Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 5, HPVPeriod=6,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 5 anys + 35-69 VPH 7 anys",
             Age_first_cyto=25, Age_first_hpv=35, Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 5, HPVPeriod=7,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 5 anys + 35-69 VPH 8 anys",
             Age_first_cyto=25, Age_first_hpv=35, Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 5, HPVPeriod=8,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 5 anys + 35-69 VPH 9 anys",
             Age_first_cyto=25, Age_first_hpv=35, Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 5, HPVPeriod=9,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="25-34 cito 5 anys + 35-69 VPH 10 anys",
             Age_first_cyto=25, Age_first_hpv=35, Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 5, HPVPeriod=10,
             screenCoverage =provaCoverage(c(4,5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        
        ### Bloc 11
        
        list(sim.name="30-34 cito 5 anys + 35-64 VPH 5 anys",
             Age_first_cyto=30, Age_first_hpv=35, Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 5, HPVPeriod=5,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 5 anys + 35-64 VPH 6 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 5, HPVPeriod=6,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 5 anys + 35-64 VPH 7 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 5, HPVPeriod=7,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 5 anys + 35-64 VPH 8 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 5, HPVPeriod=8,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 5 anys + 35-64 VPH 9 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 5, HPVPeriod=9,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 5 anys + 35-64 VPH 10 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 5, HPVPeriod=10,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        
        ### Bloc 12
        
        list(sim.name="30-34 cito 5 anys + 35-69 VPH 5 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 5, HPVPeriod=5,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 5 anys + 35-69 VPH 6 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 5, HPVPeriod=6,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 5 anys + 35-69 VPH 7 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 5, HPVPeriod=7,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 5 anys + 35-69 VPH 8 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 5, HPVPeriod=8,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 5 anys + 35-69 VPH 9 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 5, HPVPeriod=9,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-34 cito 5 anys + 35-69 VPH 10 anys",
             Age_first_cyto=30, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 5, HPVPeriod=10,
             screenCoverage =provaCoverage(c(5),Coverage),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        
        ### Bloc 13
        
        list(sim.name="30-64 VPH 5 anys",
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=5,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),

             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-64 VPH 6 anys",
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=6,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-64 VPH 7 anys",
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=7,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-64 VPH 8 anys",
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=8,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),

             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-64 VPH 9 anys",
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=9,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-64 VPH 10 anys",
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=64,
             dnaScAgeGroups = 5:11, screenPeriod = 3, HPVPeriod=10,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        
        ### Bloc 14
        
        list(sim.name="35-64 VPH 5 anys",
             Age_first_cyto=0, Age_first_hpv=35, Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=5,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="35-64 VPH 6 anys",
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=6,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="35-64 VPH 7 anys",
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=7,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="35-64 VPH 8 anys",
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=8,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="35-64 VPH 9 anys",
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=9,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="35-64 VPH 10 anys",
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=64,
             dnaScAgeGroups = 6:11, screenPeriod = 3, HPVPeriod=10,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:11,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        
        ### Bloc 15
        
        list(sim.name="30-69 VPH 5 anys",
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=5,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-69 VPH 6 anys",
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=6,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-69 VPH 7 anys",
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=7,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-69 VPH 8 anys",
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=8,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-69 VPH 9 anys",
             Age_first_cyto=0, Age_first_hpv=30, Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=9,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="30-69 VPH 10 anys",
             Age_first_cyto=0, Age_first_hpv=30,Age_final_hpv=69,
             dnaScAgeGroups = 5:12, screenPeriod = 3, HPVPeriod=10,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(5:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        
        ### Bloc 16
        
        list(sim.name="35-69 VPH 5 anys",
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=5,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="35-69 VPH 6 anys",
             
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=6,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="35-69 VPH 7 anys",
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=7,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="35-69 VPH 8 anys",
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=8,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="35-69 VPH 9 anys",
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=9,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        ),
        list(sim.name="35-69 VPH 10 anys",
             Age_first_cyto=0, Age_first_hpv=35,Age_final_hpv=69,
             dnaScAgeGroups = 6:12, screenPeriod = 3, HPVPeriod=10,
             screenCoverage = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             
             dnaScreenCoverage=provaCoverage(6:12,Coverage) ,
             vaccineCov=cobertura_vacuna
             
        )
      )
      return(ALL.PARAMETERS)
}
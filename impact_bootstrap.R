sbs.fun<-function(y){
  sbs.taxon<-taxon_cube$data %>%
    filter(year==y) %>% 
    dplyr::select(scientificName,cellCode,obs) %>%
    group_by(scientificName,cellCode) %>%
    summarise(across(obs, sum), .groups = "drop") %>%
    pivot_wider(names_from = scientificName, values_from = obs) %>%
    arrange(cellCode) %>% 
    column_to_rownames(var = "cellCode")  #%>% 
    #mutate(across(all_of(everything()), ~ ifelse(is.na(.),0,.)))
    
  #colnames(sbs.taxon)<-NULL
  sbs.taxon<-as.matrix(sbs.taxon)
  return(sbs.taxon)
}


my_boot_statistic <- function(data, indices, fun) {
  d <- data[indices]
  y<-indices
  return(fun(d,y))
}

my_fun<-function(x,y){
  
  species_list<-colnames(x)

  if(!exists("taxon_status_list")){
    full_species_list<-sort(unique(taxon_cube$data$scientificName))
    taxon_status_list<-taxon_status(species_list = full_species_list,
                                    source = "WCVP",
                                    region = "South Africa")
  }
  


  intro.sf<-taxon_cube$data %>%
    filter(year==period[-c(1:12)][y]) %>%
    left_join(taxa_list_status,
              by = c("scientificName" = "taxon"))


  status.sf <- intro.sf %>%
    group_by(cellCode) %>%
    summarise(
      total_intro_obs = sum(obs[introduction_status == "introduced"], na.rm = TRUE),
      total_native_obs = sum(obs[introduction_status == "native"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(c(total_intro_obs, total_native_obs), ~ ifelse(.==0,NA,.))) %>%
    mutate(intro_native=total_intro_obs/total_native_obs) %>%
    arrange(cellCode)
  if (!exists("eicat_score_list")){
    full_species_list<-sort(unique(taxon_cube$data$scientificName))
    eicat_score_list=eicat_impact(eicat_data = eicat_data,species_list = full_species_list,
                                  fun="max")
  }
  
  eicat_score<-eicat_score_list[species_list,]
  
  siteScore<-status.sf$intro_native
  abdundance_impact = sweep(as.matrix(x),2,eicat_score,FUN = "*")
  impactScore = siteScore*abdundance_impact
  impact<-sum(impactScore,na.rm = TRUE) 
  
  
   return(indices)
}



#sbs.taxon_list<-map(period[-c(1:12)],sbs.fun)





fun=my_fun
samples<-10
bootstrap_list <- map(period[-c(1:12)],sbs.fun) %>% 
  # Perform bootstrapping
  purrr::map(~boot::boot(
    data = .,
    statistic = my_boot_statistic,
    R = samples,
    fun = fun))
bootstrap_list

sweep(as.matrix(sbs.taxon_list[[2]]), 2, NA, FUN = "*")


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
    
  colnames(sbs.taxon)<-NULL
  sbs.taxon<-as.matrix(sbs.taxon)
  return(sbs.taxon)
}


my_boot_statistic <- function(data, indices, fun) {
  d <- data[indices]
  return(fun(d))
}

my_fun<-function(x){
  
  m=mean(x,na.rm=TRUE)

  if (is.nan(m)) {
    m <- NA
  }
   return(m)
}

sbs.taxon_list<-map(period[-c(1:12)],sbs.fun)





fun=my_fun
samples<-1000
bootstrap_list <- map(period[-c(1:12)],sbs.fun) %>% 
  # Perform bootstrapping
  purrr::map(~boot::boot(
    data = .,
    statistic = my_boot_statistic,
    R = samples,
    fun = fun))
bootstrap_list

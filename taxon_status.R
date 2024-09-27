taxon_status<-function(species_list,
                       source="manual",
                       status_data=NULL,
                       region=NULL,
                       col_scientificName=NULL,
                       col_introductionStatus=NULL){
  
  if(source=="WCVP"){
    if(is.null(region)){
      stop("region needs to be specified for WCVP")
    } else 
    {
      # Download WCVP for native taxa in area of interest
      native_list <- rWCVP::wcvp_checklist() %>%
        filter(area_code_l3 %in% rWCVP::get_wgsrpd3_codes(region)) %>%
        filter((accepted_name %in% species_list) & occurrence_type=="native")
      
      # create taxa list
      species_list<-data.frame(taxon=species_list)
      
      # create new dataframe with introduction status
      taxa_list_status<-species_list%>%
        mutate(introduction_status = ifelse(taxon%in%native_list$accepted_name,
                                            "native","introduced"))
    }
    
  } else if (source=="manual"){
    if(is.null(status_data)){
      stop("status_data needs to be provided for `manual` source")
    } else if ("data.frame" %in% class(status_data)){
      
      if(all(c("scientificName","introductionStatus")%in%names(status_data))){
        data<-status_data %>% 
          filter((scientificName %in% species_list) & introductionStatus=="native")
        # create species list
        species_list<-data.frame(taxon=species_list)
        
        # create new dataframe with introduction status
        taxa_list_status<-species_list%>%
          mutate(introduction_status = ifelse(taxon%in%data$scientificName,
                                              "native","introduced"))
        
      }
      else if(!(is.null(col_scientificName)|is.null(col_introductionStatus))){
        data<-status_data %>% 
          rename(c(scientificName=col_scientificName,
                   introductionStatus=col_introductionStatus)) %>% 
          filter((scientificName %in% species_list) & introductionStatus=="native")
        # create species list
        species_list<-data.frame(taxon=species_list)
        
        # create new dataframe with introduction status
        taxa_list_status<-species_list%>%
          mutate(introduction_status = ifelse(taxon%in%data$scientificName,
                                              "native","introduced"))
        
      }
      else {stop("both `col_scientificName` and `col_introductionStatus` not found in `status_data` nor specified ")}
    } else{ stop("status_data must be a dataframe")}
  } else {stop("`source` should either be WCVP or manual")}
  
  
  return(taxa_list_status)
}

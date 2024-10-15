sbtFun<-function(tryfile,taxa_cube){
  # read the try data if path is given
  if("character" %in% class(tryfile)){
    trydata<-rtry::rtry_import(
      input=tryfile,
      separator = "\t",
      encoding = "Latin-1",
      quote = "",
      showOverview = TRUE
    )
  } else if("data.frame" %in% class(tryfile)){ #Check if tryfile contains necessary columns
    if(any(!c("AccSpeciesName","TraitID","TraitName","OrigValueStr") %in% colnames(tryfile))){
      requiredcol<-c("AccSpeciesName","TraitID","TraitName","OrigValueStr")
      missingcol<-requiredcol[!c("AccSpeciesName","TraitID","TraitName","OrigValueStr") %in% colnames(tryfile)]
      cli::cli_abort(c("{missingcol} is/are not in the {.var tryfile} column "))
    }
    
    trydata<-tryfile
  } else { # stop and report if tryfile is not a file path or dataframe
    cli::cli_abort(c("{.var tryfile} is not a file path or dataframe"))
  }
  
  
  # extract unique species name from GBIF occurrence data
  species_list <- sort(unique(taxa_cube$taxa$data$scientificName))
  
  SpeciesbyTrait<-trydata %>%
    # drop rows which contains no trait
    tidyr::drop_na(TraitID) %>%
    # select species name, trait and trait value
    dplyr::select(AccSpeciesName,TraitID,OrigValueStr) %>%
    #group by Species and trait
    dplyr::group_by(AccSpeciesName,TraitID) %>%
    #choose the first trait value if there are multiples trait for a species
    dplyr::summarise(across(OrigValueStr, first), .groups = "drop") %>%
    # reshape to wide format to have specie by trait dataframe
    tidyr::pivot_wider(names_from = TraitID, values_from = OrigValueStr) %>%
    # select species that are only present in gbif data
    dplyr::filter(AccSpeciesName %in% species_list) %>%
    # convert species names to row names
    tibble::column_to_rownames(var = "AccSpeciesName") %>% 
    # select traits that contain at least a single value
    dplyr::select_if(~ !all(is.na(.)))
  # add the rows of the remaining species without traits from TRY
  na.df<-as.data.frame(matrix(NA,nrow = length(setdiff(species_list,rownames(SpeciesbyTrait))),
                              ncol = ncol(SpeciesbyTrait)))
  row.names(na.df)<-setdiff(species_list,rownames(SpeciesbyTrait))
  names(na.df)<-names(SpeciesbyTrait) 
  SpeciesbyTrait<-rbind(SpeciesbyTrait,na.df)
  #collect traitID
  trait<-colnames(SpeciesbyTrait)
  #sort rows according to unique species list
  SpeciesbyTrait<-as.data.frame(SpeciesbyTrait[species_list,])
  
  # Download WCVP for native taxa in area of interest
  native_list <- rWCVP::wcvp_checklist() %>%
    filter(area_code_l3 %in% rWCVP::get_wgsrpd3_codes("South Africa")) %>%
    filter((accepted_name %in% species_list) & occurrence_type=="native")
  
  # create taxa list
  taxa_list<-data.frame(taxon=species_list)
  
  # create new dataframe with introduction status
  taxa_list_status<-taxa_list%>%
    mutate(introduction_status = ifelse(taxon%in%native_list$accepted_name,
                                        "native","introduced"))
  # add introduction status to trait column
  SpeciesbyTrait$introduction_status<-taxa_list_status$introduction_status
  
  # create species by trait matrix
  sbtM<-as.matrix(SpeciesbyTrait)
  
  #remove column and row names 
  rownames(sbtM)<-NULL
  colnames(sbtM)<-NULL
  
  
  #extract the trait names
  traitname<-trydata %>%
    # drop rows which contains no trait
    tidyr::drop_na(TraitID) %>%
    # select the trait ID and trait name
    dplyr::select(TraitID,TraitName) %>%
    dplyr::group_by(TraitID) %>%
    dplyr::summarise(across(TraitName, first), .groups = "drop") %>%
    tibble::column_to_rownames("TraitID")
  #create trait name to align with the sbt column
  traitname<-c(traitname[trait,],"Introduction status")
  traitname<-data.frame('TraitID'=c(trait,"Introduction status"),'TraitName'=traitname)
  return(list("sbt"=sbtM,"traitname"=traitname))
  
}

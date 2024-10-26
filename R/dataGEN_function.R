# >> dataGEN function
# ------------------------------------------------------------------------------


#' @title Prepare data for invasib()
#' 
#' @description
#' Prepare site-by-species (sbs), species-by-trait (sbt) and 
#' site-by-environment (sbe) matrices to calculate invasibility cube. 
#' The function `dataGEN` computes and returns appendix data if specified
#' by the user.
#' 
#' @param taxa Character or dataframe. The character should be the scientific 
#' name of the focal taxa while the dataframe is the GBIF occurrences data which must 
#' contain "decimalLatitude","decimalLongitude","species","speciesKey",
#' "coordinateUncertaintyInMeters","dateIdentified", and "year".
#' @param region.sf sf object. The shapefile of the region of study
#' @param option Numeric. 1 to compute and return sbs, sbt and sbe. 2 to compute
#' and return sbs and sbe.
#' @param limit Integer. Number of records to return from GBIF download.
#' Default is 500
#' @param ref Character or dataframe.The character should be the scientific name 
#' of the reference taxa while the dataframe is the GBIF occurrences data which 
#' must contain "decimalLatitude","decimalLongitude","species","speciesKey",
#' "coordinateUncertaintyInMeters","dateIdentified", and "year".
#' @param country Character. Country for which the GBIP occurrences data should 
#' be downloaded. Country should be 2-letter country code (ISO-3166-1).
#' Default is ZA
#' @param res Numeric. The resolution of grid cells to be used. Default is 0.25
#' @param tryfile Character. File path of the *.txt of the TRY 'released' data or 
#' dataframe imported by rtry::rtry_import()
#' #' @param rastfile Raster. should be path for storing the worldclim data 
#' (see geodata::worldclim_global) or any other raster data of environmental 
#' variable
#'
#' @return A list containing the containing sbs, sbt, sbe matrices and their 
#' appendices if specified.
#'
#'#' @examples
#' \dontrun{
#' KZN_sf<-readRDS(paste0(getwd(),"/KZN_sf.rds"))
#' taxa_Fabaceae_KZN <- readRDS(paste0(getwd(),"/Fabaceae_KZN.rds"))
#' species_trait<-readRDS(paste0(getwd(),"/TRY_traits.rds"))
#' rast_path <- "C:/Users/mukht/Documents"
#' datalist <- dataGEN(taxa=taxa_Fabaceae_KZN,
#'                          region.sf=KZN_sf,
#'                          option=1,
#'                          limit=500,
#'                          ref=NULL,
#'                          res=0.25,
#'                          tryfile=species_trait,
#'                          country='ZA',
#'                          rastfile=rast_path,
#'                          appendix=TRUE)
#' }
#' 

# # install.packages("devtools")
# #devtools::install_github("b-cubed-eu/b3gbi")

# library(b3gbi)
# library(readr)
# library(tidyverse)
# library(rgbif) 
# library(terra)
# library(sf)
# library(rtry) 
# library(tidyverse)
# library(geodata)


# KZN_sf<-readRDS(paste0(getwd(),"/KZN_sf.rds"))
# taxa_Fabaceae_KZN <- readRDS(paste0(getwd(),"/Fabaceae_KZN.rds"))
# species_trait<-readRDS(paste0(getwd(),"/TRY_traits.rds"))
# rast_path <- "C:/Users/mukht/Documents" 

# datalist <- dataGEN(taxa=taxa_Fabaceae_KZN,
#                     region.sf=KZN_sf,
#                     option=1,
#                     limit=500,
#                     ref=NULL,
#                     res=0.25,
#                     tryfile=species_trait,
#                     country='ZA',
#                     rastfile=rast_path,
#                     appendix=TRUE)



dataGEN = function(taxa,
                   region.sf,
                   option,
                   limit=500,
                   ref=NULL,
                   res=0.25,
                   tryfile,
                   country='ZA',
                   stateProvince=NULL,
                   rastfile,
                   appendix=FALSE){
  

  
  if(option==1){
    
    taxa_cube <- taxaFun(taxa=taxa,
                         region.sf=region.sf,
                         limit=limit, 
                         ref=ref,
                         country=country,
                         res=res,
                         stateProvince=stateProvince)
    
    sbs <- sbsFun(taxa_cube = taxa_cube,
                  region.sf = region.sf, 
                  res = res,
                  appendix = appendix)
    
    sbt <- sbtFun(tryfile = tryfile,
                  taxa_cube = taxa_cube,
                  appendix=appendix)
    
    siteID <- sbs$siteID
    sbe <- sbeFun(rastfile=rastfile,
                  region.sf=region.sf,
                  res=res,siteID=siteID,
                  appendix=appendix)
    
    return(list("sbs"=sbs,"sbt"=sbt,"sbe"=sbe))
  } else if(option==2){
    
    taxa_cube <- taxaFun(taxa=taxa,
                         region.sf=region.sf,
                         limit=limit, 
                         ref=ref,
                         country=country,
                         res=res,
                         stateProvince=stateProvince)
    
    sbs <- sbsFun(taxa_cube = taxa_cube,
                  region.sf = region.sf, 
                  res = res,
                  appendix = appendix)
    
    siteID <- sbs$siteID
    sbe <- sbeFun(rastfile=rastfile,
                  region.sf=region.sf,
                  res=res,
                  siteID=siteID)
    
    return(list("sbs"=sbs,"sbe"=sbe))
  } else{
    message("option has to be 1 (returns sbs,sbt,sbe) or 2 (returns sbs and sbe)")}
  
  
}

# >> taxaFun function
# ------------------------------------------------------------------------------



#' @title Prepare Data Cubes
#'
#'@description Prepare data cube to calculate Species-by-site (sbs), 
#' Species-by-trait (sbt), and Site-by-environment (sbe). The function `taxaFun`
#' can take in the scientific name of the taxa of interest as in character or 
#' GBIF occurrences data containing necessary columns. The GBIF occurrences is 
#' downloaded if scientific names is given. The function returns data cubes of 
#' focal taxa and reference taxa. If the reference taxa is not specified, the 
#' function returns the focal taxa as reference taxa.
#'@param taxa Character or dataframe. The character should be the scientific 
#'name of the focal taxa while the dataframe is the GBIF occurrences data which must 
#'contain "decimalLatitude","decimalLongitude","species","speciesKey",
#'"coordinateUncertaintyInMeters","dateIdentified", and "year".
#'@param region.sf sf object. The shapefile of the region of study
#'@param limit Integer. Number of records to return from GBIF download.
#'Default is set to 500
#'@param ref Character or dataframe.The character should be the scientific name 
#'of the reference taxa while the dataframe is the GBIF occurrences data which 
#'must contain "decimalLatitude","decimalLongitude","species","speciesKey",
#'"coordinateUncertaintyInMeters","dateIdentified", and "year".
#' @param country Character. Country for which the GBIP occurrences data should 
#' be downloaded. Country should be 2-letter country code (ISO-3166-1).
#' Default is ZA
#' @param res Numeric. The resolution of grid cells to be used. Default is 0.25
#'
#' @return A list containing the `sim_cubes` of focal and reference taxa.
#'
#'@import dplyr
#'
#' @examples
#' \dontrun{
#' KZN_sf<-readRDS(paste0(getwd(),"/KZN_sf.rds"))
#' taxa_Fabaceae_KZN <- readRDS(paste0(getwd(),"/Fabaceae_KZN.rds"))
#' taxa_cube <- taxaFun(taxa=taxa_Fabaceae_KZN, region.sf=KZN_sf)
#' }
#' 





taxaFun <- function(taxa,
                    region.sf,
                    limit=500, 
                    ref=NULL,
                    country=NULL,
                    res=0.25,
                    stateProvince=NULL){
  
  grid <- region.sf %>%
    sf::st_make_grid(cellsize = c(res,res),
                     offset = c(sf::st_bbox(region.sf)$xmin,
                                sf::st_bbox(region.sf)$ymin)) %>%
    sf::st_sf() %>%
    dplyr::mutate(cellid = dplyr::row_number())
  
  # download taxa (target taxa) if the scientific name is given as character
  if("character" %in% class(taxa)){
    taxa.gbif_download = rgbif::occ_data(scientificName=taxa, # download data from gbif
                                         country=country,
                                         stateProvince = stateProvince,
                                         hasCoordinate=TRUE,
                                         hasGeospatialIssue=FALSE,
                                         limit = limit)
    
    taxa.df = as.data.frame(taxa.gbif_download$data) #extract data from the downloaded file
  } else if("data.frame" %in% class(taxa)){ #check if data fame contains the required columns
    if(any(!c("decimalLatitude","decimalLongitude",
              "species","speciesKey","coordinateUncertaintyInMeters","dateIdentified","year") %in% colnames(taxa))){
      requiredcol<-c("decimalLatitude","decimalLongitude","species","speciesKey","coordinateUncertaintyInMeters","dateIdentified","year")
      missingcol<-requiredcol[!c("decimalLatitude","decimalLongitude","species","speciesKey","coordinateUncertaintyInMeters","dateIdentified","year") %in% colnames(taxa)]
      cli::cli_abort(c("{missingcol} is/are not in the {.var taxa} column ", 
                       "x" = "{.var taxa} should be a data of GBIF format "))
    }
    # take taxa data frame if accurate
    taxa.df<-taxa
  } else { # stop and report if taxa is not a scientific name or dataframe
    cli::cli_abort(c("{.var taxa} is not a character or dataframe"))
  }
  
  
  taxa.sf = taxa.df %>%
    dplyr::select(decimalLatitude,decimalLongitude,
                  species,speciesKey,coordinateUncertaintyInMeters,dateIdentified,year) %>% #select occurrence data
    dplyr::filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
    dplyr::filter(coordinateUncertaintyInMeters<=res*1000) %>% 
    #dplyr::mutate(coordinateUncertaintyInMeters = coordinateUncertaintyInMeters/(res*1000)^2) %>% 
    dplyr::mutate(dateIdentified = as.Date(dateIdentified)) %>%  # convert date to date format
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
                 crs = 4326) %>% 
    sf::st_join(grid) %>% 
    as.data.frame() %>% 
    dplyr::select(-geometry) %>% 
    dplyr::mutate(occurrences=1) %>% 
    dplyr::group_by(species,coordinateUncertaintyInMeters,year,speciesKey,
                    #iucnRedListCategory,
                    cellid) %>% 
    dplyr::summarise(across(occurrences, sum), .groups = "drop")
  
  
  taxa_cube<-b3gbi::process_cube(taxa.sf,grid_type = "custom",
                                 cols_cellCode = "cellid", cols_year = "year",
                                 cols_occurrences = "occurrences",
                                 cols_species = "species",cols_speciesKey = "speciesKey",
                                 cols_minCoordinateUncertaintyInMeters = "coordinateUncertaintyInMeters")
  # download reference taxa if provided or use the target taxa if otherwise
  if(!is.null(ref)){
    if("character" %in% class(ref)){
      ref.gbif_download = rgbif::occ_data(scientificName=ref, # download data from gbif
                                          country=country,
                                          stateProvince = stateProvince,
                                          hasCoordinate=TRUE,
                                          hasGeospatialIssue=FALSE,
                                          limit = limit)
      
      ref.df = as.data.frame(ref.gbif_download$data)
      
      ref.sf = ref.df %>%
        dplyr::select(decimalLatitude,decimalLongitude,
                      species,speciesKey,coordinateUncertaintyInMeters,dateIdentified,year) %>% #select occurrence data
        dplyr::filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
        dplyr::filter(coordinateUncertaintyInMeters<=res*1000) %>% 
        dplyr::mutate(dateIdentified = as.Date(dateIdentified)) %>% # convert date to date format
        sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
                     crs = 4326) %>% 
        sf::st_join(grid) %>% 
        as.data.frame() %>% 
        dplyr::select(-geometry) %>% 
        dplyr::mutate(occurrences=1) %>% 
        dplyr::group_by(species,coordinateUncertaintyInMeters,year,speciesKey,
                        #iucnRedListCategory,
                        cellid) %>% 
        dplyr::summarise(across(occurrences, sum), .groups = "drop")
      
      
      ref_cube<-b3gbi::process_cube(ref.sf,grid_type = "custom",
                                    cols_cellCode = "cellid", cols_year = "year",
                                    cols_occurrences = "occurrences",
                                    cols_species = "species",cols_speciesKey = "speciesKey",
                                    cols_minCoordinateUncertaintyInMeters = "coordinateUncertaintyInMeters")
    } else if("data.frame" %in% class(ref)){ # check columns of ref if provided
      if(any(!c("decimalLatitude","decimalLongitude",
                "species","speciesKey","coordinateUncertaintyInMeters","dateIdentified","year") %in% colnames(ref))){
        requiredcol<-c("decimalLatitude","decimalLongitude","species","speciesKey","coordinateUncertaintyInMeters","dateIdentified","year")
        missingcol<-requiredcol[!c("decimalLatitude","decimalLongitude","species","speciesKey","coordinateUncertaintyInMeters","dateIdentified","year") %in% colnames(ref)]
        cli::cli_abort(c("{missingcol} is/are not in the {.var ref} column ", 
                         "x" = "{.var ref} should be a data of GBIF format "))
      }
      ref.sf = ref %>%
        dplyr::select(decimalLatitude,decimalLongitude,
                      species,speciesKey,coordinateUncertaintyInMeters,dateIdentified,year) %>% #select occurrence data
        dplyr::filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
        dplyr::filter(coordinateUncertaintyInMeters<=res*1000) %>%
        dplyr::mutate(dateIdentified = as.Date(dateIdentified)) %>% # convert date to date format
        sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
                     crs = 4326) %>% 
        sf::st_join(grid) %>% 
        as.data.frame() %>% 
        dplyr::select(-geometry) %>% 
        dplyr::mutate(occurrences=1) %>% 
        dplyr::group_by(species,coordinateUncertaintyInMeters,year,speciesKey,
                        #iucnRedListCategory,
                        cellid) %>% 
        dplyr::summarise(across(occurrences, sum), .groups = "drop")
      
      
      ref_cube<-b3gbi::process_cube(ref.sf,grid_type = "custom",
                                    cols_cellCode = "cellid", cols_year = "year",
                                    cols_occurrences = "occurrences",
                                    cols_species = "species",cols_speciesKey = "speciesKey",
                                    cols_minCoordinateUncertaintyInMeters = "coordinateUncertaintyInMeters")
      
    } else { # stop and report if taxa is not a scientific name or dataframe
      cli::cli_abort(c("{.var ref} is not a character or dataframe"))
    }
    
  } else {
    ref_cube <- taxa_cube
  }
  
  return(list("taxa"=taxa_cube,"ref"=ref_cube))
}

# >> sbsFUN function
# ------------------------------------------------------------------------------

#' @title Create species-by-site matrix
#' 
#' @description Create species-by-site (sbs) matrix over the given region and
#' time. The rows and columns of the matrix represents the sites and species
#' respectively. 
#' @param taxa_cube List. A list containing `sim_cubes`of the focal and reference
#' taxa.
#' @param region.sf sf object. The shapefile of the region of study
#' @param res Numeric. The resolution of grid cells to be used. Default is 0.25
#' @param col_temporal Character. The name of the column containing the temporal
#' dimension, e.g., "year". Default is NULL
#'
#' 
#' #' @examples
#' \dontrun{
#' sbs<-sbsFun(taxa_cube=taxa_cube,
#'              region.sf=KZN_sf,
#'              appendix=TRUE)
#' }


sbsFun <- function(taxa_cube,
                   region.sf,
                   res=0.25,
                   col_temporal=NULL,
                   appendix=FALSE){
  
  
  if(appendix){
    # extract unique species name from GBIF occurrence data
    species_list<-sort(unique(taxa_cube$taxa$data$scientificName))
    
    sbsM<-taxa_cube$taxa$data %>%
      dplyr::select(scientificName,cellCode,obs) %>%
      dplyr::group_by(scientificName,cellCode) %>%
      dplyr::summarise(across(obs, sum), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
      dplyr::arrange(cellCode) %>% 
      tibble::column_to_rownames(var = "cellCode") %>% 
      dplyr::mutate_all(~ replace(., is.na(.), 0))
    
    #create species-by-site-by-time (sbsbt) if col_temporal is provided
    sbsbtM<-NULL
    if(!is.null(col_temporal)){
      sbsbtM<-list()
      period<-unlist(unique(taxa_cube$taxa$data[col_temporal]))
      for(t in period){
        sbsbtM[[as.character(t)]]<-taxa_cube$taxa$data %>%
          filter(.data[[col_temporal]]==t) %>% 
          dplyr::select(scientificName,cellCode,obs) %>%
          dplyr::group_by(scientificName,cellCode) %>%
          dplyr::summarise(across(obs, sum), .groups = "drop") %>%
          tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
          dplyr::arrange(cellCode) %>% 
          tibble::column_to_rownames(var = "cellCode") 
      }
      
    }
    
    #create grid for region
    grid <- region.sf %>%
      sf::st_make_grid(cellsize = c(res,res),
                       offset = c(sf::st_bbox(region.sf)$xmin,
                                  sf::st_bbox(region.sf)$ymin)) %>%
      sf::st_sf() %>%
      dplyr::mutate(cellid = dplyr::row_number())
    
    # get coordinates of the occurrence sites
    coords <- sf::st_coordinates(sf::st_centroid(grid))
    coords <- coords[as.integer(rownames(sbsM)),]
    colnames(coords)<-c("x","y")
    
    
    # create binary matrix
    sbsM.binary<-sbsM %>% 
      dplyr::mutate(across(everything(), ~ ifelse(. >= 1, 1, .)))
    
    
    #create site uncertainty
    site_uncertainty <- taxa_cube$taxa$data %>% 
      dplyr::select(all_of(c("cellCode","minCoordinateUncertaintyInMeters"))) %>% 
      dplyr::group_by(cellCode) %>% 
      dplyr::summarise(across(minCoordinateUncertaintyInMeters, sum), .groups = "drop")
    
    
    # compute sbs for ref if it is different from taxa
    if(identical(taxa_cube$taxa,taxa_cube$ref)){
      sbsM.ref<-sbsM
    } else{
      
      sbsM.ref<-taxa_cube$ref$data %>%
        dplyr::select(scientificName,cellCode,obs) %>%
        dplyr::group_by(scientificName,cellCode) %>%
        dplyr::summarise(across(obs, sum), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
        dplyr::arrange(cellCode) %>% 
        tibble::column_to_rownames(var = "cellCode") %>% 
        dplyr::mutate_all(~ replace(., is.na(.), 0))
      
    }
    
    return(list("sbs"=sbsM,"sbs.ref"=sbsM.ref,"sbs.binary"=sbsM.binary,
                "coords"=coords,"species_list"=species_list,
                "siteID"=rownames(sbsM),"site_unc"=site_uncertainty,
                "sbsbt"=sbsbtM))
  } else{
    
    sbsM<-taxa_cube$taxa$data %>%
      dplyr::select(scientificName,cellCode,obs) %>%
      dplyr::group_by(scientificName,cellCode) %>%
      dplyr::summarise(across(obs, sum), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
      dplyr::arrange(cellCode) %>% 
      tibble::column_to_rownames(var = "cellCode") %>% 
      dplyr::mutate_all(~ replace(., is.na(.), 0))
    
    #create grid for region
    grid <- region.sf %>%
      sf::st_make_grid(cellsize = c(res,res),
                       offset = c(sf::st_bbox(region.sf)$xmin,
                                  sf::st_bbox(region.sf)$ymin)) %>%
      sf::st_sf() %>%
      dplyr::mutate(cellid = dplyr::row_number())
    
    # get coordinates of the occurrence sites
    coords <- sf::st_coordinates(sf::st_centroid(grid))
    coords <- coords[as.integer(rownames(sbsM)),]
    colnames(coords)<-c("x","y")
    
    #create species-by-site-by-time (sbsbt) if col_temporal is provided
    sbsbtM<-NULL
    if(!is.null(col_temporal)){
      sbsbtM<-list()
      period<-unlist(unique(taxa_cube$taxa$data[col_temporal]))
      for(t in period){
        sbsbtM[[as.character(t)]]<-taxa_cube$taxa$data %>%
          filter(.data[[col_temporal]]==t) %>% 
          dplyr::select(scientificName,cellCode,obs) %>%
          dplyr::group_by(scientificName,cellCode) %>%
          dplyr::summarise(across(obs, sum), .groups = "drop") %>%
          tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
          dplyr::arrange(cellCode) %>% 
          tibble::column_to_rownames(var = "cellCode")
      }
      
    }
    
    return(list("sbs"=sbsM,
                "coords"=coords,
                "siteID"=rownames(sbsM)))
  }
  
}

# >> sbtFUN function
# ------------------------------------------------------------------------------

#' @title Create species-by-trait matrix
#' 
#' @description Create species-by-trait (sbt) matrix for the given species.
#' The rows and columns of the matrix represents the species and trait
#' variables respectively. 
#' @param tryfile Character. File path of the *.txt of the TRY 'released' data or 
#' dataframe imported by rtry::rtry_import()
#' @param taxa_cube List. A list containing `sim_cubes`of the focal and reference
#' taxa.
#' #' @examples
#' \dontrun{
#' species_trait<-readRDS(paste0(getwd(),"/TRY_traits.rds"))
#' taxa_cube <- taxaFun(taxa=taxa_Fabacae, region.sf=KZN_sf)
#' sbt <- sbtFun(tryfile=species_trait, taxa_cube=taxa_cube)
#' }
sbtFun<-function(tryfile,
                 taxa_cube,
                 appendix=FALSE){
  
  
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
  SpeciesbyTrait<-rbind(SpeciesbyTrait,na.df) %>% 
    dplyr::mutate_if(is.character, as.numeric) %>% 
    dplyr::select_if(~ !all(is.na(.))) %>% 
    dplyr::mutate_all(~ replace(., is.na(.), 0))
  
  
  
  #collect traitID
  trait<-colnames(SpeciesbyTrait)
  #sort rows according to unique species list
  SpeciesbyTrait<-as.data.frame(SpeciesbyTrait[species_list,])
  
  # # Download WCVP for native taxa in area of interest
  # native_list <- rWCVP::wcvp_checklist() %>%
  #   filter(area_code_l3 %in% rWCVP::get_wgsrpd3_codes("South Africa")) %>%
  #   filter((accepted_name %in% species_list) & occurrence_type=="native")
  # 
  # # create taxa list
  # taxa_list<-data.frame(taxon=species_list)
  # 
  # # create new dataframe with introduction status
  # taxa_list_status<-taxa_list%>%
  #   mutate(introduction_status = ifelse(taxon%in%native_list$accepted_name,
  #                                       "native","introduced"))
  # # add introduction status to trait column
  # SpeciesbyTrait$introduction_status<-taxa_list_status$introduction_status
  # 
  # # create species by trait matrix
  # sbtM<-as.matrix(SpeciesbyTrait)
  
  # #remove column and row names 
  # rownames(sbtM)<-NULL
  # colnames(sbtM)<-NULL
  
  
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
  traitname<-c(traitname[trait,])
  traitname<-data.frame('TraitID'=paste0("trait_",1:ncol(SpeciesbyTrait)),
                        'TraitName'=traitname)
  
  names(SpeciesbyTrait)<-paste0("trait_",1:ncol(SpeciesbyTrait))
  
  if (appendix){
    
    return(list("sbt"=SpeciesbyTrait,"traitname"=traitname))
  } else {
    return(list("sbt"=SpeciesbyTrait))
  }
  
}

#' @title Create site-by-environment matrix
#' 
#' @description Create site-by-environment (sbe) matrix over the given region.
#' The rows and columns of the matrix represents the sites and environmental
#' variables respectively. 
#' @param rastfile Raster. should be path for storing the worldclim data 
#' (see geodata::worldclim_global) or any other raster data of environmental 
#' variable
#' @param region.sf sf object. The shapefile of the region of study
#' @param res Numeric. The resolution of grid cells to be used. Default is 0.25
#' @param siteID Vector. The cellids of site with occurrences. An output of `sbs`
#' 
#' 
#' 
#' #' @examples
#' \dontrun{
#' rast_path <- "C:/Users/mukht/Documents"
#' KZN_sf<-readRDS(paste0(getwd(),"/KZN_sf.rds"))
#' siteID <- sbs$siteID
#' sbe<-sbeFun(rastfile=rast_path,
#'                  region.sf=KZN_sf,
#'                  res=0.25,
#'                  siteID=siteID,
#'                  appendix=FALSE)
#' }




sbeFun <- function(rastfile,
                   region.sf,
                   res=0.25,
                   siteID,
                   appendix=FALSE){
  # read the rast data if path is given
  if("character" %in% class(rastfile)){
    # Download the WorldClim Bioclimatic variables for the world at a 10 arc-minute resolution
    env <- geodata::worldclim_global(var='bio',
                                     res=2.5, path=rastfile,
                                     version="2.1")
    names(env) <- c('AnnTemp','DiurRange','Isotherm','TempSeas',
                    
                    'MaxTemp','MinTemp','TempRange',
                    
                    'MeanTWQ','MeanTDQ','MeanTWaQ','MeanTCQ','AnnPrec',
                    
                    'PrecWetM','PrecDrM','PrecSeas','PrecWetQ',
                    
                    'PrecDrQ','PrecWaQ','PrecCQ')
  } else if("SpatRaster" %in% class(rastfile)){
    env<-rastfile
  } else { # stop and report if rastfile is not a file path or SpatRaster
    cli::cli_abort(c("{.var rastfile} is not a file path or SpatRaster"))
  }
  
  # Crop environmental variables to extent of of the country's boundary
  env <- terra::crop(env, terra::ext(region.sf))
  # Convert the cropped raster to points, removing NA values
  env_points <- as.data.frame(env, xy = TRUE, na.rm = TRUE)
  
  # Define grid cell for sites and environmental variables
  gridQDS <- terra::rast(region.sf,res=c(res,res), crs="EPSG:4326")
  envQDS <- terra::rast(region.sf,res=c(res,res), crs="EPSG:4326")
  # create raster for environmental data
  for (v in names(env_points[,-c(1,2)])) {
    idw_model <- gstat::gstat(id = v ,formula = as.formula(paste(v,"~1")), data = env_points,
                              locations = ~x+y,nmax = 7, set = list(idp = 2))
    interpolated_raster <- terra::interpolate(gridQDS, idw_model)
    envQDS <- c(envQDS,interpolated_raster[[1]])
  }
  
  #create grid for region
  grid <- region.sf %>%
    sf::st_make_grid(cellsize = c(res,res),
                     offset = c(sf::st_bbox(region.sf)$xmin,
                                sf::st_bbox(region.sf)$ymin)) %>%
    sf::st_sf() %>%
    dplyr::mutate(cellid = dplyr::row_number())
  
  # create site by environment
  sbeM <- as.data.frame(envQDS, xy = TRUE, na.rm = TRUE) %>% 
    sf::st_as_sf(coords = c("x", "y"),crs = 4326) %>% 
    sf::st_join(grid) %>% 
    dplyr::filter(cellid %in% as.integer(siteID)) %>% 
    as.data.frame() %>% 
    dplyr::arrange(cellid) %>% 
    dplyr::select(-any_of(c("cellid","geometry")))
  
  # collect environmental variable name
  variable.name<-colnames(sbeM)
  
  if(appendix){
    return(list("sbe"=sbeM,"variable.name"=variable.name))
  } else {return(list("sbe"=sbeM))}
  
  
}

# >> sbeFUN function
# ------------------------------------------------------------------------------


#' @title Create site-by-environment matrix
#' 
#' @description Create site-by-environment (sbe) matrix over the given region.
#' The rows and columns of the matrix represents the sites and environmental
#' variables respectively. 
#' @param rastfile Raster. should be path for storing the worldclim data 
#' (see geodata::worldclim_global) or any other raster data of environmental 
#' variable
#' @param region.sf sf object. The shapefile of the region of study
#' @param res Numeric. The resolution of grid cells to be used. Default is 0.25
#' @param siteID Vector. The cellids of site with occurrences. An output of `sbs`
#' 
#' 
#' 
#' #' @examples
#' \dontrun{
#' rast_path <- "C:/Users/mukht/Documents"
#' KZN_sf<-readRDS(paste0(getwd(),"/KZN_sf.rds"))
#' siteID <- sbs$siteID
#' sbe<-sbeFun(rastfile=rast_path,
#'                  region.sf=KZN_sf,
#'                  res=0.25,
#'                  siteID=siteID,
#'                  appendix=FALSE)
#' }




sbeFun <- function(rastfile,
                   region.sf,
                   res=0.25,
                   siteID,
                   appendix=FALSE){
  # read the rast data if path is given
  if("character" %in% class(rastfile)){
    # Download the WorldClim Bioclimatic variables for the world at a 10 arc-minute resolution
    env <- geodata::worldclim_global(var='bio',
                                     res=2.5, path=rastfile,
                                     version="2.1")
    names(env) <- c('AnnTemp','DiurRange','Isotherm','TempSeas',
                    
                    'MaxTemp','MinTemp','TempRange',
                    
                    'MeanTWQ','MeanTDQ','MeanTWaQ','MeanTCQ','AnnPrec',
                    
                    'PrecWetM','PrecDrM','PrecSeas','PrecWetQ',
                    
                    'PrecDrQ','PrecWaQ','PrecCQ')
  } else if("SpatRaster" %in% class(rastfile)){
    env<-rastfile
  } else { # stop and report if rastfile is not a file path or SpatRaster
    cli::cli_abort(c("{.var rastfile} is not a file path or SpatRaster"))
  }
  
  # Crop environmental variables to extent of of the country's boundary
  env <- terra::crop(env, terra::ext(region.sf))
  # Convert the cropped raster to points, removing NA values
  env_points <- as.data.frame(env, xy = TRUE, na.rm = TRUE)
  
  # Define grid cell for sites and environmental variables
  gridQDS <- terra::rast(region.sf,res=c(res,res), crs="EPSG:4326")
  envQDS <- terra::rast(region.sf,res=c(res,res), crs="EPSG:4326")
  # create raster for environmental data
  for (v in names(env_points[,-c(1,2)])) {
    idw_model <- gstat::gstat(id = v ,formula = as.formula(paste(v,"~1")), data = env_points,
                              locations = ~x+y,nmax = 7, set = list(idp = 2))
    interpolated_raster <- terra::interpolate(gridQDS, idw_model)
    envQDS <- c(envQDS,interpolated_raster[[1]])
  }
  
  #create grid for region
  grid <- region.sf %>%
    sf::st_make_grid(cellsize = c(res,res),
                     offset = c(sf::st_bbox(region.sf)$xmin,
                                sf::st_bbox(region.sf)$ymin)) %>%
    sf::st_sf() %>%
    dplyr::mutate(cellid = dplyr::row_number())
  
  # create site by environment
  sbeM <- as.data.frame(envQDS, xy = TRUE, na.rm = TRUE) %>% 
    sf::st_as_sf(coords = c("x", "y"),crs = 4326) %>% 
    sf::st_join(grid) %>% 
    dplyr::filter(cellid %in% as.integer(siteID)) %>% 
    as.data.frame() %>% 
    dplyr::arrange(cellid) %>% 
    dplyr::select(-any_of(c("cellid","geometry")))
  
  # collect environmental variable name
  variable.name<-colnames(sbeM)
  
  if(appendix){
    return(list("sbe"=sbeM,"variable.name"=variable.name))
  } else {return(list("sbe"=sbeM))}
  
  
}


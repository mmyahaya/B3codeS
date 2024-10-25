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





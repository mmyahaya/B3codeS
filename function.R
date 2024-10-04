library(dplyr)
library(readr)
library(tidyverse)
library(rgbif) # for occ_download
library(terra)
library(sf)
library(rtry) # for processing try data
library(rasterVis)
library(rWCVP)
library(rWCVPdata)

taxaFun <- function(taxa,country.sf,limit=500, ref=NULL,country='ZA',res=0.25){
  
  grid <- country.sf %>%
    sf::st_make_grid(cellsize = c(res,res),
                     offset = c(sf::st_bbox(country.sf)$xmin,
                                sf::st_bbox(country.sf)$ymin)) %>%
    sf::st_sf() %>%
    dplyr::mutate(cellid = dplyr::row_number())
  
  # download taxa (target taxa) if the scientific name is given as character
  if("character" %in% class(taxa)){
    taxa.gbif_download = rgbif::occ_data(scientificName=taxa, # download data from gbif
                                  country=country,
                                  hasCoordinate=TRUE,
                                  hasGeospatialIssue=FALSE,
                                  limit = limit)

    taxa.df = as.data.frame(taxa.gbif_download$data) #extract data from the downloaded file
  } else if("data.frame" %in% class(taxa)){ #check if data fame contains the required columns
    if(any(!c("decimalLatitude","decimalLongitude",
              "species","speciesKey","coordinateUncertaintyInMeters","dateIdentified","year") %in% colnames(taxa))){
      requiredcol<-c("decimalLatitude","decimalLongitude","species","coordinateUncertaintyInMeters","dateIdentified","year")
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
    dplyr::mutate(occurrenceStatus=1)
    
    
    taxa_cube<-b3gbi::process_cube(taxa.sf,grid_type = "custom",
                                    cols_cellCode = "cellid", cols_year = "year",
                                    cols_occurrences = "occurrenceStatus",
                                    cols_species = "species",cols_speciesKey = "speciesKey",
                                    cols_minCoordinateUncertaintyInMeters = "coordinateUncertaintyInMeters")
  # download reference taxa if provided or use the target taxa if otherwise
  if(!is.null(ref)){
    if("character" %in% class(ref)){
      ref.gbif_download = rgbif::occ_data(scientificName=ref, # download data from gbif
                                   country='ZA',
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
        dplyr::mutate(occurrenceStatus=1)
      
      
        ref_cube<-b3gbi::process_cube(ref.sf,grid_type = "custom",
                                     cols_cellCode = "cellid", cols_year = "year",
                                     cols_occurrences = "occurrenceStatus",
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
        dplyr::mutate(occurrenceStatus=1)
      
      
      ref_cube<-b3gbi::process_cube(ref.sf,grid_type = "custom",
                                    cols_cellCode = "cellid", cols_year = "year",
                                    cols_occurrences = "occurrenceStatus",
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

# Specie by site
sbsFun <- function(taxa_cube,country.sf,res){

  # extract unique species name from GBIF occurrence data
  species_list<-sort(unique(taxa_cube$taxa$data$scientificName))
  
  sbsM<-taxa_cube$taxa$data %>%
    dplyr::select(scientificName,cellCode,obs) %>%
    dplyr::group_by(scientificName,cellCode) %>%
    dplyr::summarise(across(obs, sum), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
    dplyr::arrange(cellCode) %>% 
    tibble::column_to_rownames(var = "cellCode") %>% 
    as.matrix()
  
  #create grid for region
  grid <- country.sf %>%
    sf::st_make_grid(cellsize = c(res,res),
                     offset = c(sf::st_bbox(country.sf)$xmin,
                                sf::st_bbox(country.sf)$ymin)) %>%
    sf::st_sf() %>%
    dplyr::mutate(cellid = dplyr::row_number())
  
  # get coordinates of the occurrence sites
  coords <- sf::st_coordinates(sf::st_centroid(grid))
  coords <- coords[as.integer(rownames(sbsM)),]
  colnames(coords)<-c("Longitude","Latitude")
 
  colnames(sbsM)<-NULL # remove column names
  # create binary matrix
  sbsM.binary<-sbsM
  sbsM.binary[sbsM.binary>0]<-1
  
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
      as.matrix()
    
    colnames(sbsM.ref)<-NULL # remove column names
  }


  return(list("sbs"=sbsM,"sbs.ref"=sbsM.ref,"sbs.binary"=sbsM.binary,
              "coords"=coords,"species_list"=species_list,
              "siteID"=rownames(sbsM),"site_unc"=site_uncertainty))
}


#specie by environment
sbeFun <- function(rastfile,country.sf,res=0.25,siteID){

  # read the rastfile if path is given
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
  env <- terra::crop(env, ext(country.sf))
  # Convert the cropped raster to points, removing NA values
  env_points <- as.data.frame(env, xy = TRUE, na.rm = TRUE)
  
  # Define grid cell for sites and environmental variables
  gridQDS <- terra::rast(country.sf,res=c(res,res), crs="EPSG:4326")
  envQDS <- terra::rast(country.sf,res=c(res,res), crs="EPSG:4326")
  # create raster for environmental data
  for (v in names(env_points[,-c(1,2)])) {
    idw_model <- gstat::gstat(id = v ,formula = as.formula(paste(v,"~1")), data = env_points,
                              locations = ~x+y,nmax = 7, set = list(idp = 2))
    interpolated_raster <- terra::interpolate(gridQDS, idw_model)
    envQDS <- c(envQDS,interpolated_raster[[1]])
  }
  
  #create grid for region
  grid <- country.sf %>%
    sf::st_make_grid(cellsize = c(res,res),
                     offset = c(sf::st_bbox(country.sf)$xmin,
                                sf::st_bbox(country.sf)$ymin)) %>%
    sf::st_sf() %>%
    dplyr::mutate(cellid = dplyr::row_number())
  
  # create site by environment
  sbeM <- as.data.frame(envQDS, xy = TRUE, na.rm = TRUE) %>% 
    sf::st_as_sf(coords = c("x", "y"),crs = 4326) %>% 
    sf::st_join(grid) %>% 
    dplyr::filter(cellid %in% as.integer(siteID)) %>% 
    as.data.frame() %>% 
    dplyr::arrange(cellid) %>% 
    dplyr::select(-any_of(c("cellid","geometry"))) %>% 
    as.matrix()

  
  # collect environmental variable name
  variable.name<-colnames(sbeM)
  colnames(sbeM)<-NULL
  return(list("sbe"=sbeM,"variable.name"=variable.name))
}


#species by trait
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
    drop_na(TraitID) %>%
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

dataGEN = function(taxa,country.sf,country='ZA',limit=500,ref=NULL,
                   res=0.25,tryfile,rastfile){
  taxa_cube <- taxaFun(taxa=taxa,country.sf=country.sf,limit=limit, ref=ref,country=country,res=res)
  sbs <- sbsFun(taxa_cube = taxa_cube,country.sf = country.sf, res = res)
  sbt <- sbtFun(tryfile = tryfile,taxa_cube = taxa_cube)
  siteID <- sbs$siteID
  sbe <- sbeFun(rastfile=rastfile,country.sf=country.sf,res=res,siteID=siteID)
  return(list("sbs"=sbs,"sbt"=sbt,"sbe"=sbe))
}





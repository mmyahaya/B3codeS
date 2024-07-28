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
library(stringr)

taxaFun <- function(taxa,limit=500, ref=NULL,country='ZA'){
  
  # download taxa (target taxa) if the scientific name is given as character
  if("character" %in% class(taxa)){
    taxa.gbif_download = rgbif::occ_data(scientificName=taxa, # download data from gbif
                                  country=country,
                                  hasCoordinate=TRUE,
                                  hasGeospatialIssue=FALSE,
                                  limit = limit)

    taxa.df = as.data.frame(taxa.gbif_download$data) #extract data from the downloaded file
  } else if("data.frame" %in% class(taxa)){ #check if data fame contains the required colums
    if(any(!c("decimalLatitude","decimalLongitude",
              "species","dateIdentified") %in% colnames(taxa))){
      requiredcol<-c("decimalLatitude","decimalLongitude","species","dateIdentified")
      missingcol<-requiredcol[!c("decimalLatitude","decimalLongitude","species","dateIdentified") %in% colnames(taxa)]
      cli::cli_abort(c("{missingcol} is/are not in the {.var taxa} column ", 
                      "x" = "{.var taxa} should be a data of GBIF format "))
    }
    # take taxa data frame if accurate
    taxa.df<-taxa
  } else { # stop and report if taxa is not a scientific name or dataframe
    cli::cli_abort(c("{.var taxa} is not a character or dataframe"))
  }

  
  taxa.df = taxa.df %>%
    dplyr::select(decimalLatitude,decimalLongitude,
                  species,dateIdentified) %>% #select occurrence data
    filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
    mutate(dateIdentified = as.Date(dateIdentified)) # convert date to date format
  taxa.sf<-sf::st_as_sf(taxa.df,coords = c("decimalLongitude", "decimalLatitude"),
                    crs = 4326)
  # download reference taxa if provided or use the target taxa if otherwise
  if(!is.null(ref)){
    if("character" %in% class(ref)){
      ref.gbif_download = rgbif::occ_data(scientificName=ref, # download data from gbif
                                   country='ZA',
                                   hasCoordinate=TRUE,
                                   hasGeospatialIssue=FALSE,
                                   limit = limit)
      
      ref.df = as.data.frame(ref.gbif_download$data)
      
      ref.df = ref.df %>%
        dplyr::select(decimalLatitude,decimalLongitude,
                      species,dateIdentified) %>% #select occurrence data
        filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
        mutate(dateIdentified = as.Date(dateIdentified)) # convert date to date format
      ref.sf<-sf::st_as_sf(ref.df,coords = c("decimalLongitude", "decimalLatitude"),
                       crs = 4326)
    } else if("data.frame" %in% class(ref)){ # check columns of ref if provided
      if(any(!c("decimalLatitude","decimalLongitude",
                "species","dateIdentified") %in% colnames(ref))){
        requiredcol<-c("decimalLatitude","decimalLongitude","species","dateIdentified")
        missingcol<-requiredcol[!c("decimalLatitude","decimalLongitude","species","dateIdentified") %in% colnames(ref)]
        cli::cli_abort(c("{missingcol} is/are not in the {.var ref} column ", 
                         "x" = "{.var ref} should be a data of GBIF format "))
      }
      ref = ref %>%
        dplyr::select(decimalLatitude,decimalLongitude,
                      species,dateIdentified) %>% #select occurrence data
        filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
        mutate(dateIdentified = as.Date(dateIdentified)) # convert date to date format
      ref.sf<-sf::st_as_sf(ref,coords = c("decimalLongitude", "decimalLatitude"),
                       crs = 4326)
      
      } else { # stop and report if taxa is not a scientific name or dataframe
        cli::cli_abort(c("{.var ref} is not a character or dataframe"))
      }

  } else {
    ref.sf=taxa.sf
  }
  return(list("taxa"=taxa.sf,"ref"=ref.sf))
}

# Specie by site
sbsFun <- function(taxa.sf,country.sf,res=0.25){
  # extract unique species name from GBIF occurrence data
  uN<-sort(unique(taxa.sf$taxa$species))
  # Create grid cells with extent of country shapefile and layers for siteID and species
  gridQDS = rast(country.sf,res=c(res,res), crs="EPSG:4326",nlyrs=length(uN)+1)
  # specify name for each layers for site ID and species
  names(gridQDS)<-c("siteID",uN)
  # Assign ID for each cell
  gridQDS[["siteID"]]<-1:ncell(gridQDS)
  # create layer for species occurrence in each cell
  for(n in uN){
    # create raster of species
    speciesQDS = terra::rasterize(dplyr::filter(taxa.sf$taxa, species==n),
                           gridQDS,
                           field=1,
                           fun="sum",
                           background = 0)
    # insert occurrence layer for each species to it assigned layer
    gridQDS[[n]] <- speciesQDS[]
  }
 # mask grid cells to country shape file
  gridQDS = mask(gridQDS, country.sf)

  # create data frame of site by species
  sbs <- as.data.frame(gridQDS)
  # select sites which have occurrences
  sbs<-sbs[rowSums(sbs[,-1])!=0,]
  # collect site ID
  siteID<-sbs$siteID
  
  
  # get coordinates of the occurrence sites
  coords <- terra::xyFromCell(gridQDS, sbs$siteID)
  colnames(coords)<-c("Longitude","Latitude")
  sbsM<-as.matrix(sbs[,-1]) # remove siteID and convert to matrix
  colnames(sbsM)<-NULL # remove column names
  # create binary matrix
  sbsM.binary<-sbsM
  sbsM.binary[sbsM.binary>0]<-1
  
  # compute sbs for ref if it is different from taxa
  if(identical(taxa.sf$taxa,taxa.sf$ref)){
    sbsM.ref<-sbsM
  } else{
    # extract unique species name from GBIF occurrence data
    uN.ref<-sort(unique(taxa.sf$ref$species))
    # Create grid cells with extent of country shapefile and layers for siteID and species
    gridQDS = terra::rast(country.sf,res=c(res,res), crs="EPSG:4326",
                   nlyrs=length(uN.ref)+1)
    # specify name for each layers of site ID and individual species
    names(gridQDS)<-c("siteID",uN.ref)
    # Assign ID for each cell
    gridQDS[["siteID"]]<-1:ncell(gridQDS)
    # create layer for species occurrence in each cell
    for(n in uN.ref){
      # create raster of species
      speciesQDS = rasterize(dplyr::filter(taxa.sf$ref, species==n),
                             gridQDS,
                             field=1,
                             fun="sum",
                             background = 0)
      # insert occurrence layer for each species to it assigned layer
      gridQDS[[n]] <- speciesQDS[]
    }
    # mask grid cells to country shape file
    gridQDS = mask(gridQDS, country.sf)

    # create data frame of site by species
    sbs.ref <- as.data.frame(gridQDS)
    # select sites base on taxa occurrence site
    sbs.ref <- sbs.ref[sbs$siteID,]
   

    sbsM.ref<-as.matrix(sbs.ref[,-1]) # remove siteID and convert to matrix
    colnames(sbsM.ref)<-NULL # remove column names
  }


  return(list("sbs"=sbsM,"sbs.ref"=sbsM.ref,"sbs.binary"=sbsM.binary,
              "coords"=coords,"species.name"=uN,"siteID"=sbs$siteID))
}


#specie by environment
sbeFun<- function(rastfile,taxa.sf,country.sf,res=0.25,fun="mean", siteID){

  # read the rastfile if path is given
  if("character" %in% class(rastfile)){
    # Download the WorldClim Bioclimatic variables for the world at a 10 arc-minute resolution
    env = geodata::worldclim_global(var='bio',
                                        res=2.5, path=rastfile,
                                        version="2.1")
    names(env) = c('AnnTemp','DiurRange','Isotherm','TempSeas',
                   
                   'MaxTemp','MinTemp','TempRange',
                       
                       'MeanTWQ','MeanTDQ','MeanTWaQ','MeanTCQ','AnnPrec',
                       
                       'PrecWetM','PrecDrM','PrecSeas','PrecWetQ',
                   
                   'PrecDrQ','PrecWaQ','PrecCQ')
  } else if("SpatRaster" %in% class(rastfile)){
    env<-rastfile
  } else { # stop and report if rastfile is not a file path or SpatRaster
    cli::cli_abort(c("{.var rastfile} is not a file path or SpatRaster"))
  }

  # Crop Bioclimatic variables to extent of of the country's boundary
  env = crop(env, ext(country.sf))
  # Extract environmental data for species occurrence point
  taxa.env.sf = extract(env, taxa.sf$taxa, bind=TRUE)
  
  # Define grid cell for sites
  gridQDS = rast(country.sf,res=c(res,res), crs="EPSG:4326")
  
  # create raster for environmental data
  envQDS = rasterize(taxa.env.sf,
                     gridQDS,
                     field=names(env),
                     fun=fun,
                     background = NA)
 

  # extract site by environment from the QDS layers
  sitebyEnv = as.data.frame(envQDS[])
  # select occurrence sites
  sitebyEnv<-sitebyEnv[siteID,]
  # create species by site matrix and drop siteID
  sbeM <-as.matrix(sitebyEnv)
  variable.name<-colnames(sbeM)
  colnames(sbeM)<-NULL
  return(list("sbe"=sbeM,"variable.name"=variable.name))
}


#species by trait
sbtFun<-function(tryfile,taxa.sf){
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
  uN<-sort(unique(taxa.sf$taxa$species))

  SpeciesbyTrait<-trydata %>%
    # drop rows which contains no trait
    drop_na(TraitID) %>%
    # select species name, trait and trait value
    select(AccSpeciesName,TraitID,OrigValueStr) %>%
    #group by Species and trait
    group_by(AccSpeciesName,TraitID) %>%
    #choose the first trait value if there are multiples trait for a species
    summarise(across(OrigValueStr, first), .groups = "drop") %>%
    # reshape to wide format to have specie by trait dataframe
    pivot_wider(names_from = TraitID, values_from = OrigValueStr) %>%
    # select species that are only present in gbif data
    filter(AccSpeciesName %in% uN) %>%
    # convert species names to row names
    column_to_rownames(var = "AccSpeciesName")
  # add the rows of the remaining species without traits from TRY
  na.df<-as.data.frame(matrix(NA,nrow = length(setdiff(uN,rownames(SpeciesbyTrait))),
                               ncol = ncol(SpeciesbyTrait)))
    row.names(na.df)<-setdiff(uN,rownames(SpeciesbyTrait))
    names(na.df)<-names(SpeciesbyTrait) 
    SpeciesbyTrait<-rbind(SpeciesbyTrait,na.df)
    #sort according to unique species vector
    SpeciesbyTrait<-SpeciesbyTrait[uN,]
    
    # Download WCVP for native taxa in area of interest
    native_list <- rWCVP::wcvp_checklist(taxon = stringr::word(uN[1],1), taxon_rank = "genus") %>%
      filter(area_code_l3 %in% get_wgsrpd3_codes("South Africa")) %>%
      filter((accepted_name %in% uN) & occurrence_type=="native")

    # create taxa list
    taxa_list<-data.frame(taxon=uN)

    # create new dataframe with introduction status
    taxa_list_status<-taxa_list%>%
      mutate(introduction_status = ifelse(taxon%in%native_list$accepted_name,
                                          "native","introduced"))
    # add introduction status to trait column
    SpeciesbyTrait$introduction_status<-taxa_list_status$introduction_status
    
    
    sbtM<-as.matrix(SpeciesbyTrait)
   #collect traitID
    trait<-colnames(sbtM[,-ncol(sbtM)])
    #remove column and row names 
    rownames(sbtM)<-NULL
    colnames(sbtM)<-NULL
    
    
    #extract the trait names
    traitname<-trydata %>%
      # drop rows which contains no trait
      drop_na(TraitID) %>%
      # select the trait ID and trait name
      select(TraitID,TraitName) %>%
      group_by(TraitID) %>%
      summarise(across(TraitName, first), .groups = "drop") %>%
      column_to_rownames("TraitID")
    #create trait name to align with the sbt column
    traitname<-c(traitname[trait,],"Introduction status")
    traitname<-data.frame('TraitID'=c(trait,"Introduction status"),'TraitName'=traitname)
  return(list("sbt"=sbtM,"traitname"=traitname))

}

dataGEN = function(taxa,country.sf,country='ZA',limit=500,ref=NULL,
                   res=0.25,tryfile,rastfile,fun="mean"){
  taxa.sf <- taxaFun(taxa = taxa,limit = limit,ref = ref, country = country)
  sbs <- sbsFun(taxa.sf = taxa.sf,country.sf = country.sf,res = res )
  sbt <- sbtFun(tryfile = tryfile,taxa.sf = taxa.sf)
  siteID <- sbs$siteID
  sbe <- sbeFun(rastfile = rastfile, taxa.sf=taxa.sf, country.sf = country.sf,res = res,
  fun=fun,siteID = siteID)
  return(list("sbs"=sbs,"sbt"=sbt,"sbe"=sbe))
}





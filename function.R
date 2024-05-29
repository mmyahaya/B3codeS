library(dplyr)
library(readr)
library(tidyverse)
library(rgbif) # for occ_download
library(terra)
library(sf)
library(rtry) # for processing try data
library(rasterVis)
library(lubridate)
rsa_country_sf = st_read("C:/Users/26485613/OneDrive - Stellenbosch University/Downloads/Code_Data/Code_Data/boundary_south_africa_land_geo.shp")
# "C:/Users/mukht/Documents/boundary_SA/boundary_south_africa_land_geo.shp"

taxaFun <- function(taxa,limit=500, ref=NULL,country='ZA'){
  
  # download taxa if scientific name is given
  if("character" %in% class(taxa)){
    taxa.gbif_download = occ_data(scientificName=taxa, # download data from gbif
                                  country=country,
                                  hasCoordinate=TRUE,
                                  hasGeospatialIssue=FALSE,
                                  limit = limit)

    taxa.df = as.data.frame(taxa.gbif_download$data) #extract data from the downloaded file
  } else if("data.frame" %in% class(taxa)){
    if(any(!c("decimalLatitude","decimalLongitude",
              "species","dateIdentified") %in% colnames(taxa))){
      requiredcol<-c("decimalLatitude","decimalLongitude","species","dateIdentified")
      missingcol<-requiredcol[!c("decimalLatitude","decimalLongitude","species","dateIdentified") %in% colnames(taxa)]
      cli::cli_abort(c("{missingcol} is/are not in the {.var taxa} column ", 
                      "x" = "{.var taxa} should be a data of GBIF format "))
    }
    
    taxa.df<-taxa
  } else { # stop and report if taxa is not a scientific name or dataframe
    cli::cli_abort(c("{.var taxa} is not a character or dataframe"))
  }

  
  taxa.df = taxa.df %>%
    dplyr::select(decimalLatitude,decimalLongitude,
                  species,dateIdentified) %>% #select occurrence data
    filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
    mutate(dateIdentified = as.Date(dateIdentified)) # convert date to date format
  taxa.sf<-st_as_sf(taxa.df,coords = c("decimalLongitude", "decimalLatitude"),
                    crs = 4326)
  # download  or use ref dataframe taxa if provided
  if(!is.null(ref)){
    if("character" %in% class(ref)){
      ref.gbif_download = occ_data(scientificName=ref, # download data from gbif
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
      ref.sf<-st_as_sf(ref.df,coords = c("decimalLongitude", "decimalLatitude"),
                       crs = 4326)
    } else if("data.frame" %in% class(ref)){
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
      ref.sf<-st_as_sf(ref,coords = c("decimalLongitude", "decimalLatitude"),
                       crs = 4326)
      
      } else { # stop and report if taxa is not a scientific name or dataframe
        cli::cli_abort(c("{.var ref} is not a character or dataframe"))
      }

  } else {
    ref.sf=taxa.sf
  }


  return(list("taxa"=taxa.sf,"ref"=ref.sf))
}

# Save taxa.sf to drive to avoid redownloading same taxa with same limit next time
#  write.csv(taxa.sf,"taxa.sf.csv",row.names = FALSE)
#add res
sbsFun <- function(taxa.sf,country.sf,res=0.25){
  # extract unique species name from GBIF occurrence data
  uN<-sort(unique(taxa.sf$taxa$species))
  # Create grid cells with extent of country shapefile and layers for siteID and species
  gridQDS = rast(country.sf,res=c(res,res), crs="EPSG:4326",nlyrs=length(uN)+1)
  # specify name for each layers of site ID and individual species
  names(gridQDS)<-c("siteID",uN)
  # Assign ID for each cell
  gridQDS[["siteID"]]<-1:ncell(gridQDS)
  # create layer for species occurrence in each cell
  for(n in uN){
    # create raster of species
    speciesQDS = rasterize(dplyr::filter(taxa.sf$taxa, species==n),
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
  sbs<-drop_na(sbs,siteID)
  # get coordinates of sites on the country shapefile
  coords <- xyFromCell(gridQDS, sbs$siteID)
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
    gridQDS = rast(country.sf,res=c(res,res), crs="EPSG:4326",
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
    sbs.ref<-drop_na(sbs.ref,siteID)

    sbsM.ref<-as.matrix(sbs.ref[,-1]) # remove siteID and convert to matrix
    colnames(sbsM.ref)<-NULL # remove column names
  }


  return(list("sbs"=sbsM,"sbs.ref"=sbsM.ref,"sbs.binary"=sbsM.binary,
              "coords"=coords,"species.name"=uN))
}



sbeFun<- function(rastfile,country.sf,res=0.25){

  # read the rastfile if path is given
  if("character" %in% class(rastfile)){
    # Download the WorldClim Bioclimatic variables for the world at a 10 arc-minute resolution
    env = geodata::worldclim_global(var='bio',
                                        res=10, path=rastfile,
                                        version="2.1") # Set your own path directory
  } else if("SpatRaster" %in% class(rastfile)){
    env<-rastfile
  } else { # stop and report if rastfile is not a file path or SpatRaster
    cli::cli_abort(c("{.var rastfile} is not a file path or SpatRaster"))
  }


  # Crop Bioclimatic variables to extent of of the country's boundary
  env = crop(env, country.sf)

  gridQDS = rast(country.sf,res=c(res,res), crs="EPSG:4326")


  # Transfer values from rastfile data to QDS
  envQDS<-resample(env,gridQDS) # bilinear interpolation
  envQDS[["siteID"]]<-1:ncell(envQDS)


  # mask envQDS to the country map
  envQDS<-mask(envQDS,country.sf)

  # extract site by environment from the bioQDS layers
  {sitebyEnv <- as.data.frame(envQDS[])
    sitebyEnv <- drop_na(sitebyEnv,siteID)
    sbeM<-as.matrix(dplyr::select(sitebyEnv, !siteID))
    variable.name<-colnames(sbeM)
    colnames(sbeM)<-NULL}
  return(list("sbe"=sbeM,"variable.name"=variable.name))
}




sbtFun<-function(tryfile,taxa.sf){
  # read the try data if path is given
  if("character" %in% class(tryfile)){
    trydata<-rtry_import(
      input=tryfile,
      separator = "\t",
      encoding = "Latin-1",
      quote = "",
      showOverview = TRUE
    )
  } else if("data.frame" %in% class(tryfile)){
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
  uN<-sort(unique(taxa.sf$species))

  SpeciesbyTrait<-trydata %>%
    # drop rows which contains no trait
    drop_na(TraitID) %>%
    # select species name, trait and trait value
    select(AccSpeciesName,TraitID,TraitName,OrigValueStr) %>%
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
    names(na.df)<-names(SpeciesbyTrait) # column names
    SpeciesbyTrait<-rbind(SpeciesbyTrait,na.df)
    SpeciesbyTrait<-SpeciesbyTrait[order(rownames(SpeciesbyTrait)),]
    sbtM<-as.matrix(SpeciesbyTrait)
    rownames(sbtM)<-NULL
    colnames(sbtM)<-NULL

    traitname<-trydata %>%
      # drop rows which contains no trait
      drop_na(TraitID) %>%
      # select species name, trait and trait value
      select(TraitID,TraitName) %>%
      group_by(TraitID) %>%
      summarise(across(TraitName, first), .groups = "drop") %>%
      column_to_rownames("TraitID")
    #create Trait name to align with the sbt column
    traitname<-traitname[colnames(SpeciesbyTrait),]
    traitname<-data.frame('TraitID'=colnames(SpeciesbyTrait),'TraitName'=traitname)
  return(list("sbt"=sbtM,"traitname"=traitname))

}

dataGEN = function(taxa,country.sf,country='ZA',limit=500,ref=NULL,
                   res=0.25,tryfile,rastfile){
  taxa.sf <- taxaFun(taxa = taxa,limit = limit,ref = ref, country = country)
  sbs <- sbsFun(taxa.sf = taxa.sf,country.sf = country.sf,res = res )
  sbt <- sbtFun(tryfile = tryfile,taxa.sf = taxa.sf$taxa)
  sbe <- sbeFun(rastfile = rastfile,country.sf = country.sf,res = res)
  return(list("sbs"=sbs,"sbt"=sbt,"sbe"=sbe))
}


datalist<-dataGEN(taxa=taxa_df,country.sf=rsa_country_sf,ref=ref.data1,
                  tryfile=try33852,rastfile=precdata)

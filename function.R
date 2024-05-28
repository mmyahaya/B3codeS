library(dplyr)
library(readr)
library(tidyverse)
library(rgbif) # for occ_download
library(terra)
library(sf)
library(rtry) # for processing try data
library(rasterVis)
library(lubridate)
rsa_country_sf = st_read("C:/Users/mukht/Documents/boundary_SA/boundary_south_africa_land_geo.shp")
# add reference and get it sbs
taxaFun <- function(taxa,limit){
  gbif_download = occ_data(scientificName=taxa, # download data from gbif
                           country='ZA',
                           hasCoordinate=TRUE,
                           hasGeospatialIssue=FALSE,
                           limit = limit)
  
  taxa.df = as.data.frame(gbif_download$data) #extract data from the downloaded file
  
  
  taxa.df = taxa.df %>%
    dplyr::select(decimalLatitude,decimalLongitude,
                  species,dateIdentified) %>% #select occurrence data
    filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
    mutate(dateIdentified = as.Date(dateIdentified)) # convert date to date format
  taxa.sf<-st_as_sf(taxa.df,coords = c("decimalLongitude", "decimalLatitude"),
                    crs = 4326)
  return(taxa.sf)
}

# Save taxa.sf to drive to avoid redownloading same taxa with same limit next time
#  write.csv(taxa.sf,"taxa.sf.csv",row.names = FALSE)
#add res
sbsFun <- function(taxa.sf,country.shp){
  # extract unique species name from GBIF occurrence data
  uN<-sort(unique(taxa.sf$species))
  # Create grid cells with extent of country shapefile and layers for siteID and species
  gridQDS = rast(country.shp,res=c(0.25,0.25), crs="EPSG:4326",nlyrs=length(uN)+1)
  # specify name for each layers of site ID and individual species
  names(gridQDS)<-c("siteID",uN)
  # Assign ID for each cell
  gridQDS[["siteID"]]<-1:ncell(gridQDS)
  # create layer for species occurrence in each cell
  system.time(for(n in uN){
    # create raster of species
    speciesQDS = rasterize(dplyr::filter(taxa.sf, species==n),
                           gridQDS,
                           field=1,
                           fun="sum",
                           background = 0)
    # insert occurrence layer for each species to it assigned layer
    gridQDS[[n]] <- speciesQDS[]
  })
 # mask grid cells to country shape file
  gridQDS = mask(gridQDS, country.shp)
  
  # create data frame of site by species
  sbs <- as.data.frame(gridQDS)
  sbs<-drop_na(sbs,siteID)
  # get coordinates of sites on the country shapefile
  coords <- xyFromCell(gridQDS, sbs$siteID)
  colnames(coords)<-c("Longitude","Latitude")
  sbsM<-as.matrix(sbs[,-1])
  colnames(sbsM)<-NULL
  sbsM.binary<-sbsM
  sbsM.binary[sbsM.binary>0]<-1
  return(list("sbs"=sbsM,"sbs.binary"=sbsM.binary,"coords"=coords,"species.name"=uN))
}



sbeFun<- function(path,country.shp){
  # Download the WorldClim Bioclimatic variables for the world at a 10 arc-minute resolution
  bio_10m = geodata::worldclim_global(var='bio',
                                      res=10, path=path,
                                      version="2.1") # Set your own path directory
  
  # Crop Bioclimatic variables to extent of of the country's boundary
  bio_10m = crop(bio_10m, rsa_country_sf)
 
  gridQDS = rast(country.shp,res=c(0.25,0.25), crs="EPSG:4326")
  
  
  # Transfer values from worldclim raster data to QDS
  bioQDS<-resample(bio_10m,gridQDS) # bilinear interpolation 
  bioQDS[["siteID"]]<-1:ncell(bioQDS)
  
  
  # mask bioQDS to the country map
  bioQDS<-mask(bioQDS,country.shp)
  
  # extract site by environment from the bioQDS layers
  {sitebyEnv <- as.data.frame(bioQDS[])
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
    trydata<-tryfile
  } else { # stop and report if tryfile is not a file path or dataframe
    cli::cli_abort(c("{.var tryfile} is not a file path or dataframe"))
    }
  
  if(any(!c("AccSpeciesName","TraitID","TraitName","OrigValueStr") %in% colnames(tryfile))){
    requiredcol<-c("AccSpeciesName","TraitID","TraitName","OrigValueStr")
    missingcol<-requiredcol[!c("AccSpeciesName","TraitID","TraitName","OrigValueStr") %in% colnames(tryfile)]
    cli::cli_abort(c("{missingcol} is/are not in the {.var tryfile} column "))
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


taxa.sf <- taxaFun('Acacia',500)
sbs<-sbsFun(taxa.sf = taxa.sf, country.shp = rsa_country_sf)

path = "C:/Users/mukht/Documents" #path for worldclim

sbe<-sbeFun(path = path, country.shp = rsa_country_sf )

try_path<-"33852.txt" # path for trydata

sbt<-sbtFun(tryfile = trytest,taxa.sf = taxa.sf)

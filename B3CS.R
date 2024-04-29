

library(dplyr)
library(readr)
library(tidyverse)
library(rgbif) # for occ_download
library(terra)
library(sf)
library(rtry) # for processing try data
library(rasterVis)
##### gbif data ####
taxa = 'Tracheophyta' # scientific name

gbif_download = occ_data(scientificName=taxa, # download data from gbif
                         country='ZA',
                         hasCoordinate=TRUE,
                         hasGeospatialIssue=FALSE,
                         limit = 2000)

taxa.df = as.data.frame(gbif_download$data) #extract data from the downloaded file


taxa.occ = taxa.df %>%
  dplyr::select(key,decimalLatitude,decimalLongitude,
                species,dateIdentified) %>% #select occurrence data
  filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
  mutate(dateIdentified = as.Date(dateIdentified)) # convert date to date format

taxa.sf<-st_as_sf(taxa.occ,coords = c("decimalLongitude", "decimalLatitude"),
                  crs = 4326) # convert long and lat point to geometry
taxa.sf$count=1
# Define a grid over spatial extend
gridQDS = rast(ext(taxa.sf),res=c(0.25,0.25), crs="EPSG:4326")

countQDS = rasterize(taxa.sf,
                     gridQDS,
                     field='count',
                     fun=sum,
                     background = 0)


levelplot(countQDS) # improve plots
plot(grid)
plot(taxa.sf[1], add=TRUE)
##### TRY data ####
library(rtry)

path_to_data <- system.file("testdata", "data_TRY_15160.txt", package = "rtry")
path_to_data

TRYdata1 <- rtry_import(path_to_data)
input_path<-"C://Users//26485613//OneDrive - Stellenbosch University//Documents//Practice space"

tryFULL.df<-rtry_import(input_path,
                         separator = ",",
                         encoding = "UTF-8",
                         quote = "\"",
                         showOverview = TRUE)




##### Junks ####
taxa.df %>%
  drop_na(species) %>%
  count(species, sort = TRUE) %>%
  filter(n>3)




plot(countQDS)# plot grid cell
plot(taxa.sf[1], add=TRUE) #add occurence points

uN<-unique(taxa.sf$species)

##### TRY data ####


input_path<-"C://Users//26485613//OneDrive - Stellenbosch University//Documents//Practice space//33312.txt"
# import data
tryFULL.df<-rtry_import(
  input=input_path,
  separator = "\t",
  encoding = "Latin-1",
  quote = "",
  showOverview = TRUE
)

uN.df<-as.data.frame(uN)

TRYcat <- readxl::read_excel("TRYcat.xlsx", sheet = 1)
# create data of species with traits
CATtrait.df<-inner_join(uN.df,TRYcat[,1:22],
                        by=join_by("uN"=="AccSpeciesName"))

##### Junks ####
taxa.df %>%
  drop_na(species) %>%
  count(species, sort = TRUE) %>%
  filter(n>3)

summary(CATtrait.df)

CATtrait.df<- CATtrait.df %>%
  mutate(decade=ifelse(dateIdentified>as.Date("2024-02-28"),"First","Second"))


dataGEN = function(arg1,TaxaName..){


  # prepare data from other cubes and accessible datasets


  # spatial polygons, shapefiles, squares (corner coordinates) for spatial extent
  # spatial resolution, specified and default
  # alien status of identified species list
  # [we need to make a decision on species list completeness (this needs to consider other cubes and SDM cubes, etc.)]
  # temporal range
  # option = 1, 2; 1 for dissim() and 2 for invasib()
  # Return, site by species, site by xyt, site by env, site by site distance; species by trait, species by species phylogenetic distance, matrices (and default setting if null)
  return(arg1)
}





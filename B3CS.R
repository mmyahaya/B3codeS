

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


rasterVis::levelplot(countQDS)->countPlot # improve plots
plot(countQDS)
plot(taxa.sf[1], add=TRUE)


##### TRY data ####
library(rtry)

path_to_data <- system.file("testdata", "data_TRY_15160.txt", package = "rtry")
path_to_data

data_TRY_15160 <- rtry_import(path_to_data)
input_path<-"C://Users//26485613//OneDrive - Stellenbosch University//Documents//Practice space"

tryFULL.df<-rtry_import(input_path,
                         separator = ",",
                         encoding = "UTF-8",
                         quote = "\"",
                         showOverview = TRUE)




##### Junks ####
taxa.sf %>%
  drop_na(species) %>%
  count(species, sort = TRUE) %>%
  nrow()

taxa.sf %>%
  filter(species=="Vachellia karroo") %>%
  count(dateIdentified) %>%
  nrow()


plot(countQDS)# plot grid cell
plot(taxa.sf[1], add=TRUE) #add occurence points



##### TRY data ####
uniqueName<-data.frame("uN"=unique(taxa.sf$species))

#TRYcat <- readxl::read_excel("TRYcat.xlsx", sheet = 1)
# create data of species with traits
speciesID<-inner_join(uniqueName,TryAccSpecies,
                        by=join_by("uN"=="AccSpeciesName")) %>%
  select(AccSpeciesID)


# get traits with more than 10000 observations
tTable<- table_traits %>%
  filter(ObsNum>10000)


speciesID<-as.numeric(speciesID$AccSpeciesID)
traitID<-as.numeric(tTable$TraitID)

# print in TRY input format
dput(speciesID)
dput(traitID)



input_path<-"C://Users//26485613//OneDrive - Stellenbosch University//Documents//Practice space//33312.txt"
# import data
try33576.df<-rtry_import(
  input="C:/Users/26485613/OneDrive - Stellenbosch University/Documents/Practice space/33576.txt",
  separator = "\t",
  encoding = "Latin-1",
  quote = "",
  showOverview = TRUE
)


try1<-try33576.df %>%
  # drop rows which contains no trait
  drop_na(TraitID) %>%
  # select species name, trait and trait value
  select(AccSpeciesName,TraitID,OrigValueStr) %>%
  #group by Species and trait
  group_by(AccSpeciesName,TraitID) %>%
  #choose the first trait value if there are multiples trait for a species
  summarise(across(OrigValueStr, first), .groups = "drop") %>%
  # reshape to wide format to have specie by trait dataframe
  pivot_wider(names_from = TraitID, values_from = OrigValueStr)






##### Junks ####
taxa.df %>%
  drop_na(species) %>%
  count(species, sort = TRUE) %>%
  filter(n>3)

summary(CATtrait.df)

CATtrait.df<- CATtrait.df %>%
  mutate(decade=ifelse(dateIdentified>as.Date("2024-02-28"),"First","Second"))

try33576.df %>%
  dplyr::group_by(AccSpeciesName, TraitID) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)

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





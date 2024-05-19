

library(dplyr)
library(readr)
library(tidyverse)
library(rgbif) # for occ_download
library(terra)
library(sf)
library(rtry) # for processing try data
library(rasterVis)
library(lubridate)
##### gbif data ####
taxa = 'Tracheophyta' # scientific name

gbif_download = occ_data(scientificName=taxa, # download data from gbif
                         country='ZA',
                         hasCoordinate=TRUE,
                         hasGeospatialIssue=FALSE,
                         limit = 2000)

taxa.df = as.data.frame(gbif_download$data) #extract data from the downloaded file


taxa.occ = taxa.df %>%
  dplyr::select(speciesKey,decimalLatitude,decimalLongitude,
                species,dateIdentified) %>% #select occurrence data
  filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
  mutate(dateIdentified = as.Date(dateIdentified)) # convert date to date format

taxa.sf<-st_as_sf(taxa.occ,coords = c("decimalLongitude", "decimalLatitude"),
                  crs = 4326) # convert long and lat point to geometry


#### Site by Species #####
# extract unique species name from GBIF occurrence data
uN<-unique(taxa.sf$species)
# Read RSA land area shapefile
rsa_country_sf = st_read("C:/Users/mukht/Documents/boundary_SA/boundary_south_africa_land_geo.shp")


# Create grid cells with extent of the data and layers for siteID and species
gridQDS = rast(rsa_country_sf,res=c(0.25,0.25), crs="EPSG:4326",nlyrs=length(uN)+1)
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
                         fun="max",
                         background = 0)
  # insert occurrence layer for each species to it assigned layer
  gridQDS[[n]] <- speciesQDS[]
})
# Make QDS Mask. Remember to mask the background to NA
rsa_mask = rasterize(rsa_country_sf, gridQDS, background=NA)

gridQDS_mask = mask(gridQDS, rsa_mask)

# create data frame of site by species
SitebySpecies <- as.data.frame(gridQDS_mask[])

##### speciebyxyt ####

taxa.sf$day <- yday(taxa.sf$dateIdentified)
taxa.sf <- taxa.sf %>% 
  mutate(period = case_when(
    day >= 1 & day <= 14 ~ 1,
    day >= 15 & day <= 28 ~ 2,
    day >= 29 & day <= 42 ~ 3,
    day >= 43 & day <= 56 ~ 4,
    day >= 57 & day <= 70 ~ 5,
    day >= 71 & day <= 84 ~ 6,
    day >= 85 & day <= 98 ~ 7,
    day >= 99 & day <= 112 ~ 8,
    day >= 113 & day <= 126 ~ 9,
    TRUE ~ NA_integer_  # Default case, if any value falls outside the specified ranges
  ))
# create an empty list for species by xyt
Speciesbyxyt <- list() 
for (t in sort(unique(taxa.sf$period))) {
  for(n in uN){
    # create raster of species
    speciesQDS = rasterize(dplyr::filter(taxa.sf, species==n, period==t),
                           gridQDS,
                           field=1,
                           fun="max",
                           background = 0)
    # insert occurrence layer for each species to it assigned layer
    gridQDS[[n]] <- speciesQDS[]
  }
  
  # mask the site
  gridQDS_mask.t = mask(gridQDS, rsa_mask)
  
  # create data frame of site by species
  Speciesbyxyt[[paste0("P",t)]] <- as.data.frame(t(gridQDS_mask.t[]))
  
}

##### Specie by trait  #####
# Make the rows
uniqueName<-data.frame("uN"=uN)

# get TRY species ID for gbif species data
speciesID<-inner_join(uniqueName,TryAccSpecies,
                      by=join_by("uN"=="AccSpeciesName")) %>%
  select(AccSpeciesID)


# get traits with more than 10000 observations
tTable<- table_traits %>%
  filter(ObsNum>10000)
# create specieID and traitID for request from TRY database
speciesID<-as.numeric(speciesID$AccSpeciesID)
traitID<-as.numeric(tTable$TraitID)

# print in TRY input format
dput(speciesID)
dput(traitID)



input_path<-"C:/Users/mukht/Downloads/33576.txt"
  #"C:/Users/26485613/OneDrive - Stellenbosch University/Documents/Practice space/33576.txt"
# 
try33576<-rtry_import(
  input=input_path,
  separator = "\t",
  encoding = "Latin-1",
  quote = "",
  showOverview = TRUE
)


SpeciesbyTrait<-try33576 %>%
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
{na.df<-as.data.frame(matrix(NA,nrow = length(setdiff(uN,rownames(SpeciesbyTrait))),
                            ncol = ncol(SpeciesbyTrait)))
row.names(na.df)<-setdiff(uN,rownames(SpeciesbyTrait))
names(na.df)<-names(SpeciesbyTrait)
SpeciesbyTrait<-rbind(SpeciesbyTrait,na.df)}
#### Species by xyt ####
# ID = StringJoin('Long','Lat','Time')
# create longitude, latitude and time vectors
lon <- taxa.occ$decimalLongitude # x
lat <- taxa.occ$decimalLatitude # y
Time <- taxa.occ$dateIdentified # t
# concatenate longitude, latitude and time
xyt <- paste(lon, lat, Time, sep = ",")
# bind xyt to occurrence dataframe
taxa.occ$xyt<-xyt
# create value for presence count
taxa.occ$count <- 1

Speciesbyxyt <- taxa.occ %>%
  select(species,count,xyt) %>%
  #group by Species and xyt
  group_by(species,xyt) %>%
  # count and add species occurrence in each xyt
  summarise(across(count, sum), .groups = "drop") %>%
  # reshape to wide format to have specie by xyt dataframe
  pivot_wider(names_from = xyt, values_from = count)



#####
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


uN[7]
plot(gridQDS[[uN[7]]])
plot(filter(taxa.sf, species==uN[7]),add=TRUE)

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
  filter(n>10)




plot(countQDS)# plot grid cell
plot(taxa.sf[1], add=TRUE) #add occurrence points







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




gridQDS
plot(gridQDS[[uN[90:93]]])


r <- rast(ext(taxa.sf),nrows=4, ncols=3, nlyrs=3)
# Assign unique IDs to each cell
names(r)<- c("ID","sp1","sp2")
r$ID <- 1:ncell(r)
r$sp1<-data.frame(rfield)
r$sp2<-rbinom(12,1,0.3)
names(gridQDS)[1:10]

r[]
# Print the raster to see the unique IDs
print(r)
plot(r$sp1)
plot(sp1, add=T)
# Get coordinates of each cell
coords <- xyFromCell(r, 1:12)

# Print the coordinates
print(coords)
data.frame(r)


rfield = rasterize(dplyr::filter(taxa.sf, species==uN[67]),
                     r,
                     field=1,
                     fun="max",
                     background = 0)

data.frame(rfield)
plot(rfield)
plot(dplyr::filter(taxa.sf, species==uN[67]), add=T)

sp1<-taxa.sf$geometry[1:12]
sp2<-taxa.sf$geometry[101:112]

taxa.sf %>%
  group_by(species) %>%
  filter(species=="Clerodendrum ternatum") #%>%
  #count(species)
taxa.sf[,2]
dplyr::filter(taxa.sf, species==uN[67])
sum(speciesQDS[])


chelsaA18 <- terra::rast('CHELSA_swb_2018_V.2.1.tif')
plot(chelsaA18)
chelsa.SA<-crop(chelsaA18,ext(taxa.sf))
plot(chelsa.SA)
chelsa.SA[]
envSA<-rasterize(chelsa.SA,
                 gridQDS,
                 field=1,
                 fun=mean,
                 background = 0)


juli<-data.frame("origDate"=taxa.occ$dateIdentified,
                 "JuliDate"=julian(taxa.occ$dateIdentified))
lubri<-data.frame("origDate"=taxa.occ$dateIdentified,
                  "lubriDate"=yday(taxa.occ$dateIdentified))



t=7
system.time(for(n in uN){
  # create raster of species
  speciesQDS = rasterize(dplyr::filter(taxa.sf, species==n, period==7),
                         gridQDS,
                         field=1,
                         fun="max",
                         background = 0)
  # insert occurrence layer for each species to it assigned layer
  gridQDS[[n]] <- speciesQDS[]
}
)

# mask the site
gridQDS_mask.t = mask(gridQDS, rsa_mask)

# create data frame of site by species
Speciesbyxyt[[paste0("P",t)]] <- as.data.frame(t(gridQDS_mask.t[]))


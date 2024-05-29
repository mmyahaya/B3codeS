

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
taxa = 'Acacia' # scientific name

gbif_download = occ_data(scientificName=taxa, # download data from gbif
                         country='ZA',
                         hasCoordinate=TRUE,
                         hasGeospatialIssue=FALSE,
                         limit = 17726)

taxa.df = as.data.frame(gbif_download$data) #extract data from the downloaded file


taxa.occ = taxa.df %>%
  dplyr::select(decimalLatitude,decimalLongitude,
                species,dateIdentified) %>% #select occurrence data
  filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
  mutate(dateIdentified = as.Date(dateIdentified)) # convert date to date format
taxa.sf<-st_as_sf(taxa.occ,coords = c("decimalLongitude", "decimalLatitude"),
                  crs = 4326) # convert long and lat point to geometry


#### Site by Species #####
# extract unique species name from GBIF occurrence data
uN<-sort(unique(taxa.sf$species))
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
                         fun="sum",
                         background = 0)
  # insert occurrence layer for each species to it assigned layer
  gridQDS[[n]] <- speciesQDS[]
})
# Mask the grid cell to country shape file

gridQDS = mask(gridQDS, rsa_country_sf)
#plot(gridQDS[[21:40]])
#lines(rsa_country_sf['Land'])

# create data frame of site by species
sbs <- as.data.frame(gridQDS)
sbs<-drop_na(sbs,siteID)
sbsM<-as.matrix(sbs[,-1])
colnames(sbsM)<-NULL
##### speciebyxyt ####

#taxa.sf$day <- yday(taxa.sf$dateIdentified)
taxa.sf$year<-year(taxa.sf$dateIdentified)
taxa.sf <- taxa.sf %>% 
  mutate(period = case_when(
    year == 2024 ~ 1,
    year == 2023  ~ 2,
    year == 2022  ~ 3,
    year == 2021 ~ 4,
    year == 2020  ~ 5,
    year == 2019  ~ 6,
    year == 2018  ~ 7,
    year <= 2017  ~ 8,
    TRUE ~ NA_integer_  # Default case, if any value falls outside the specified ranges
  ))
# create an empty list for species by xyt
Speciesbyxyt <- list() 
system.time(for (t in sort(unique(taxa.sf$period))) {
  for(n in uN){
    # create raster of species
    speciesQDS = rasterize(dplyr::filter(taxa.sf, species==n, period==t),
                           gridQDS,
                           field=1,
                           fun="sum",
                           background = 0)
    # insert occurrence layer for each species to it assigned layer
    gridQDS[[n]] <- speciesQDS[]
  }
  
  # mask the site
  gridQDS.t = mask(gridQDS, rsa_country_sf)
  
  sbs.t <- as.data.frame(gridQDS.t)
  sbs.t<-drop_na(sbs.t,siteID)
  sbsM.t<-as.matrix(sbs[,-1])
  colnames(sbsM.t)<-NULL
  
  # create data frame of site by species
  
  
  Speciesbyxyt[[paste0("T",t)]] <- sbsM.t
  
})


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



input_path<-"33852.txt"
  #"C:/Users/mukht/Downloads/33576.txt"
  #"C:/Users/26485613/OneDrive - Stellenbosch University/Documents/Practice space/33576.txt"

try33852<-rtry_import(
  input=input_path,
  separator = "\t",
  encoding = "Latin-1",
  quote = "",
  showOverview = TRUE
)


SpeciesbyTrait<-try33852 %>%
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
names(na.df)<-names(SpeciesbyTrait) # column names
SpeciesbyTrait<-rbind(SpeciesbyTrait,na.df)
SpeciesbyTrait<-SpeciesbyTrait[order(rownames(SpeciesbyTrait)),]
sbtM<-as.matrix(SpeciesbyTrait)
rownames(sbtM)<-NULL
colnames(sbtM)<-NULL}



traitname<-try33852 %>%
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
#### Site by Environment ####
path = "C:/Users/mukht/Documents"

# Download the WorldClim Bioclimatic variables for the world at a 10 arc-minute resolution
bio_10m = geodata::worldclim_global(var='bio',
                                    res=10, path=path,
                                    version="2.1") # Set your own path directory

# Define 'extent' of boundary
# South Africa
rsa_ext = raster::extent(rsa_country_sf)

# Crop Bioclimatic variables to extent of South African boundary
rsa_bio_10m = crop(bio_10m, rsa_ext)



# Transfer values from worldclim raster data to QDS
bioQDS<-resample(rsa_bio_10m,gridQDS) # bilinear interpolation 
bioQDS[["siteID"]]<-1:ncell(bioQDS)
# mask bioQDS to rsa land
bioQDS<-mask(bioQDS,rsa_country_sf)

# extract site by environment from the bioQDS layers
{sitebyEnv <- as.data.frame(bioQDS[])
sitebyEnv <- drop_na(sitebyEnv,siteID)
sbeM<-as.matrix(dplyr::select(sitebyEnv, !siteID)) 
colnames(sbeM)<-NULL}








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
coords <- raster::xyFromCell(r,c(1,4,5,11))






chelsaA18 <- terra::rast('CHELSA_swb_2018_V.2.1.tif')
plot(chelsaA18)
chelsa.SA<-crop(chelsaA18,ext(taxa.sf))
plot(chelsa.SA)
chelsa.SA[]
envSA<-rasterize(chelsa.SA, #ERROR
                 gridQDS,
                 field=1,
                 fun=mean,
                 background = 0)



df <- apply(taxa.df,2,as.character)
write.table(df,"taxa(Acacia).csv",row.names = F)
length(base::setdiff(unique(taxa.df$species),uN))


plot(gridQDS[[uN[24]]])
plot(dplyr::filter(taxa.sf, species==uN[24]),add=T)
lines(rsa_country_sf['Land'])


path = "C:/Users/mukht/Documents" #path for worldclim

precdata <- terra::rast('prec_2021-2040.tif')
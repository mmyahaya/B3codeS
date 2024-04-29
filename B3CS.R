

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



library(sf)

# Define the bounding box coordinates
bbox <- st_bbox(rsa_ext, crs = 4326)

# Create a grid with 1 km x 1 km cells
grid <- st_make_grid(bbox, cellsize = c(.020, .020),n=20, what = "polygons")

# Plot the grid
plot(grid)


lonlat<-cbind(taxa.occ$decimalLongitude,taxa.occ$decimalLatitude)

pts<-vect(lonlat)
pts <- vect(lonlat, crs="+proj=longlat +datum=WGS84")
rsa_ext = extent(16, 33, -35, -22)


plot(st_make_grid(what = "polygons"), axes = TRUE)
plot(st_make_grid(what = "corners"), add = TRUE, col = 'green', pch=3)
sfc = st_sfc(st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,0)))))
plot(st_make_grid(sfc, cellsize = .1, square = FALSE))
points(pts, add = TRUE)
st_make_grid(what = "centers")
grid

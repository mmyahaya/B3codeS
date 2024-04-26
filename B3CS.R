

library(dplyr)
library(readr)
library(tidyverse)
library(rgbif) # for occ_download
library(terra)
library(sf)

##### gbif data ####
taxa = 'Tracheophyta' # scientific name

gbif_download = occ_data(scientificName=taxa,
                         country='ZA',
                         hasCoordinate=TRUE,
                         hasGeospatialIssue=FALSE,
                         limit = 2000)

taxa.df = as.data.frame(gbif_download$data)


taxa.occ = taxa.df %>%
  dplyr::select(key,decimalLatitude,decimalLongitude,occurrenceStatus,genus,
                species,genericName,dateIdentified) %>%
  filter_all(all_vars(!is.na(.))) %>%
  mutate(dateIdentified = as.Date(dateIdentified))

taxa.sf<-st_as_sf(taxa.occ,coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Define a grid over spatial extend
grid<- st_make_grid(
  taxa.sf,
  square = TRUE,
  cellsize = c(0.2, 0.2)
)
# create intersect
intersected <- st_intersection(grid, taxa.sf)

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





rsa_ext = extent(16, 33, -35, -22)

# head(lepidop.df)

lepidop.sf = st_as_sf(lepidop.df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
plot(lepidop.sf['stateProvince'])

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

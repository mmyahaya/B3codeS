

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
                         limit = 2000)

taxa.df = as.data.frame(gbif_download$data) #extract data from the downloaded file


taxa.occ = taxa.df %>%
  dplyr::select(decimalLatitude,decimalLongitude,
                species,dateIdentified,year,month,day) %>% #select occurrence data
  filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
  mutate(dateIdentified = as.Date(dateIdentified)) # convert date to date format
taxa.sf<-st_as_sf(taxa.occ,coords = c("decimalLongitude", "decimalLatitude"),
                  crs = 4326) # convert long and lat point to geometry


#### Site by Species #####
# extract unique species name from GBIF occurrence data
uN<-sort(unique(taxa.sf$species))
# Read RSA land area shapefile
rsa_country_sf = st_read("C:/Users/mukht/Documents/boundary_SA/boundary_south_africa_land_geo.shp")

plot(rsa_country_sf)
# Create grid cells with extent of the data and layers for siteID and species
gridQDS = rast(rsa_country_sf,res=c(0.25,0.25), crs="EPSG:4326",nlyrs=length(uN)+1)
# specify name for each layers of site ID and individual species
names(gridQDS)<-c("siteID",uN)
# Assign ID for each cell
gridQDS[["siteID"]]<-1:ncell(gridQDS)
# create layer for species occurrence in each cell

system.time(species_stack <- lapply(uN, function(n) {
  # Create raster of species occurrences
  speciesQDS <- rasterize(dplyr::filter(taxa.sf$taxa, species == n), 
                          gridQDS, field = 1, fun = "sum", background = 0)
  names(speciesQDS)<-n
  return(speciesQDS)
}))
#convert each species to layers of raster
species_stack <- rast(species_stack)
# marge with site ID 
gridQDS <- c(gridQDS[["siteID"]], species_stack)

# Mask the grid cell to country shape file

gridQDS = mask(gridQDS, rsa_country_sf)
library(RColorBrewer)
plot(gridQDS[[5:10]],col=brewer.pal(5,"OrRd"))
lines(rsa_country_sf['Land'])

# create data frame of site by species
sbs <- as.data.frame(gridQDS)
sbs<-drop_na(sbs,siteID)
sbsM<-as.matrix(sbs[,-1])
colnames(sbsM)<-NULL
##### speciebyxyt ####

#taxa.sf$day <- yday(taxa.sf$dateIdentified)
#taxa.sf$year<-year(taxa.sf$dateIdentified)
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






remotes::install_github('matildabrown/rWCVPdata')
rWCVP::wcvp_checklist(taxon = "Acacia", taxon_rank = c("genus"))->A

rWCVP::wcvp_occ_mat(taxon = "Acacia", taxon_rank = c("genus"))->B
dplyr::filter(A,taxon_name==uN)->C

pivot_wider(names_from = TraitID, values_from = OrigValueStr) 
tidyr::pivot_wider(data,names_from ="B", values_froml = "B")

# Add new columns "simple" and "complex"
data <- data %>%
  mutate(simple = ifelse(B == "simple", 1, NA),
         complex = ifelse(B == "complex", 1, NA))
# Get the unique non-NA values from column B
unique_values <- na.omit(unique(data$B))

# Create the new columns dynamically
data_with_new_cols <- data %>%
  mutate(across(all_of(unique_values), 
                ~ ifelse(B == cur_column(), 1, NA), 
                .names = "{col}"))


data_with_new_cols <- data %>%
  mutate(value = 1) %>%  # Create a temporary column with value 1 for pivoting
  pivot_wider(names_from = B, values_from = value, values_fill = NA)





To convert categorical values in columns to new columns (often referred to as one-hot encoding), you can use the dummies package or the fastDummies package, which provide convenient functions for this purpose. Here’s how you can do it using both packages:
  
  Using fastDummies Package
Install and load the fastDummies package.
Use the dummy_cols function to convert categorical columns to new columns.
r
Copy code
# Install and load the fastDummies package
install.packages("fastDummies")
library(fastDummies)

# Create a sample dataframe
data <- data.frame("A" = 1:5, 
                   "B" = c("simple", NA, "complex", "simple", NA), 
                   "C" = factor(c("low", "high", "medium", "medium", "high")))

# Convert categorical columns to new columns
data_with_dummies <- dummy_cols(data, select_columns = c("B", "C"), remove_first_dummy = FALSE, remove_selected_columns = TRUE)

# View the resulting dataframe
print(data_with_dummies)
Using dummies Package
Install and load the dummies package.
Use the dummy.data.frame function to convert categorical columns to new columns.
r
Copy code
# Install and load the dummies package
install.packages("caret")
library(caret)

# Create a sample dataframe
data <- data.frame("A" = 1:5, 
                   "B" = c("simple", NA, "complex", "simple", NA), 
                   "C" = factor(c("low", "high", "medium", "medium", "high")))

# Convert categorical columns to new columns
data_with_dummies <- dummy.data.frame(data, names = c("B", "C"), sep = "_")

# View the resulting dataframe
encoded_data <- data.frame(data, dummyVars(data$B))
?dummyVars


data <- data.frame(
  A = 1:5, 
  B = c("simple", NA, "complex", "simple", NA), 
  C = factor(c("low", "high", "medium", "medium", NA)),
  D = c("cat", NA, "cat", "dog", "fish"),
  E = c(1, 2, 3, NA, 5)  # This is a numeric column
)

replace(data, is.na(data), "NA")
# Load necessary library
library(dplyr)

# Create a sample dataframe
data <- data.frame("A" = 1:5, 
                   "B" = c("simple", NA, "complex", "simple", NA), 
                   "C" = factor(c("low", "high", "medium", "medium", "high")))

# Replace NAs in categorical columns with a placeholder (e.g., "NA")
data <- lapply(data, function(x) if(is.character(x)) replace(x, is.na(x), "NA") else x)

# Identify categorical columns
categorical_cols <- sapply(data, is.factor) | sapply(data, is.character)
categorical_cols <- names(categorical_cols[categorical_cols])

# Initialize a list to store dummy variables
dummy_list <- list()

# Loop through each categorical column and create dummy variables
for (col in categorical_cols) {
  # Create dummy variables for the column
  dummy_matrix <- model.matrix(as.formula(paste("~", col, "- 1")), data)
  # Add the dummy variables to the list
  dummy_list[[col]] <- dummy_matrix
}

# Combine the original dataframe with the dummy variables
data_with_dummies <- cbind(data, do.call(cbind, dummy_list))

# Optionally remove the original categorical columns
data_with_dummies <- data_with_dummies %>%
  select(-all_of(categorical_cols))

# View the resulting dataframe
print(data_with_dummies)




library(dplyr)
library(readr)
library(tidyverse)
library(rgbif) # for occ_download
library(terra)
library(sf)
library(rtry) # for processing try data
library(rasterVis)


##### gbif data ####
taxa = 'Fabaceae' # scientific name
#Tracheophyta Fabaceae
gbif_download = occ_data(scientificName=taxa, # download data from gbif
                         country='ZA',
                         hasCoordinate=TRUE,
                         hasGeospatialIssue=FALSE,
                         limit = 70000)

#taxa.df1 = as.data.frame(gbif_download[["Acacia"]]$data) %>%  dplyr::select(decimalLatitude,decimalLongitude,
  #                                                                 species,coordinateUncertaintyInMeters,dateIdentified,year,month,day)
#taxa.df2 = as.data.frame(gbif_download[["Vachellia"]]$data) |> dplyr::select(decimalLatitude,decimalLongitude,
    #                                                                 species,coordinateUncertaintyInMeters,dateIdentified,year,month,day)
#extract data from the downloaded file

#taxa.df=rbind(taxa.df1,taxa.df2)
taxa.occ =gbif_download$data %>%
  dplyr::select(decimalLatitude,decimalLongitude,#"coordinateUncertaintyInMeters"
                species,coordinateUncertaintyInMeters,dateIdentified,year,month,day) %>% #select occurrence data
  filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
  mutate(dateIdentified = as.Date(dateIdentified)) # convert date to date format
taxa.sf<-st_as_sf(taxa.occ,coords = c("decimalLongitude", "decimalLatitude"),
                  crs = 4326) # convert long and lat point to geometry


#### Site by Species #####
# extract unique species name from GBIF occurrence data
specie_list<-sort(unique(taxa.sf$species))
# Read RSA land area shapefile
rsa_country_sf = st_read("C:/Users/mukht/Documents/boundary_SA/boundary_south_africa_land_geo.shp")

plot(rsa_country_sf)
# Create grid cells with extent of the data and layers for siteID and species
gridQDS = rast(rsa_country_sf,res=c(0.25,0.25), crs="EPSG:4326",nlyrs=length(specie_list)+1)
# specify name for each layers of site ID and individual species
names(gridQDS)<-c("siteID",specie_list)
# Assign ID for each cell
gridQDS[["siteID"]]<-1:ncell(gridQDS)
# create layer for species occurrence in each cell

system.time(species_stack <- lapply(specie_list, function(n) {
  # Create raster of species occurrences
  speciesQDS <- rasterize(dplyr::filter(taxa.sf, species == n), 
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
plot(gridQDS[[20:29]],col=brewer.pal(5,"OrRd"))
lines(rsa_country_sf)

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
  for(n in specie_list){
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
uniqueName<-data.frame("specie_list"=specie_list)

# get TRY species ID for gbif species data
speciesID<-inner_join(uniqueName,TryAccSpecies,
                      by=join_by("specie_list"=="AccSpeciesName")) %>%
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



input_path<-"TRY_Vascular.txt"
#"C:/Users/mukht/Downloads/33576.txt"
#"C:/Users/26485613/OneDrive - Stellenbosch University/Documents/Practice space/33576.txt"

try_Vascular<-rtry::rtry_import(
  input=input_path,
  separator = "\t",
  encoding = "Latin-1",
  quote = "",
  showOverview = TRUE
)


SpeciesbyTrait<-try_Vascular %>%
  # drop rows which contains no trait
  drop_na(TraitID) %>%
  # select species name, trait and trait value
  select(AccSpeciesName,TraitID,TraitName,OrigValueStr,Comment) %>%
  #group by Species and trait
  group_by(AccSpeciesName,TraitID) %>%
  #choose the first trait value if there are multiples trait for a species
  summarise(across(OrigValueStr, first), .groups = "drop") %>%
  # reshape to wide format to have specie by trait dataframe
  pivot_wider(names_from = TraitID, values_from = OrigValueStr) %>% 
  # select species that are only present in gbif data
  filter(AccSpeciesName %in% species_list) %>% 
  # convert species names to row names
  column_to_rownames(var = "AccSpeciesName") 
# add the rows of the remaining species without traits from TRY
{na.df<-as.data.frame(matrix(NA,nrow = length(setdiff(species_list,rownames(SpeciesbyTrait))),
                             ncol = ncol(SpeciesbyTrait)))
  row.names(na.df)<-setdiff(species_list,rownames(SpeciesbyTrait))
  names(na.df)<-names(SpeciesbyTrait) # column names
  SpeciesbyTrait<-rbind(SpeciesbyTrait,na.df)
  SpeciesbyTrait<-SpeciesbyTrait[order(rownames(SpeciesbyTrait)),]}

  {sbtM<-as.matrix(SpeciesbyTrait)
  rownames(sbtM)<-NULL
  colnames(sbtM)<-NULL}



traitname<-try_all %>%
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
                                    res=2.5, path=path,
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




# World checklist for vascular plants

remotes::install_github('matildabrown/rWCVPdata')
library(rWCVP)
library(rWCVPdata)

# Download WCVP for acacia within South African area
AcaciaWCVP <- rWCVP::wcvp_checklist() %>% 
  filter(area_code_l3 %in% get_wgsrpd3_codes("South Africa"))
  

# List of my species
my_Acacia_list<-data.frame( taxon = c("Acacia acinacea", "Acacia adunca", "Acacia baileyana", "Acacia binervata", 
  "Acacia crassiuscula", "Acacia cultriformis", "Acacia cyclops", 
  "Acacia dealbata", "Acacia decurrens", "Acacia elata", "Acacia falciformis", 
  "Acacia implexa", "Acacia longifolia", "Acacia mearnsii", "Acacia melanoxylon", 
  "Acacia paradoxa", "Acacia piligera", "Acacia podalyriifolia", 
  "Acacia provincialis", "Acacia pubescens", "Acacia pycnantha", 
  "Acacia retinodes", "Acacia saligna", "Acacia schinoides", "Acacia stricta", 
  "Acacia ulicifolia", "Acacia viscidula"))

# Check list for my species
SA_native_list <- AcaciaWCVP %>% 
  filter((accepted_name %in% uSpecies) & occurrence_type=="native") %>% 
  select(accepted_name,occurrence_type,area,region)


# New dataframe with introduction status

my_Acacia_list_status<-
  my_Acacia_list%>%
  mutate(introduction_status = ifelse(taxon%in%SA_native_list$accepted_name,
                                      "native","introduced"))





points_in_sea <- taxa.sf[!st_intersects(taxa.sf, rsa_country_sf, sparse = FALSE), ]


lines(rsa_country_sf['Land'])


specie_list= sbs$species.name

# Confidence score for pre cautious country level impact
confidence_M = eicat_data %>% 
  mutate(impact_category=substr(impact_category,1,2)) %>% 
  mutate(impact_category=factor(impact_category,levels=c("DD","MC","MN","MO","MR","MV"))) %>% 
  filter(impact_category!="NA") %>% 
  filter(impact_country=="South Africa") %>% 
  select(scientific_name,
         impact_category,confidence_rating) %>% 
  mutate(confidence_rating = case_when(
    confidence_rating %in% c("low","Low") ~ 1,
    confidence_rating %in% c("medium", "Medium") ~ 2,
    confidence_rating %in% c("high", "High") ~ 3,
    TRUE ~ 0  # Default case, if any value falls outside the specified ranges
  )) %>% 
  group_by(scientific_name,impact_category) %>%
  #choose the first trait value if there are multiples trait for a species
  summarise(across(confidence_rating, max), .groups = "drop") %>% 
  #reshape to wide format to have specie by trait dataframe
  pivot_wider(names_from = impact_category, values_from = confidence_rating) %>% 
  #select species that are only present in gbif data
  filter(scientific_name %in% specie_list) %>% 
  #convert species names to row names
  column_to_rownames(var = "scientific_name") 
  #add the rows of the remaining species without traits from TRY
  {na.df<-as.data.frame(matrix(NA,nrow = length(setdiff(specie_list,rownames(confidence_M))),
                               ncol = ncol(confidence_M)))
    row.names(na.df)<-setdiff(specie_list,rownames(confidence_M))
    names(na.df)<-names(confidence_M) # column names
    confidence_M<-rbind(confidence_M,na.df)
    confidence_M<-confidence_M[order(rownames(confidence_M)),]
    cbiM<-as.matrix(confidence_M)
    #rownames(sbtM)<-NULL
    #colnames(sbtM)<-NULL
    }


category_M = eicat_data %>% 
  mutate(impact_category=substr(impact_category,1,2)) %>% 
  filter(impact_category %in% c("MC","MN","MO","MR","MV")) %>% 
  #filter(impact_country=="South Africa") %>% 
  select(scientific_name,
         impact_category) %>% 
  mutate(category_value = case_when(
    impact_category == "MC" ~ 0,
    impact_category == "MN" ~1,
    impact_category == "MO" ~2,
    impact_category == "MR" ~3,
    impact_category == "MV" ~4,
    TRUE ~ 0  # Default case, if any value falls outside the specified ranges
  )) %>% 
  #mutate(impact_category=factor(impact_category,levels=c("MC","MN","MO","MR","MV"))) %>% 
  group_by(scientific_name,impact_category) %>%
  #choose the first trait value if there are multiples trait for a species
  summarise(across(category_value, first), .groups = "drop") %>% 
  # reshape to wide format to have specie by trait dataframe
  pivot_wider(names_from = impact_category, values_from = category_value) %>% 
  # select species that are only present in gbif data
  filter(scientific_name %in% specie_list) %>% 
  #convert species names to row names
  column_to_rownames(var = "scientific_name")
  
frequency_M
{na.df<-as.data.frame(matrix(NA,nrow = length(setdiff(specie_list,rownames(category_M))),
                               ncol = ncol(category_M)))
    row.names(na.df)<-setdiff(specie_list,rownames(category_M))
    names(na.df)<-names(category_M) # column names
    category_M<-rbind(category_M,na.df)
    category_M<-category_M[order(rownames(category_M)),]
    category_M<-as.matrix(category_M)
    #rownames(sbtM)<-NULL
    #colnames(sbtM)<-NULL
    }


frequency_M<-eicat_data %>% 
  mutate(impact_category=substr(impact_category,1,2)) %>% 
  filter(impact_category %in% c("MC","MN","MO","MR","MV")) %>% 
  #filter(impact_country=="South Africa") %>% 
  select(scientific_name,
         impact_category) %>% 
  mutate(category_value = case_when(
    impact_category == "MC" ~ 0,
    impact_category == "MN" ~1,
    impact_category == "MO" ~2,
    impact_category == "MR" ~3,
    impact_category == "MV" ~4,
    TRUE ~ 0  # Default case, if any value falls outside the specified ranges
  )) %>% 
  #mutate(impact_category=factor(impact_category,levels=c("MC","MN","MO","MR","MV"))) %>% 
  group_by(scientific_name,impact_category) %>%
  #choose the first trait value if there are multiples trait for a species
  summarise(across(category_value, length), .groups = "drop") %>% 
  # reshape to wide format to have specie by trait dataframe
  pivot_wider(names_from = impact_category, values_from = category_value) %>% 
  # select species that are only present in gbif data
  filter(scientific_name %in% specie_list) %>% 
  #convert species names to row names
  column_to_rownames(var = "scientific_name") %>% 
  rowwise() %>%
  mutate(Total = sum(c_across(everything()), na.rm = TRUE))


category_Mm<-apply(category_M,1, function(x) max(x,na.rm = T))
sbs$sbs%*%diag(category_Mm)->impact_score

impact_M = eicat_data %>% 
  mutate(impact_category=substr(impact_category,1,2)) %>% 
  filter(impact_category!="NA") %>% 
  select(scientific_name,
         impact_category,confidence_rating) %>% 
  mutate(confidence_rating = case_when(
    confidence_rating %in% c("low","Low") ~ 1,
    confidence_rating %in% c("medium", "Medium") ~ 2,
    confidence_rating %in% c("high", "High") ~ 3,
    TRUE ~ 0  # Default case, if any value falls outside the specified ranges
  )) %>% 
  group_by(scientific_name,impact_category) %>%
  #choose the first trait value if there are multiples trait for a species
  summarise(across(confidence_rating, length), .groups = "drop") %>% 
  # reshape to wide format to have specie by trait dataframe
  pivot_wider(names_from = impact_category, values_from = confidence_rating) %>% 
  # select species that are only present in gbif data
  filter(scientific_name %in% specie_list)
  



data.frame("name"=specie_list, "introduction_status"=sbt$sbt[,ncol(sbt$sbt)])


intro.sf<-taxon_cube$data %>% 
  left_join(taxa_list_status,
            by = c("scientificName" = "taxon")) %>% 
  as.data.frame()





gridQDS = rast(rsa_country_sf,res=c(0.25,0.25), crs="EPSG:4326")
# specify name for each layers of site ID and individual species

# create layer for species occurrence in each cell

system.time(species_stack <- lapply(c("introduced","native"), function(n) {
  # Create raster of species occurrences
  speciesQDS <- rasterize(dplyr::filter(intro.sf, introduction_status == n), 
                          gridQDS, field = 1, fun = "sum", background = NA)
  names(speciesQDS)<-n
  return(speciesQDS)
}))
#convert each species to layers of raster
species_stack <- rast(species_stack)
species_stack$intro_native = sqrt(species_stack$introduced/species_stack$native)

plot(species_stack)
lines(rsa_country_sf['Land'])
A=as.numeric(A)

gridQDS$lyr.1<-B
B=as.data.frame(B[])
B[A,]<-rowSums(sbs$sbs.binary)
plot(gridQDS[[2]])


intro.sf <- intro.sf %>% 
  filter(coordinateUncertaintyInMeters<=250)

intro.sf$coordinateUncertaintyInMeters = (intro.sf$coordinateUncertaintyInMeters/(res*1000))^2
uncertaintyDS <- rasterize(intro.sf, 
                        gridQDS, field = "coordinateUncertaintyInMeters", fun = mean, background = NA)
plot(uncertaintyDS)




northEU<-st_read("C:/Users/mukht/Downloads/world-administrative-boundaries/world-administrative-boundaries.shp")
ZAsf<-filter(northEU,name=="South Africa") %>% select(name,geometry)
plot(ZAsf)


# Save the sf object to a .rds file
saveRDS(northEU, file = "C:/Users/mukht/Documents/B3codeS/countries_shapefile.rds")
countries_sf<-readRDS("C:/Users/mukht/Documents/B3codeS/countries_shapefile.rds")
SA.sf<-filter(countries_sf,name=="South Africa") %>% select(name,geometry)
plot(SA.sf)

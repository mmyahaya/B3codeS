taxa.data <- function(taxa,limit){
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
#taxa.sf <- taxa.data('Acacia',500)
# Save taxa.sf to drive to avoid redownloading same taxa with same limit next time
#  write.csv(taxa.sf,"taxa.sf.csv",row.names = FALSE)

sbsMatrix <- function(taxa.sf,country.shp){
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
  # Make QDS Mask with NA as background 
  rsa_mask = rasterize(country.shp, gridQDS, background=NA)
  gridQDS = mask(gridQDS, rsa_mask)
  
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
sbs<-sbsMatrix(taxa.sf = taxa.sf, country.shp = rsa_country_sf)

#' @title Create site-by-environment matrix
#' 
#' @description Create site-by-environment (sbe) matrix over the given region.
#' The rows and columns of the matrix represents the sites and environmental
#' variables respectively. 
#' @param rastfile Raster. should be path for storing the worldclim data 
#' (see geodata::worldclim_global) or any other raster data of environmental 
#' variable
#' @param country.sf sf object. The shapefile of the region of study
#' @param res Numeric. The resolution of grid cells to be used. Default is 0.25
#' @param siteID Vector. The cellids of site with occurrences. An output of `sbs`

sbeFun <- function(rastfile,country.sf,res=0.25,siteID){
  
  # read the rast data if path is given
  if("character" %in% class(rastfile)){
    # Download the WorldClim Bioclimatic variables for the world at a 10 arc-minute resolution
    env <- geodata::worldclim_global(var='bio',
                                     res=2.5, path=rastfile,
                                     version="2.1")
    names(env) <- c('AnnTemp','DiurRange','Isotherm','TempSeas',
                    
                    'MaxTemp','MinTemp','TempRange',
                    
                    'MeanTWQ','MeanTDQ','MeanTWaQ','MeanTCQ','AnnPrec',
                    
                    'PrecWetM','PrecDrM','PrecSeas','PrecWetQ',
                    
                    'PrecDrQ','PrecWaQ','PrecCQ')
  } else if("SpatRaster" %in% class(rastfile)){
    env<-rastfile
  } else { # stop and report if rastfile is not a file path or SpatRaster
    cli::cli_abort(c("{.var rastfile} is not a file path or SpatRaster"))
  }
  
  # Crop environmental variables to extent of of the country's boundary
  env <- terra::crop(env, terra::ext(country.sf))
  # Convert the cropped raster to points, removing NA values
  env_points <- as.data.frame(env, xy = TRUE, na.rm = TRUE)
  
  # Define grid cell for sites and environmental variables
  gridQDS <- terra::rast(country.sf,res=c(res,res), crs="EPSG:4326")
  envQDS <- terra::rast(country.sf,res=c(res,res), crs="EPSG:4326")
  # create raster for environmental data
  for (v in names(env_points[,-c(1,2)])) {
    idw_model <- gstat::gstat(id = v ,formula = as.formula(paste(v,"~1")), data = env_points,
                              locations = ~x+y,nmax = 7, set = list(idp = 2))
    interpolated_raster <- terra::interpolate(gridQDS, idw_model)
    envQDS <- c(envQDS,interpolated_raster[[1]])
  }
  
  #create grid for region
  grid <- country.sf %>%
    sf::st_make_grid(cellsize = c(res,res),
                     offset = c(sf::st_bbox(country.sf)$xmin,
                                sf::st_bbox(country.sf)$ymin)) %>%
    sf::st_sf() %>%
    dplyr::mutate(cellid = dplyr::row_number())
  
  # create site by environment
  sbeM <- as.data.frame(envQDS, xy = TRUE, na.rm = TRUE) %>% 
    sf::st_as_sf(coords = c("x", "y"),crs = 4326) %>% 
    sf::st_join(grid) %>% 
    dplyr::filter(cellid %in% as.integer(siteID)) %>% 
    as.data.frame() %>% 
    dplyr::arrange(cellid) %>% 
    dplyr::select(-any_of(c("cellid","geometry"))) %>% 
    as.matrix()
  
  
  # collect environmental variable name
  variable.name<-colnames(sbeM)
  colnames(sbeM)<-NULL
  return(list("sbe"=sbeM,"variable.name"=variable.name))
}

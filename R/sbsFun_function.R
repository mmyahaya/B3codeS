#' @title Create species-by-site matrix
#' 
#' @description Create species-by-site (sbs) matrix over the given region and
#' time. The rows and columns of the matrix represents the sites and species
#' respectively. 
#' @param taxa_cube List. A list containing `sim_cubes`of the focal and reference
#' taxa.
#' @param region.sf sf object. The shapefile of the region of study
#' @param res Numeric. The resolution of grid cells to be used. Default is 0.25
#' @param col_temporal Character. The name of the column containing the temporal
#' dimension, e.g., "year". Default is NULL
#'
#' 
#' 
#' sbs<-sbsFun(taxa_cube,KZN,appendix=TRUE)

sbsFun <- function(taxa_cube,
                   region.sf,
                   res=0.25,
                   col_temporal=NULL,
                   appendix=FALSE){
  
  
  if(appendix){
    # extract unique species name from GBIF occurrence data
    species_list<-sort(unique(taxa_cube$taxa$data$scientificName))
    
    sbsM<-taxa_cube$taxa$data %>%
      dplyr::select(scientificName,cellCode,obs) %>%
      dplyr::group_by(scientificName,cellCode) %>%
      dplyr::summarise(across(obs, sum), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
      dplyr::arrange(cellCode) %>% 
      tibble::column_to_rownames(var = "cellCode") %>% 
      mutate_all(~ replace(., is.na(.), 0))
    
    #create species-by-site-by-time (sbsbt) if col_temporal is provided
    sbsbtM<-NULL
    if(!is.null(col_temporal)){
      sbsbtM<-list()
      period<-unlist(unique(taxa_cube$taxa$data[col_temporal]))
      for(t in period){
        sbsbtM[[as.character(t)]]<-taxa_cube$taxa$data %>%
          filter(.data[[col_temporal]]==t) %>% 
          dplyr::select(scientificName,cellCode,obs) %>%
          dplyr::group_by(scientificName,cellCode) %>%
          dplyr::summarise(across(obs, sum), .groups = "drop") %>%
          tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
          dplyr::arrange(cellCode) %>% 
          tibble::column_to_rownames(var = "cellCode") %>% 
          as.matrix()
      }
      
    }
    
    #create grid for region
    grid <- region.sf %>%
      sf::st_make_grid(cellsize = c(res,res),
                       offset = c(sf::st_bbox(region.sf)$xmin,
                                  sf::st_bbox(region.sf)$ymin)) %>%
      sf::st_sf() %>%
      dplyr::mutate(cellid = dplyr::row_number())
    
    # get coordinates of the occurrence sites
    coords <- sf::st_coordinates(sf::st_centroid(grid))
    coords <- coords[as.integer(rownames(sbsM)),]
    colnames(coords)<-c("Longitude","Latitude")
    
  
    # create binary matrix
    sbsM.binary<-sbsM %>% 
      mutate(across(everything(), ~ ifelse(. >= 1, 1, .)))
    
    
    #create site uncertainty
    site_uncertainty <- taxa_cube$taxa$data %>% 
      dplyr::select(all_of(c("cellCode","minCoordinateUncertaintyInMeters"))) %>% 
      dplyr::group_by(cellCode) %>% 
      dplyr::summarise(across(minCoordinateUncertaintyInMeters, sum), .groups = "drop")
    
    
    # compute sbs for ref if it is different from taxa
    if(identical(taxa_cube$taxa,taxa_cube$ref)){
      sbsM.ref<-sbsM
    } else{
      
      sbsM.ref<-taxa_cube$ref$data %>%
        dplyr::select(scientificName,cellCode,obs) %>%
        dplyr::group_by(scientificName,cellCode) %>%
        dplyr::summarise(across(obs, sum), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
        dplyr::arrange(cellCode) %>% 
        tibble::column_to_rownames(var = "cellCode") %>% 
        mutate_all(~ replace(., is.na(.), 0))
      
      colnames(sbsM.ref)<-NULL # remove column names
    }
    
    return(list("sbs"=sbsM,"sbs.ref"=sbsM.ref,"sbs.binary"=sbsM.binary,
                "coords"=coords,"species_list"=species_list,
                "siteID"=rownames(sbsM),"site_unc"=site_uncertainty,
                "sbsbt"=sbsbtM))
  } else{
    
    sbsM<-taxa_cube$taxa$data %>%
      dplyr::select(scientificName,cellCode,obs) %>%
      dplyr::group_by(scientificName,cellCode) %>%
      dplyr::summarise(across(obs, sum), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
      dplyr::arrange(cellCode) %>% 
      tibble::column_to_rownames(var = "cellCode") %>% 
      as.matrix()
    
    #create grid for region
    grid <- region.sf %>%
      sf::st_make_grid(cellsize = c(res,res),
                       offset = c(sf::st_bbox(region.sf)$xmin,
                                  sf::st_bbox(region.sf)$ymin)) %>%
      sf::st_sf() %>%
      dplyr::mutate(cellid = dplyr::row_number())
    
    # get coordinates of the occurrence sites
    coords <- sf::st_coordinates(sf::st_centroid(grid))
    coords <- coords[as.integer(rownames(sbsM)),]
    colnames(coords)<-c("Longitude","Latitude")
    
    #create species-by-site-by-time (sbsbt) if col_temporal is provided
    sbsbtM<-NULL
    if(!is.null(col_temporal)){
      sbsbtM<-list()
      period<-unlist(unique(taxa_cube$taxa$data[col_temporal]))
      for(t in period){
        sbsbtM[[as.character(t)]]<-taxa_cube$taxa$data %>%
          filter(.data[[col_temporal]]==t) %>% 
          dplyr::select(scientificName,cellCode,obs) %>%
          dplyr::group_by(scientificName,cellCode) %>%
          dplyr::summarise(across(obs, sum), .groups = "drop") %>%
          tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
          dplyr::arrange(cellCode) %>% 
          tibble::column_to_rownames(var = "cellCode") %>% 
          as.matrix()
      }
      
    }

    return(list("sbs"=sbsM,
                "coords"=coords,
                "siteID"=rownames(sbsM)))
  }
  
  
  
}


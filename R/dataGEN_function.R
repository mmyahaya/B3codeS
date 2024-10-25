#' @title Prepare data for invasib()
#' 
#' @description
#' Prepare site-by-species (sbs), species-by-trait (sbt) and 
#' site-by-environment (sbe) matrices to calculate invasibility cube. 
#' The function `dataGEN` computes and returns appendix data if specified
#' by the user.
#' 
#' @param taxa Character or dataframe. The character should be the scientific 
#' name of the focal taxa while the dataframe is the GBIF occurrences data which must 
#' contain "decimalLatitude","decimalLongitude","species","speciesKey",
#' "coordinateUncertaintyInMeters","dateIdentified", and "year".
#' @param region.sf sf object. The shapefile of the region of study
#' @param option Numeric. 1 to compute and return sbs, sbt and sbe. 2 to compute
#' and return sbs and sbe.
#' @param limit Integer. Number of records to return from GBIF download.
#' Default is 500
#' @param ref Character or dataframe.The character should be the scientific name 
#' of the reference taxa while the dataframe is the GBIF occurrences data which 
#' must contain "decimalLatitude","decimalLongitude","species","speciesKey",
#' "coordinateUncertaintyInMeters","dateIdentified", and "year".
#' @param country Character. Country for which the GBIP occurrences data should 
#' be downloaded. Country should be 2-letter country code (ISO-3166-1).
#' Default is ZA
#' @param res Numeric. The resolution of grid cells to be used. Default is 0.25
#' @param tryfile Character. File path of the *.txt of the TRY 'released' data or 
#' dataframe imported by rtry::rtry_import()
#'
#' @return A list containing the containing sbs, sbt, sbe matrices and their 
#' appendices if specified.
#'
#'#' @examples
#' \dontrun{
#' KZN_sf<-readRDS(paste0(getwd(),"/KZN_sf.rds"))
#' taxa_Fabaceae_KZN <- readRDS(paste0(getwd(),"/Fabaceae_KZN.rds"))
#' species_trait<-readRDS(paste0(getwd(),"/TRY_traits.rds"))
#' rast_path <- "C:/Users/mukht/Documents"
#' datalist <- dataGEN(taxa=taxa_Fabaceae_KZN,
#'                          region.sf=KZN_sf,
#'                          option=1,
#'                          limit=500,
#'                          ref=NULL,
#'                          res=0.25,
#'                          tryfile=species_trait,
#'                          country='ZA',
#'                          rastfile=rast_path,
#'                          appendix=TRUE)
#' }
#' 

dataGEN = function(taxa,
                   region.sf,
                   option,
                   limit=500,
                   ref=NULL,
                   res=0.25,
                   tryfile,
                   country='ZA',
                   stateProvince=NULL,
                   rastfile,
                   appendix=FALSE){
  

  
  if(option==1){
    
    taxa_cube <- taxaFun(taxa=taxa,
                         region.sf=region.sf,
                         limit=limit, 
                         ref=ref,
                         country=country,
                         res=res,
                         stateProvince=stateProvince)
    
    sbs <- sbsFun(taxa_cube = taxa_cube,
                  region.sf = region.sf, 
                  res = res,
                  appendix = appendix)
    
    sbt <- sbtFun(tryfile = tryfile,
                  taxa_cube = taxa_cube,
                  appendix=appendix)
    
    siteID <- sbs$siteID
    sbe <- sbeFun(rastfile=rastfile,
                  region.sf=region.sf,
                  res=res,siteID=siteID,
                  appendix=appendix)
    
    return(list("sbs"=sbs,"sbt"=sbt,"sbe"=sbe))
  } else if(option==2){
    
    taxa_cube <- taxaFun(taxa=taxa,
                         region.sf=region.sf,
                         limit=limit, 
                         ref=ref,
                         country=country,
                         res=res,
                         stateProvince=stateProvince)
    
    sbs <- sbsFun(taxa_cube = taxa_cube,
                  region.sf = region.sf, 
                  res = res,
                  appendix = appendix)
    
    siteID <- sbs$siteID
    sbe <- sbeFun(rastfile=rastfile,
                  region.sf=region.sf,
                  res=res,
                  siteID=siteID)
    
    return(list("sbs"=sbs,"sbe"=sbe))
  } else{
    message("option has to be 1 (returns sbs,sbt,sbe) or 2 (returns sbs and sbe)")}
  
  
}


datalist <- dataGEN(taxa=taxa_Fabaceae_KZN,
                         region.sf=KZN_sf,
                         option=1,
                         limit=500,
                         ref=NULL,
                         res=0.25,
                         tryfile=species_trait,
                         country='ZA',
                         rastfile=rast_path,
                         appendix=TRUE)

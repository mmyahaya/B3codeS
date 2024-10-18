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
#' @param country.sf sf object. The shapefile of the region of study
#' @param option Numeric. 1 to compute and return sbs, sbt and sbe. 2 to compute
#' and return sbs and sbe.
#' @param limit Integer. Number of records to return from GBIF download.
#' Default is 500
#' @param ref Character or dataframe.The character should be the scientific name 
#' of the reference taxa while the dataframe is the GBIF occurrences data which 
#' must contain "decimalLatitude","decimalLongitude","species","speciesKey",
#' "coordinateUncertaintyInMeters","dateIdentified", and "year".
#' @param res Numeric. The resolution of grid cells to be used. Default is 0.25
#' #' @param tryfile Character. File path of the *.txt of the TRY 'released' data or 
#' dataframe imported by rtry::rtry_import()
#'
#' @return A list containing the containing sbs, sbt, sbe matrices and their 
#' appendices if specified.
#'

dataGEN = function(taxa,
                   country.sf,
                   option,
                   limit=500,
                   ref=NULL,
                   res=0.25,
                   tryfile,
                   country='ZA',
                   rastfile,
                   appendix=FALSE){
  

  
  if(option==1){
    
    taxa_cube <- taxaFun(taxa=taxa,
                         country.sf=country.sf,
                         limit=limit, 
                         ref=ref,
                         country=country,
                         res=res)
    
    sbs <- sbsFun(taxa_cube = taxa_cube,
                  country.sf = country.sf, 
                  res = res,
                  appendix = appendix)
    
    sbt <- sbtFun(tryfile = tryfile,
                  taxa_cube = taxa_cube,
                  appendix=appendix)
    
    siteID <- sbs$siteID
    sbe <- sbeFun(rastfile=rastfile,
                  country.sf=country.sf,
                  res=res,siteID=siteID,
                  appendix=appendix)
    
    return(list("sbs"=sbs,"sbt"=sbt,"sbe"=sbe))
  } else if(option==2){
    
    taxa_cube <- taxaFun(taxa=taxa,
                         country.sf=country.sf,
                         limit=limit, 
                         ref=ref,
                         country=country,
                         res=res)
    
    sbs <- sbsFun(taxa_cube = taxa_cube,
                  country.sf = country.sf, 
                  res = res,
                  appendix = appendix)
    
    siteID <- sbs$siteID
    sbe <- sbeFun(rastfile=rastfile,
                  country.sf=country.sf,
                  res=res,
                  siteID=siteID)
    
    return(list("sbs"=sbs,"sbe"=sbe))
  } else{
    message("option has to be 1 (returns sbs,sbt,sbe) or 2 (returns sbs and sbe)")}
  
  
}

datalist <- dataGEN(taxa=taxa_Fabacae,
                   country.sf=SA.sf,
                   option=1,
                   limit=500,
                   ref=NULL,
                   res=0.25,
                   tryfile=try_path,
                   country='ZA',
                   rastfile=rast_path,
                   appendix=TRUE)

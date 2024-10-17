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
                   appendix=FALSE)

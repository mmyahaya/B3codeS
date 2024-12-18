
#' @taxa (target taxa) can be scientific name(chracter) or GBIF occurrence data(dataframe) 
#' @limit is number of records to return, 500 is the default.
# ref (reference taxa) can be scientific name or GBIF occurrence data. 
#' ref is set to taxa if not provided
# country should be 2-letter country code (ISO-3166-1). Default is ZA
# country.sf should be sf
# res is the resolution of grid cells to be used. Default is 0.25

#' @tryfile should be file path of the *.txt of the TRY 'released' data of 
#' dataframe imported by rtry::rtry_import()

#' @rastfile should be path for storing the worldclim data 
#' (see geodata::worldclim_global) or any other raster data of environmental variable





#precdata <- rast('prec_2021-2040.tif')
#rastpath<-"C:\\Users\\mukht\\Documents"
#taxa_Acacia<-read.csv("taxa_Acacia.csv",sep="")
#taxa_Indigofera<-read.csv("taxa_Indigofera.csv")
#trydata see ?rtry_import for Usage
#rsa_country_sf<- st_read()


rsa_country_sf = sf::st_read("C:/Users/mukht/Documents/boundary_SA/boundary_south_africa_land_geo.shp")


datalist<-dataGEN(taxa=taxa_Fabacae,country.sf=SA.sf,
                  tryfile=try_path,rastfile=rastpath,limit = 1000)

datalist2<-dataGEN(taxa="Vachellia",country.sf=rsa_country_sf,
                  tryfile="TRY_All.txt",rastfile=rastpath,limit = 2000)


# datalist is a nested list.
# datalist[["sbs"]]$sbs  : site by species
# datalist[["sbt"]]$traitname ; Name of trait

taxa.sf<-taxaFun(taxa=taxa_Fabacae,country.sf = rsa_country_sf, ref = taxa_Acacia)
sbs <- sbsFun(taxa_cube = A)
sbt<-sbtFun(tryfile =  "TRY_Acacia.txt",taxa_cube = taxa_cube)
sbe<-sbeFun(rastfile =  precdata,country.sf =  rsa_country_sf,siteID = sbs$siteID)


# taxa can be specie name or GBIF occurrence data 
# Limit is number of records to return, 500 is the default.
# ref can be scientific name or GBIF occurrence data. ref is set to taxa if not provided
# Country should be 2-letter country code (ISO-3166-1).


taxa.sf <- taxaFun(taxa = "Acacia")



sbs<-sbsFun(taxa.sf = taxa.sf, country.shp = rsa_country_sf)

path = "C:/Users/mukht/Documents" #path for worldclim

precdata <- rast('prec_2021-2040.tif')


sbe<-sbeFun(rastfile= chelsaA18, country.shp = rsa_country_sf )

try_path<-"33852.txt" # path for trydata

sbt<-sbtFun(tryfile = try33576,taxa.sf = taxa.sf$taxa)



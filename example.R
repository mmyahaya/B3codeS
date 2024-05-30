
# taxa (target taxa) can be scientific name(chracter) or GBIF occurrence data(dataframe) 
# limit is number of records to return, 500 is the default.
# ref (reference taxa) can be scientific name or GBIF occurrence data. ref is set to taxa if not provided
# country should be 2-letter country code (ISO-3166-1). Default is ZA
# country.sf should be sf
# res is the resolution of grid cells to be used. Default is 0.25
# tryfile should be file path of the *.txt of the TRY 'released' data of dataframe imported by rtry::rtry_import()
# rastfile should be path for storing the worldclim data (see geodata::worldclim_global) or any other raster data of environmental variable


#

#precdata <- rast('prec_2021-2040.tif')
#taxa_indi<-readr::read_csv("taxa_Indigofera.csv")
#trydata see ?rtry_import for Usage
#rsa_country_sf<- st_read()
datalist<-dataGEN(taxa=taxa_Acaia,country.sf=rsa_country_sf,ref=taxa_Indigofera,
                  tryfile=trydata,rastfile=precdata)








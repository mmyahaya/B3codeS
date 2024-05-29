taxa.sf <- taxaFun(taxa = taxa_df, ref = ref.data)



sbs<-sbsFun(taxa.sf = taxa.sf, country.shp = rsa_country_sf)

path = "C:/Users/mukht/Documents" #path for worldclim

precdata <- rast('prec_2021-2040.tif')


sbe<-sbeFun(rastfile= chelsaA18, country.shp = rsa_country_sf )

try_path<-"33852.txt" # path for trydata

sbt<-sbtFun(tryfile = try33576,taxa.sf = taxa.sf$taxa)



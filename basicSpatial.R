
# Load libraries
library(terra)
library(sf)
library(sp)

library(ambient)

library(rasterVis)
library(colorspace)
library(RColorBrewer)

# Read RSA land area shapefile
(rsa_country_sf = st_read("C:/Users/mukht/OneDrive/Documents/boundary_SA/boundary_south_africa_land_geo.shp"))
plot(rsa_country_sf['Land'])

# Create a QDS raster
(gridQDS = rast(rsa_country_sf,res=c(0.25,0.25), crs='+init=epsg:4326'))
gridQDS$qdsID = 1:ncell(gridQDS)
plot(gridQDS)

# Make QDS Mask. Remember to mask the background to NA
rsa_mask = rasterize(rsa_country_sf, gridQDS, background=NA)
plot(rsa_mask)

# Mask the QDS raster to check it's working
gridQDS_mask = mask(gridQDS, rsa_mask)
plot(gridQDS_mask[[uN[1:4]]])
plot(dplyr::filter(taxa.sf, species==uN[4]),add=TRUE)
# Create 100 random points across South Africa
random_pts = st_sample(rsa_country_sf, size=1000, type="random")    
random_pts_sf = vect(random_pts)
points(random_pts_sf, pch=20)
lines(rsa_country_sf['Land'])

# Extract data to points
qds_values = extract(gridQDS_mask, random_pts_sf, xy=TRUE, bind=TRUE)
qds_values

# Create an empty raster with finer scale than QDS
# egEnv1 = rast(rsa_country_sf['Land'],res=c(0.0625,0.0625), crs='+init=epsg:4326')
egEnv1 = disagg(gridQDS,fact=4)

# Fill the raster with random values
set.seed(123)  # for reproducibility
egEnv1$egBio = runif(ncell(egEnv1), min=0, max=10)  # generate random values and assign the values to the raster cells

# Check the raster
plot(egEnv1$egBio)
lines(rsa_country_sf['Land'])

# INSTEAD CREATE A RASTER WITH PATTERN/PATCHES
# Generate a Perlin noise based matrix
perlin_matrix = noise_perlin(dim = dim(egEnv1$egBio), frequency = 0.1)

# Normalize the matrix to range [0, 1]
habitat = (perlin_matrix - min(perlin_matrix)) / (max(perlin_matrix) - min(perlin_matrix))

# Convert the matrix to a new layer in your SpatRaster
egEnv1$egBioP = rast(egEnv1$egBio, vals=habitat)

# Plot the resulting raster with patches
plot(egEnv1$egBioP, main="Raster with Patches")
lines(rsa_country_sf['Land'],lwd=3, col="blue")

# Extract data to points (SpatVector)
env1_values = extract(egEnv1, qds_values, bind=TRUE)
env1_values
points(env1_values, pch=20)

# Check frequency of records in QDS
table(env1_values$qdsID)

# Better map using rasterVis
# plot(egEnv1$egBioP)
# lines(rsa_country_sf['Land'])
# points(env1_values, pch=20)

(brks = seq(0, 1, by=0.1))
(nb = length(brks)-1)
colsN = brewer.pal(n = nb, name = "Spectral")
# cols6 = brewer.pal(n = 6, name = "Spectral")
# cols10 = rev(brewer.pal(n = 10, name = "Spectral"))
# myPal1 = terrain.colors(nb)
# myPal2 = diverge_hcl(10,c=100,l=c(50,90),power=1)

levelplot(egEnv1$egBioP, at = brks, col.regions = colsN,
          main='A better raster map v1') +
  latticeExtra::layer(sp.lines(as(rsa_country_sf['Land'], 'Spatial')))



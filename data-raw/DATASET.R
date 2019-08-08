## code to prepare raster data frames in mase/data directory
library(raster)
#shape files found at https://data.fs.usda.gov/geodata/edw/datasets.php?xmlKeyword=Ecomap
province_shape <- sf::read_sf("~/summer 19/forestry-ds/data/bigDATA/S_USA.EcoMapProvinces.shp")
#image files found at https://data.fs.usda.gov/geodata/rastergateway/biomass/index.php
black_hills_pop_bio <- raster::raster("conus_forest_biomass_mg_per_ha.img")
black_hills_pop_for <- raster::raster("conus_forest_nonforest_probability.img")
black_hills_pop_fg <- raster::raster("conus_forestgroup.img")

#get Black Hills province (M334)
eco_region <- province_shape %>% dplyr::filter(MAP_UNIT_S == "M334") %>%
  sf::st_transform(crs(black_hills_pop_for))

#crop raster files and reduce resolution by 3
black_hills_pop_bio <- crop(black_hills_pop_bio, eco_region) %>%
  mask(eco_region) %>%
  aggregate(fact = 3)
black_hills_pop_for <- crop(black_hills_pop_for, eco_region) %>%
  mask(eco_region) %>%
  aggregate(fact = 3)
black_hills_pop_fg <- crop(black_hills_pop_fg, eco_region) %>%
  mask(eco_region) %>%
  aggregate(fact = 3)

#as data frame
bh_bio <- black_hills_pop_bio %>% 
  projectRaster(crs="+proj=utm +zone=13  +ellps=GRS80 +units=m +no_defs") %>%
  as("SpatialGridDataFrame") %>%
  as.data.frame() %>% rename(biomass = conus_forest_biomass_mg_per_ha)
bh_for <- black_hills_pop_for %>% 
  projectRaster(crs="+proj=utm +zone=13  +ellps=GRS80 +units=m +no_defs") %>%
  as("SpatialGridDataFrame") %>%
  as.data.frame() %>% rename(forest_nonforest_probability = conus_forest_nonforest_probability)
bh_fg <- black_hills_pop_fg %>% 
  projectRaster(crs="+proj=utm +zone=13  +ellps=GRS80 +units=m +no_defs") %>%
  as("SpatialGridDataFrame") %>%
  as.data.frame() %>% rename(forest_group = conus_forestgroup) %>%
  head(nrow(bh_for)) #data doesn't always line up, remove the 3 straggler pixels
#tidy up
black_hills_pop <- cbind(biomass = bh_bio %>% pull(biomass),
                         forest_group = bh_fg %>% pull(forest_group),
                         bh_for)
#remove(list = setdiff(ls(), c("black_hills_pop")))

#create the sample
set.seed(128)
sample <- sample(1:nrow(black_hills_pop), 232)
black_hills_samp <- black_hills_pop[sample,]

#remove response variable (biomass) from population
black_hills_pop <- black_hills_pop[,-1]

#export to /data folder
usethis::use_data(black_hills_pop, black_hills_samp, overwrite = TRUE)

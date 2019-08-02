#' Spatially indexed forest nonforest probabilities for the Black Hills in South Dakota and Wyoming.
#'
#' A dataset containing forest nonforest probabilities across the Black Hills from geospatial satellite raster data.
#'
#' @format A data frame with 24240 rows and 3 variables:
#' \describe{
#'   \item{forest_nonforest_prob}{probability of forest}
#'   \item{s1}{spatial grid index 1}
#'   \item{s2}{spatial grid index 2}
#' }
#' @source \url{https://data.fs.usda.gov/geodata/rastergateway/}
"black_hills_pop"

#' Spatially indexed biomass and forest nonforest probability samples from the Black Hills in South Dakota and Wyoming.
#'
#' A sample of biomass measures and forest nonforest probabilities from geospatial satellite raster data.
#'
#' @format A data frame with 232 rows and 4 variables:
#' \describe{
#'   \item{forest_nonforest_prob}{probability of forest}
#'   \item{biomass}{forest biomass in Megagrams per Hectare}
#'   \item{s1}{spatial grid index 1}
#'   \item{s2}{spatial grid index 2}
#' }
#' @source \url{https://data.fs.usda.gov/geodata/rastergateway/}
"black_hills_samp"
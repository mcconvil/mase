library(dplyr)
raw <- readRDS(url("https://mcconvil.github.io/fia_workshop_2021/data/IDdata.rds","rb"))

# in "means" form currently
IdahoPop <- raw$dunitzonal[ ,c("COUNTYFIPS", "tcc", "elev", "tnt.1", "tnt.2", "npixels")] |>
  filter(!(COUNTYFIPS %in% c("16027", "16047", "16053", "16063", "16067", "16075")))

usethis::use_data(IdahoPop, overwrite = TRUE)

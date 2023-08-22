raw <- readRDS(url("https://mcconvil.github.io/fia_workshop_2021/data/IDdata.rds","rb"))

# in "means" form currently
IdahoPop <- raw$dunitzonal[ ,c("COUNTYFIPS", "tcc", "elev", "tnt.1", "tnt.2", "npixels")]

usethis::use_data(IdahoPop, overwrite = TRUE)

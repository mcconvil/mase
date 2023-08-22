raw <- readRDS(url("https://mcconvil.github.io/fia_workshop_2021/data/IDdata.rds","rb"))

IdahoSamp <- raw$pltassgn[ ,c("COUNTYFIPS", "tcc", "elev", "ppt", "tmean", "tnt", "BA_TPA_ADJ")]

usethis::use_data(IdahoSamp, overwrite = TRUE)

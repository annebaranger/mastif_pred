library(targets)

Sys.setenv("TAR_PROJECT"="data")
tar_visnetwork()
# tar_make()


Sys.setenv("TAR_PROJECT"="data_gbif")
tar_visnetwork()
tar_make()

Sys.setenv("TAR_PROJECT"="data_nfi")
tar_visnetwork()
tar_make()

Sys.setenv("TAR_PROJECT"="fit")
tar_visnetwork()
# tar_make()

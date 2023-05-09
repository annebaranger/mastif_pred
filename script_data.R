# _targets.R file

#library
library(targets)
# lapply(c("stringr","ggplot2","tidyr","dplyr"),require,character.only=TRUE)


#Options
source("R/functions_data.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","tidyr","dplyr","terra"),
               error = "continue") 

#Targets
list(
  tar_target(
    mastif.data,
    get_mastif_europe(dir.data="data/fecundityMastif.rdata")
  ),
  NULL
)
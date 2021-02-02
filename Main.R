
#Install packages

#Packages to be used
packages<-c("readxl","here","tidyverse","ggplot2","gridExtra","knitr","BRugs","coda","rjags")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#Create Dir.
dir.create(here("Figure"))
dir.create(here("Results"))

#Run scripts
source(here("Data","data.R"))

source(here("Script","bayes1.R"))

source(here("Script","bayes2.R"))

source(here("Script","model.R"))

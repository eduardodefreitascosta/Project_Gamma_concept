
if(!require(here)){
  install.packages("here")
}
library(here)

dir.create(here("Figura"))
dir.create(here("Resultado"))


source(here("Dados","data.R"))

source(here("Script","gabi_bayes1.R"))

source(here("Script","gabi_bayes2.R"))

source(here("Script","modelo_gabi.R"))


if(!require(here)){
  install.packages("here")
}
library(here)


if(!require(ggplot2)){
  install.packages("ggplot2")
}
library(ggplot2)

if(!require(gridExtra)){
  install.packages("gridextra")
}
library(gridExtra)

if(!require(knitr)){
  install.packages("knitr")
}
library(knitr)


dir.create(here("Figura"))
dir.create(here("Resultado"))


source(here("Dados","data.R"))

source(here("Script","gabi_bayes1.R"))

source(here("Script","gabi_bayes2.R"))

source(here("Script","modelo_gabi.R"))

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(GSA.UN)

## ----echo=TRUE, results='hide'------------------------------------------------
library(GSA.UN)
library(stats)
library(e1071)

data("out_set")

data_Bstat <- Bstat(out_set)


## ----echo=TRUE, results='hide'------------------------------------------------
library(GSA.UN)
library(stats)
library(e1071)

data("parameters_set", "out_set", "pp_names")

CM <- Cond_Moments(parameters_set, out_set, pp_names, steps=15)

## ----echo=TRUE, results='hide'------------------------------------------------
library(GSA.UN)
library(stats)

data("data_Bstat", "CM", "pp_names")

SOBOL_indices <- SOBOL(data_Bstat[,3], CM$CM_mean, CM$CM_var , pp_names)

## ----echo=TRUE, results='hide'------------------------------------------------
library(GSA.UN)

data("data_Bstat", "CM", "pp_names")

AMA_indices <- AMA(data_Bstat, CM, pp_names, steps= 15)

## ----echo=TRUE, results='hide'------------------------------------------------
library(GSA.UN)
library(stats)
library(e1071)

data("parameters_set", "out_set", "pp_names")

GSA_results <- GSAtool(parameters_set, out_set, pp_names, steps = 15, save=FALSE)


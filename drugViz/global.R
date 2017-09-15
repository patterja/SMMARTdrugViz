library(data.table)
library(ggplot2)
library(PharmacoGx)
library(gridExtra)


load("./data/auc_erlo.Rda")
load("./data/exp89.Rda")
load("./data/expEgfr.Rda")

medianWithoutNA <- function(x) {median(x[which(!is.na(x))])}


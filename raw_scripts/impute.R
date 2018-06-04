## Script to run imputation in parallel
## Each one produces 5 datasets

library(dplyr)
library(mice)
## Define command line arguments
args = commandArgs(trailingOnly = TRUE)

## Parameters
ITERATION = args[1]

load("/home/mourikisa/data/geneProperties_final.Rdata")

geneProperties$duplicated = factor(geneProperties$duplicated)
geneProperties$WGD = factor(geneProperties$WGD)
geneProperties$hub = factor(geneProperties$hub)
geneProperties$central = factor(geneProperties$central)
geneProperties$age = factor(geneProperties$age)
geneProperties$origin = factor(geneProperties$origin)
geneProperties$exp.breadth = factor(geneProperties$exp.breadth)

geneProperties_miceImputed <- mice(geneProperties[2:15], m=10, maxit = 20)

save(geneProperties_miceImputed, file = paste0("/home/mourikisa/data/geneProperties_final_miceImputed_object_", ITERATION, ".Rdata"))





###############################################################################
##  Predicting cancer genes on the validation cohort using sysSVM best models
###############################################################################

detach(name = "package:dplyr", unload = T)
library(plyr)
library(dplyr)
source("~/rosalind_lustre/mourikisa/data/OAC/functions.R")

ns = c("nonsynonymous","stopgain","frameshift deletion","splicing","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","stoploss","frameshift substitution")
dam = c("nonsynonymous","frameshift deletion","frameshift insertion","frameshift substitution","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","splicing","stopgain","stoploss")
trunc = c("frameshift deletion","frameshift insertion","frameshift substitution","stopgain","stoploss") ## Always damaging==TRUE
non_trunc = c("nonsynonymous","splicing")

## First load the data
## muts
load("~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/Rdata/muts_21_literature_OACs_annovar_dbnsfp_oncodriveClust.Rdata")

## CNVs
load("~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/Rdata/cnvs_21_literature_OACs.Rdata")

## Get total table
total_table = createTotalTable(muts=muts, cnvs = cnvs)

## Get ML input
oac_ML_input = getMLinput(df=total_table)
write.table(oac_ML_input, file="~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/sysSVM/oac_21_literature_ML_input.tsv", row.names = F, sep = "\t", quote=F)

load("~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/sysSVM/OAC/training_set.Rdata")
load("~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/sysSVM/OAC/validation_set.Rdata")
prediction_set = rbind(training[,!names(training)%in%c("type")], validation)
save(prediction_set, file="~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/sysSVM/prediction_set.Rdata")

load("~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/sysSVM/OAC/training_set_noScale.Rdata")
load("~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/sysSVM/OAC/validation_set_noScale.Rdata")
prediction_set_ns = rbind(training_ns[,!names(training_ns)%in%c("type")], validation_ns)
save(prediction_set_ns, file="~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/sysSVM/prediction_set_noScale.Rdata")



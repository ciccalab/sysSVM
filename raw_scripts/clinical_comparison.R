rm(list=ls())

library(survival)
library(survminer)
library(dplyr)
library(OIsurv) # Aumatically loads KMsurv
library(ranger)
library(ggplot2)
library(KMsurv)
library(igraph)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

## go to the working directory
setwd("athena/data/OAC/")
source("functions.R")

######################################################
##                      Groups
######################################################

load("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/patient_stratification_onlyTop10.Rdata")
d = d %>% select(sample, group_pathway_numbers) %>% unique

######################################################
##              Get the clinical data
######################################################

clinic = read.table("~/rosalind_lustre/mourikisa/data/OAC/Clinical_data/261_OAC_clinical.txt", header = T, sep="\t")

## Unmapping the values (Lawrence gave me the format thathe uploads to the ICGC web site)
clinic = clinic %>% mutate(specimen_donor_treatment_type=ifelse(specimen_donor_treatment_type==1, "no_treatment",
                                                                ifelse(specimen_donor_treatment_type==2, "chemotherapy",
                                                                       ifelse(specimen_donor_treatment_type==3, "radiation_therapy", 
                                                                              ifelse(specimen_donor_treatment_type==4, "chemo plus radiation therapy", 
                                                                                     ifelse(specimen_donor_treatment_type==5, "immunotherapy", 
                                                                                            ifelse(specimen_donor_treatment_type==6, "chemo plus immunotherapy", 
                                                                                                   ifelse(specimen_donor_treatment_type==7, "surgery", 
                                                                                                          ifelse(specimen_donor_treatment_type==8, "other_therapy", 
                                                                                                                 ifelse(specimen_donor_treatment_type==9, "bone marrow transplant", 
                                                                                                                        ifelse(specimen_donor_treatment_type==10, "stem cell transplant", 
                                                                                                                               ifelse(specimen_donor_treatment_type==11, "monoclonal antibodies", NA))))))))))))
clinic = clinic %>% mutate(donor_sex=ifelse(donor_sex==1, "male", 
                                            ifelse(donor_sex==2, "female", NA)))
clinic = clinic %>% mutate(tobacco_smoking_history_indicator=ifelse(tobacco_smoking_history_indicator==1, "never-smoker", 
                                                                    ifelse(tobacco_smoking_history_indicator==2, "smoker", 
                                                                           ifelse(tobacco_smoking_history_indicator%in%c(3:5), "non-smoker", NA))))

clinic = clinic %>% mutate(donor_vital_status=ifelse(donor_vital_status==1, "alive", 
                                                     ifelse(donor_vital_status==2, "deceased", NA)))


clinic = clinic %>% select(sample, specimen_donor_treatment_type, tumour_grade, 
                           tumour_stage_revised, donor_sex, donor_age_at_diagnosis, donor_survival_time,
                           donor_vital_status,
                           tobacco_smoking_history_indicator, siewert_type)
clinic[clinic== -888] = NA

######################################################
##                      Comparisons
######################################################

## get group number in the table
s = clinic %>% left_join(d) %>% select(donor_survival_time, donor_vital_status, group_pathway_numbers) %>% 
    mutate(status=ifelse(donor_vital_status=="alive", 0, 1)) %>% rename(group=group_pathway_numbers)

## Survival analysis
s = s %>% subset(donor_survival_time!="NA") %>% mutate(donor_survival_time=as.numeric(donor_survival_time))

ss = s %>% mutate(group=ifelse(group=="6", "6", "rest"))

## One for all groups together
my.surv <- Surv(time=as.numeric(ss$donor_survival_time), event = ss$status, type="right")
my.fit = survfit(my.surv~group, data = ss)
pdf("~/Desktop/survival_6_vs_rest.pdf")
ggsurvplot(my.fit,
           pval = TRUE, conf.int = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("green", "blue", "red", "black", "cyan", "orange"))

#plot(my.fit)
dev.off()




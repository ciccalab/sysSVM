#!/usr/bin/Rscript

## Script to run oncodriverClust
## Separate because I need to run all the samples together

setwd("/mnt/lustre/users/k1469280/mourikisa/data/OAC")
cat("Sourcing relevant functions...")
source("functions.R")


OUT = commandArgs(trailingOnly = TRUE)[1]


## Mutations
# mainDirs = c("~/data/OAC/71_OAC/strelka/",
#              "~/data/OAC/87_OAC/66_ICGC/strelka/")
mainDirs = c("/mnt/lustre/users/k1469280/mourikisa/data/OAC/87_OAC/21_literature/strelka/")

message("Getting mutations...")
all_muts = list()
count = 0
for(b in mainDirs){
    samples = list.dirs(b, recursive = F)
    for(s in samples){
        cat(s, "\n")
        muts_fn = paste0(s, "/parsing_and_annotation/annovar/muts_annovar_dbnsfp_19014.Rdata")
        load(muts_fn)
        sname = unlist(strsplit(s, "/"))
        sname = sname[length(sname)]
        all_muts[[sname]] = muts %>% mutate(sample=sname)
        count = count +1
    }
}
cat(paste0("Samples: ", count), "\n")

muts = do.call(rbind.fill, all_muts)
save(muts, file=paste0(OUT, "/Rdata/muts_129_66_71_OACs_annovar_dbnsfp_19014.Rdata"))

## Add also the 129 from before (uncomment the following two lines if you run the main cohort)
#message("Getting 129 previous samples...")
#load("~/data/OAC/129_OAC/Rdata/mutations_annotated_19014.Rdata")

## If you want one table (uncomment the following two lines if you run the main cohort)
#muts2 = do.call(rbind.fill, muts)
#muts2 = muts2 %>% select(-oncodriveClust)

## Fix the names of the samples (uncomment the following three lines if you run the main cohort)
#load("~/data/OAC/129_OAC/Rdata/mainCohort.Rdata")
#mainCohort = mainCohort %>% select(directory, Primary_tumour_solid_tissue) %>% rename(sample=Primary_tumour_solid_tissue)
#muts2 = muts2 %>% left_join(mainCohort) %>% select(-sample) %>% rename(sample=directory)

#m = rbind.fill(muts1, muts2)
m = muts
#rm(muts1)
#rm(muts2)
rm(muts)
rm(all_muts)

muts = runOncodriveClust(muts=m, save_dir=paste0(OUT, "/oncodriveClust"))
save(muts, file=paste0(OUT, "/Rdata/muts_129_66_71_OACs_annovar_dbnsfp_19014_oncodriveClust.Rdata"))










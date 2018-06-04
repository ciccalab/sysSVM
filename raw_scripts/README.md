## This directory was created to automate the training and prediction steps of sysSVM

## ***************************************************************************
##                  Parsing of data / Mutation annotation
## ***************************************************************************

## ------------------------------
## Sample parsing and annotation
## ------------------------------

Five scripts involved in this step
functions.R ## Functions for mutation annotation, gene annotation etc
data_parsing.R ## R script that does the data manipulation and mutation annotation
data_parsing.sh ## Submission script - calls the R script
parse_dir_str.R ## To make the commands - see below
runOncodriveClust.R ## To run GOF prediction (separate because I need to aggregate mutations form all the samples together) - I usually run this via Rscript in rosalind (activating oncodroveClust pyenv first)
summary_plotting.R ## Use this to get the overall data (aggregation of all patients) and save mut, cnv, sv data for all samples


## Example command (I usually make all these commands with R and write them in a file)
qsub -N annot -pe threaded 3 -q FCgroup.q -l h_vmem=12G -o /home/mourikisa/data/OAC/71_OAC/strelka/LP6005334-DNA_B01_vs_SM-4JPUU/test.out -e /home/mourikisa/data/OAC/71_OAC/strelka/LP6005334-DNA_B01_vs_SM-4JPUU/test.err ./mutation_annotation.sh -s /home/mourikisa/data/OAC/71_OAC/strelka/LP6005334-DNA_B01_vs_SM-4JPUU -m /home/mourikisa/data/OAC/71_OAC/strelka/LP6005334-DNA_B01_vs_SM-4JPUU/pv.1.3__rg.grch37_e71__al.bwa/strelka/filtered_results/LP6005334-DNA_B01_vs_SM-4JPUU.snp.pass.vcf -i /home/mourikisa/data/OAC/71_OAC/strelka/LP6005334-DNA_B01_vs_SM-4JPUU/pv.1.3__rg.grch37_e71__al.bwa/strelka/results/passed.somatic.indels.vcf -g hg19


As soon as you have the overall tables use createTotalTable function in the functions.R to create the total table
and then the getMLinput function to bring in systems-level properties


## ***************************************************************************
##                          Main pipeline (sysSVM)
## ***************************************************************************

## ------------------------------------------------------------------------------------------------
## STEP 01: Create the training and prediction set and the command files for training (prepare.sh)
## ------------------------------------------------------------------------------------------------

## The first part of the pipeline is to prepare the training and prediction set of the 4 kernels for all the cancer types
## To do that run the following command for all the cancer types (change the name of the cancer type and the parameters accordingly)
## Commands

## More command-line arguments were added to address both input from file and from a mysql database

qsub -N OAC -pe threaded 3 -q FCgroup.q -l h_vmem=12G ./prepare.sh OAC F /home/mourikisa/data/OAC/129_OAC/ML/oac_ML_input.tsv young,reptime T /home/mourikisa/data/OAC /home/mourikisa/data/OAC/129_OAC/ML/OAC_false_positives_included /home/mourikisa/data/OAC/config.R F

qsub -N OAC -l h_vmem=12G ./prepare.sh OAC F /mnt/lustre/users/k1469280/mourikisa/data/OAC/validation_cohort_OAC_TCGA/oac_TCGA_ML_input.tsv young,reptime T /mnt/lustre/users/k1469280/mourikisa/data/OAC /mnt/lustre/users/k1469280/mourikisa/data/OAC/validation_cohort_OAC_TCGA/sysSVM/ /mnt/lustre/users/k1469280/mourikisa/data/OAC/config.R /mnt/lustre/users/k1469280/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/mean_sd_of_training_validation_features.Rdata

qsub -N OAC -l h_vmem=12G ./prepare.sh OAC F /mnt/lustre/users/k1469280/mourikisa/data/OAC/18_OAC/sysSVM/oac_18_ML_input.tsv young,reptime T /mnt/lustre/users/k1469280/mourikisa/data/OAC /mnt/lustre/users/k1469280/mourikisa/data/OAC/18_OAC/sysSVM/ /mnt/lustre/users/k1469280/mourikisa/data/OAC/config.R F

## ------------------------------
## STEP 02: Run the command files
## ------------------------------

Change to the directory where the cancer type is and submit the files to the cluster (i.e jobs_linear.txt etc...)


## ---------------------------------------------------------------------------------------
## STEP 03: Pool the result, select the best model and calculate the scores for each gene
## ---------------------------------------------------------------------------------------

## Before this stewp remember to check GTEx tissue name for the corresponding cancer type

qsub -N OAC_score -pe threaded 3 -q FCgroup.q -l h_vmem=30G ./score_genes.sh OAC esophagus Esophagus T /home/mourikisa/data/OAC/129_OAC/ML /home/mourikisa/data/OAC/config.R F

## Use this command to predict in a new cohort
qsub -N OAC_score -l h_vmem=30G ./predict_new_samples.sh OAC esophagus Esophagus T /mnt/lustre/users/k1469280/mourikisa/data/OAC/validation_cohort_OAC_TCGA /mnt/lustre/users/k1469280/mourikisa/data/OAC/config.R /mnt/lustre/users/k1469280/mourikisa/data/OAC/validation_cohort_OAC_TCGA/prediction_set.Rdata /mnt/lustre/users/k1469280/mourikisa/data/OAC/validation_cohort_OAC_TCGA/prediction_set_noScale.Rdata /mnt/lustre/users/k1469280/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/best_models.tsv /mnt/lustre/users/k1469280/mourikisa/data/OAC/validation_cohort_OAC_TCGA/sysSVM/OAC/config.txt

## ---------------------------------------
##      Get the ranks in each sample
## ---------------------------------------

## Add refined CGCs, exclude "Not expressed" genes and define patient ranks
## Then I use the getSysCans (source it from config.R) to get a data frame with the ranks. Command as follows:
syscan = getSysCans(path = "~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/", cancer_type = "OAC", tissue = "Esophagus")
## And then save it
save(syscan, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/syscan.Rdata")


## ---------------------------------------
##      gene set enrichment analysis
## ---------------------------------------

gsea.noAmp.top10.plusCGC = pathwayEnrichmentTopGenes() ## With default parameters
save(gsea.noAmp.top10.plusCGC, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.noAmp.top10.plusCGC.Rdata")
gsea.withAmp.top10.plusCGC = pathwayEnrichmentTopGenes(remove.amplifications=F) 
save(gsea.withAmp.top10.plusCGC, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.top10.plusCGC.Rdata")


## -------------------------------------------------------------
##      Sys-candidates + Refined CGCs to enriched pathways
## -------------------------------------------------------------
## This function returns the gsea list above with an extra element which is the genes mapped to pathways (pathways also annotated as Level 1/2 and also with the number of samples/genes they cover)
gsea.noAmp.top10.plusCGC.pathways = mapPathways2Samples(gsea=gsea.noAmp.top10.plusCGC)
save(gsea.noAmp.top10.plusCGC.pathways, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.noAmp.top10.plusCGC.pathways.Rdata")
gsea.withAmp.top10.plusCGC.pathways = mapPathways2Samples(gsea=gsea.withAmp.top10.plusCGC)
save(gsea.withAmp.top10.plusCGC.pathways, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.top10.plusCGC.pathways.Rdata")

## ---------------------------------------
##      Gene sets to containers
## ---------------------------------------
gsea.noAmp.top10.plusCGC.containers = getContainers(gsea=gsea.noAmp.top10.plusCGC.pathways[["genes2paths"]])
save(gsea.noAmp.top10.plusCGC.containers, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.noAmp.top10.plusCGC.containers.Rdata")
gsea.withAmp.top10.plusCGC.containers = getContainers(gsea=gsea.withAmp.top10.plusCGC.pathways[["genes2paths"]])
save(gsea.withAmp.top10.plusCGC.containers, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.top10.plusCGC.containers.Rdata")


## --------------------------------------------------------
##      Gene patient groups using shared containers
## --------------------------------------------------------
common_processes_top5_plusCGC = getPatientCommonProcesses(genes2paths=gsea.withAmp.onlyTop5.containers[["60"]][["genes2paths"]])
sample_combns = common_processes_top5_plusCGC
save(sample_combns, file="~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/common_processes_top5_plusCGC.Rdata")

## -------------------------------------------------------------------
##      GSEA/processes/patient groups on additional datasets
## -------------------------------------------------------------------
## Only CGCs
gsea.withAmp.onlyCGC = pathwayEnrichmentTopGenes(remove.amplifications = F, genes.per.patient = NULL)
save(gsea.withAmp.onlyCGC, file="~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.onlyCGC.Rdata")
gsea.withAmp.onlyCGC.pathways = mapPathways2Samples(gsea=gsea.withAmp.onlyCGC)
save(gsea.withAmp.onlyCGC.pathways, file="~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.onlyCGC.pathways.Rdata")
gsea.withAmp.onlyCGC.containers = getContainers(gsea=gsea.withAmp.onlyCGC.pathways[["genes2paths"]])
save(gsea.withAmp.onlyCGC.containers, file = "~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.onlyCGC.containers.Rdata")


## Only Top 10 syscans
gsea.withAmp.onlyTop10 = pathwayEnrichmentTopGenes(remove.amplifications = F, genes.per.patient = 10, CGC.plus = F, remove.cgc = T)
save(gsea.withAmp.onlyTop10, file="~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.onlyTop10.Rdata")
gsea.withAmp.onlyTop10.pathways = mapPathways2Samples(gsea=gsea.withAmp.onlyTop10)
save(gsea.withAmp.onlyTop10.pathways, file="~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.onlyTop10.pathways.Rdata")
gsea.withAmp.onlyTop10.containers = getContainers(gsea=gsea.withAmp.onlyTop10.pathways[["genes2paths"]])
save(gsea.withAmp.onlyTop10.containers, file = "~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.onlyTop10.containers.Rdata")
common_processes_onlyTop10 = getPatientCommonProcesses(genes2paths=gsea.withAmp.onlyTop10.containers[["60"]][["genes2paths"]])
save(common_processes_onlyTop10, file = "~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/common_processes_onlyTop10.Rdata")

## All predictions
gsea.withAmp.allPredictionSet = pathwayEnrichmentTopGenes(remove.amplifications = F, genes.per.patient = NULL, CGC.plus = NULL, given_tb = validation_ns)
save(gsea.withAmp.allPredictionSet, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.allPredictionSet.Rdata")
gsea.withAmp.allPredictionSet.pathways = mapPathways2Samples(gsea=gsea.withAmp.allPredictionSet)
save(gsea.withAmp.allPredictionSet.pathways, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.allPredictionSet.pathways.Rdata")
gsea.withAmp.allPredictionSet.containers = getContainers(gsea=gsea.withAmp.allPredictionSet.pathways[["genes2paths"]])
save(gsea.withAmp.allPredictionSet.containers, file = "~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.allPredictionSet.containers.Rdata")


## All predictions + CGCs (validation+training sets essentially)
cohort = rbind(training_ns%>%select(-type), validation_ns)
gsea.withAmp.allPredictionSet.plusCGC = pathwayEnrichmentTopGenes(remove.amplifications = F, genes.per.patient = NULL, CGC.plus = NULL, given_tb = cohort)
save(gsea.withAmp.allPredictionSet.plusCGC, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.allPredictionSet.plusCGC.Rdata")
gsea.withAmp.allPredictionSet.plusCGC.pathways = mapPathways2Samples(gsea=gsea.withAmp.allPredictionSet.plusCGC)
save(gsea.withAmp.allPredictionSet.plusCGC.pathways, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.allPredictionSet.plusCGC.pathways.Rdata")
gsea.withAmp.allPredictionSet.plusCGC.containers = getContainers(gsea=gsea.withAmp.allPredictionSet.plusCGC.pathways[["genes2paths"]])
save(gsea.withAmp.allPredictionSet.plusCGC.containers, file = "~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.allPredictionSet.plusCGC.containers.Rdata")




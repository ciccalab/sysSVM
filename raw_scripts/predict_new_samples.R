rm(list=ls())


##########################################################
##                      Score
##########################################################

## Parameters
CANCER_TYPE = "OAC"
TISSUE_NCG = "esophagus"
TISSUE_GTEX = "Esophagus"
GOF=T
GAINS=T
RES_DIR = "~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/sysSVM/"
CONFIG_PATH = "~/rosalind_lustre/mourikisa/data/OAC/config.R"
PREDICTION_SET = "~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/sysSVM/prediction_set.Rdata"
PREDICTION_SET_NS = "~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/sysSVM/prediction_set_noScale.Rdata"
STATS = "~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/best_models.tsv" ## the stats of best models
COLUMNS = "~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/sysSVM/OAC/config.txt"


source(CONFIG_PATH)

## -----------------------------
## Predict using the best models
## -----------------------------

load("~/rosalind_lustre/mourikisa/data/geneInfoNCG5.Rdata")
load("~/rosalind_lustre/mourikisa/data/cancerGenesNCG5.Rdata")

## Check point for NCG tissues
ncg_tissues = cancerGenes %>% select(primary_site) %>% unique %>% .$primary_site
if (!(TISSUE_NCG%in%ncg_tissues)){
    stop("Please provide a valid tissue for NCG5!")
}

cat("Getting prediction set...", "\n")
## Load training and validation sets
load(PREDICTION_SET)
load(PREDICTION_SET_NS)

## What features should be excluded
cat("Getting columns...", "\n")
config = read.table(COLUMNS, header = F)
config = config[10:nrow(config), 1]
s = c("cancer_type", "sample", "entrez", "no_ALL_muts","no_TRUNC_muts", "no_NTDam_muts",
      "no_GOF_muts", "Copy_number", "CNVLoss", "CNVGain")
s = s[!(s%in%config)]

prediction_set_ns = prediction_set_ns %>% tibble::rownames_to_column() %>% 
    separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>%
    mutate(entrez=as.numeric(entrez)) %>%
    select(one_of(s))

cat("Getting NCG5 data...", "\n")
## Fix gene info table from NCG
geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
## Get a cancer gene with all the associated primary sites and cancer sites
cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>% 
    group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","), 
                                   cancer_site=paste(unique(cancer_site), collapse=",")) %>%
    ungroup
geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))

best_models = read.table(STATS, header = T)

kernels = c("linear", "radial", "polynomial", "sigmoid")
for (k in kernels){
    cat(paste0("Predicting using kernel: ", k), "\n")
    
    if (k=="linear"){
        m = best_models$analysis[best_models$kernel=="linear"]
        load(paste0("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/model_", m, ".Rdata"))
    }else if (k=="radial"){
        m = best_models$analysis[best_models$kernel=="radial"]
        load(paste0("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/model_", m, ".Rdata"))
    }else if (k=="polynomial"){
        m = best_models$analysis[best_models$kernel=="polynomial"]
        load(paste0("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/model_", m, ".Rdata"))
    }else if (k=="sigmoid"){
        m = best_models$analysis[best_models$kernel=="sigmoid"]
        load(paste0("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/model_", m, ".Rdata"))
    }
    
    ## make prediction
    predictions = predict(svm.model, prediction_set, decision.values = TRUE)
    ## add probability from Platt and Binning
    predictions_platt = plattScaling(svm.model, predictions)
    colnames(predictions_platt) = c("dv", "label", "prob_platt")
    preds = predictions_platt %>% tibble::rownames_to_column()
    
    ## Take information from gene info
    preds = preds %>%
        separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>%
        mutate(entrez=as.numeric(entrez)) %>%
        left_join(geneInfo%>%rename(gene_type=cancer_type), by=c("entrez"))
    
    ## Join also with systems level properties from non-scaled set
    preds = preds %>% left_join(prediction_set_ns, by=c("cancer_type", "sample", "entrez"))
    
    save(preds, file=paste(RES_DIR,"/predictions_", k, ".Rdata", sep=""))
}

## -------------------
## Calculate the score
## -------------------

## Get the scores (based on the stats of best models)

cat("Calculating scores...", "\n")
LINEAR_CVS = round(best_models$mean[best_models$kernel=="linear"], digits = 2)
LINEAR_VAR = round(best_models$var[best_models$kernel=="linear"], digits = 2)

POLYNOMIAL_CVS = round(best_models$mean[best_models$kernel=="polynomial"], digits = 2)
POLYNOMIAL_VAR = round(best_models$var[best_models$kernel=="polynomial"], digits = 2)

RADIAL_CVS = round(best_models$mean[best_models$kernel=="radial"], digits = 2)
RADIAL_VAR = round(best_models$var[best_models$kernel=="radial"], digits = 2)

SIGMOID_CVS = round(best_models$mean[best_models$kernel=="sigmoid"], digits = 2)
SIGMOID_VAR = round(best_models$var[best_models$kernel=="sigmoid"], digits = 2)

scores = getScore(cancer_type=CANCER_TYPE, path=RES_DIR, 
                  linear_preds = "/predictions_linear.Rdata", linearCVS = LINEAR_CVS, linearVar = LINEAR_VAR,
                  radial_preds = "/predictions_radial.Rdata", radialCVS = RADIAL_CVS, radialVar = RADIAL_VAR,
                  polynomial_preds = "/predictions_polynomial.Rdata", polynomialCVS = POLYNOMIAL_CVS, polynomialVar = POLYNOMIAL_VAR,
                  sigmoid_preds = "/predictions_sigmoid.Rdata", sigmoidCVS = SIGMOID_CVS, sigmoidVar = SIGMOID_VAR)

save(scores, file=paste0(RES_DIR, "/scores.Rdata"))


## Sensitivity analysis
## Get gene rank per patient for all predictions
scores = scores[["scores"]]
load("~/rosalind_lustre/mourikisa/data/OAC/validation_cohort_OAC_TCGA/sysSVm_noPloidyCorrection/sysSVM/OAC/validation_set_noScale.Rdata")
load("~/rosalind_lustre/mourikisa/data/OAC/validation_cohort_OAC_TCGA/sysSVm_noPloidyCorrection/sysSVM/OAC/training_set_noScale.Rdata")
training_ns = training_ns%>%tibble::rownames_to_column()%>%
    separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>% 
    mutate(entrez=as.numeric(entrez)) %>% data.frame()
validation_ns = validation_ns%>%tibble::rownames_to_column()%>%
    separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")%>%
    mutate(type="P")%>%mutate(entrez=as.numeric(entrez))%>%data.frame()
cohort = rbind(training_ns, validation_ns)
cohort = cohort %>% left_join(scores)
load("~/rosalind_lustre/mourikisa/data/geneInfoNCG5.Rdata")
geneInfoNCG5 = geneInfo %>% select(entrez, symbol) %>% unique
cohort = cohort %>% left_join(geneInfoNCG5)
cohort = cohort %>% subset(!is.na(score))
cohort = cohort %>% group_by(sample) %>% mutate(patient.rank=row_number(desc(score))) %>% ungroup()
scores_plus_ranks = cohort
save(scores_plus_ranks, file=paste0(RES_DIR, "/scores_plus_ranks.Rdata"))


ggplot(sysSVM_predictions, aes(x=gene_type, y=score)) +
    geom_boxplot() +
    facet_wrap(~sample) +
    theme(
        strip.text.x = element_text(size = 5)
    ) +
    ylab("sysSVM score") +
    xlab("NCG5 gene category") +
    ggtitle("Comparison of genes scores in the validation set")


cgcs = training_ns %>% select(entrez) %>% unique %>% .$entrez

sysSVM_predictions %>% subset(patient.rank<=10 & entrez%in%cgcs) %>% select(symbol) %>% unique %>% nrow

t = sysSVM_predictions %>% select(-patient.rank) %>% 
    subset(CNVGain!=1) %>% group_by(sample) %>% mutate(patient.rank=row_number(desc(score))) %>% ungroup()

t %>% subset(patient.rank<=10 & entrez%in%cgcs) %>% select(symbol) %>% unique %>% nrow



getSysCans = function(path=NULL, cancer_type=NULL, tissue=NULL, scores_fn="scores.Rdata", 
                      fp_dir="~/rosalind_lustre/mourikisa/data/NCG_false_positives.txt",  
                      cgc_fn="~/rosalind_lustre/mourikisa/data/518_CGC_annotation_2211.xlsx",
                      previous.studies="~/rosalind_lustre/mourikisa/data/CGCs_to_be_considered_OAC.tsv",
                      exclude_expr=c("Not Expressed")){
    
    require(readxl)
    
    if(is.null(path) | is.null(cancer_type) | is.null(tissue)){
        stop("Missing arguments")
    }
    
    ## Get the data we need for the definition of the sys cans 
    ## Load geneInfo to get the symbols as well
    load("~/rosalind_lustre/mourikisa/data/geneInfoNCG5.Rdata")
    geneInfoNCG5 = geneInfo %>% select(entrez, symbol) %>% unique
    
    ## Load the CGC that we will consider drivers
    message("Getting CGC annotation...")
    #wb = loadWorkbook(cgc_fn) ## changed function
    #cgcs = readWorksheet(wb, sheet=1) ## changed funstion
    cgcs = read_excel(cgc_fn, 1)
    ## refine to those that will be retained
    retained_CGC = cgcs %>% subset(RET.DIS=="RET")
    retained_CGC_trans = retained_CGC %>% subset(grepl("Trans", Drivers.to.be.retained))
    ## Get training and validation set
    load(paste0(path, "/training_set_noScale.Rdata"))
    training_ns = training_ns%>%tibble::rownames_to_column()%>%
        separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>% 
        mutate(entrez=as.numeric(entrez)) %>% data.frame()
    ## Refine the training set to include only the ones we decided to coonsider drivers
    message("Refining CGC drivers...")
    cgc_drivers = NULL
    for (i in 1:nrow(retained_CGC)){
        en = retained_CGC$entrez[i]
        gene_type= retained_CGC$GENETICS[i]
        symbol = retained_CGC$symbol[i]
        retained_drivers = retained_CGC$Drivers.to.be.retained[i]
        #message(paste(symbol, gene_type, sep = "\t"))
        if (gene_type=="O"){
            if(grepl("Trans", retained_drivers)){
                dd = training_ns %>% subset(entrez==en) %>% subset(CNVGain==1 | no_NTDam_muts!=0 | no_GOF_muts!=0 | BND==1 | INV==1 | INS==1)
            }else{
                dd = training_ns %>% subset(entrez==en) %>% subset(CNVGain==1 | no_NTDam_muts!=0 | no_GOF_muts!=0) 
            }
            cgc_drivers = rbind(cgc_drivers, dd)
        }else if (gene_type=="TS"){
            if(grepl("Trans", retained_drivers)){
                dd = training_ns %>% subset(entrez==en) %>% subset(Copy_number==0 | no_TRUNC_muts!=0 | no_NTDam_muts!=0 | (Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)) | BND==1 | INV==1 | INS==1)
            }else{
                dd = training_ns %>% subset(entrez==en) %>% subset(Copy_number==0 | no_TRUNC_muts!=0 | no_NTDam_muts!=0 | (Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0))) 
            }
            cgc_drivers = rbind(cgc_drivers, dd)
        }else if (gene_type=="O/TS"){
            if(grepl("Trans", retained_drivers)){
                ## All drivers are considered
                dd = training_ns %>% subset(entrez==en)
            }else{
                ## Here we need to exlude all the genes with only SVs
                dd = training_ns %>% subset(entrez==en) %>% subset(!(no_TRUNC_muts==0 & no_NTDam_muts==0 & no_GOF_muts==0 & Copy_number!=0 & CNVGain==0 & (BND!=0 | INS!=0 | INV!=0)))
            }
            cgc_drivers = rbind(cgc_drivers, dd)
        }else{
            stop("Unknown gene type")
        }
    }
    
    ## Now add in the CGC drivers the CGCs identified from previous studies
    previous = read.delim(previous.studies, header = T, sep="\t")
    ## Get only those that were not considered already
    previous = previous %>% subset(to.be.considered & !already.considered)
    ## Get the entrez for these genes
    previous = previous %>% left_join(geneInfo%>%select(entrez, symbol), by=c("Gene.symbol"="symbol"))
    
    for(i in 1:nrow(previous)){
        en = previous$entrez[i]
        gene_type= previous$Type[i]
        symbol = previous$Gene.symbol[i]
        
        if(gene_type=="SVs"){
            dd = training_ns %>% subset(entrez==en) %>% subset(BND==1 | INV==1 | INS==1 | CNVLoss==1 | CNVGain==1)
        }else if(gene_type=="CNVs"){
            dd = training_ns %>% subset(entrez==en) %>% subset(CNVLoss==1 | CNVGain==1)
        }
        cgc_drivers = rbind(cgc_drivers, dd)
    }
    
    
    load(paste0(path, "/validation_set_noScale.Rdata"))
    validation_ns = validation_ns%>%tibble::rownames_to_column()%>%
        separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")%>%
        mutate(type="P")%>%mutate(entrez=as.numeric(entrez))%>%data.frame()
    cohort = rbind(cgc_drivers, validation_ns) ## Instead of the whole training, we use only the subset of the refined CGC 
    ## Get scores
    load(paste0(path, "/",scores_fn))
    scores = scores[["scores"]]
    
    cohort = cohort %>% left_join(scores)
    #cohort$gene_type[is.na(cohort$gene_type)] = "cgc" ## This is not true in the validation cohort, all genes get a score
    cohort = cohort %>% left_join(geneInfoNCG5)
    
    ## Make NAs in CGC score to be the highest score in each patient + 1 so we can sort by score in the output table - Same as above not true in the validation cohort
    #sample2maxscore = cohort %>% group_by(sample) %>% summarise(max_score=max(score, na.rm = T)+1)
    #cohort$score[is.na(cohort$score) & cohort$gene_type=="cgc"] = sample2maxscore$max_score[match(cohort$sample, sample2maxscore$sample)][which(is.na(cohort$score))]
    
    
    samples = unique(cohort$sample)
    sys_cans = NULL
    if (!is.null(tissue)){
        message("Expression is considered for the selection of sys cans...")
        ## Load expression data (GTeX)
        message("Loading GTeX expression data")
        load("~/athena/data/expressionNCG5_gtex.Rdata") ## name of data frame is gtex
        gtex = gtex %>% select(entrez, tissue, exp_level) %>% spread(tissue, exp_level)
        ## Checkpoint
        if(!(tissue%in%colnames(gtex))){
            stop("Wrong tissue supplied")
        }
        ## Select tissue of the cancer type
        gtex = gtex[,c("entrez", tissue)]
        colnames(gtex) = c("entrez", "gtex_tissue_expression")
        ## bring in the cohort table the expression
        message("Getting GTEx expression")
        cohort = cohort %>% left_join(gtex)
        
        
        total = length(samples)
        pb <- txtProgressBar(min = 0, max = total, style = 3)
        for (i in 1:length(samples)){
            Sys.sleep(0.1)
            ## When you get the candidates here you need to take into account expression in normal tissues and in the cancer sample
            
            ## Here I get the ranks for the non-CGC genes. The CGCs will be added as patient.rank==NA
            ## I also add 2 ranks, one with all the genes and one without amplifications
            d = cohort %>% subset(sample==samples[i] & !gtex_tissue_expression%in%exclude_expr & !is.na(gtex_tissue_expression)) %>% 
                arrange(desc(score)) %>% mutate(patient.rank=row_number(-score))
            d_withoutAmp = d %>% subset(sample==samples[i] & CNVGain==0 & !gtex_tissue_expression%in%exclude_expr & !is.na(gtex_tissue_expression)) %>% 
                arrange(desc(score)) %>% mutate(patient.rank.no.amp=row_number(-score))
            ## Bring in the rank without amplification
            d = d %>% left_join(d_withoutAmp)
            
            ## Add the CGCs - it is not needed for the validation cohort
            #d_cgcs = cohort %>% subset(sample==samples[i] & gene_type=="cgc" & !gtex_tissue_expression%in%exclude_expr & !is.na(gtex_tissue_expression))
            #d = rbind.fill(d, d_cgcs)
            
            ## Add all those that are not expressed (I will exclude those, I won't use them anyways) - code left here for testing
            #d_ne = cohort %>% subset(sample==samples[i] & (gtex_tissue_expression=="Not Expressed" | is.na(gtex_tissue_expression)))
            #d = rbind.fill(d, d_ne)
            
            ## Add them to the big table
            sys_cans = rbind(sys_cans, d)
            setTxtProgressBar(pb, i)
        }
    }
    
    ## Add false positive annotation
    false_positive_genes <- read.table(fp_dir, header = T, sep = "\t")
    false_positive_genes <- false_positive_genes %>% select(entrez) %>% .$entrez
    sys_cans = sys_cans %>% mutate(NCG_FP=ifelse(entrez%in%false_positive_genes, TRUE, FALSE))
    
    ## In the final output we need everything (also the CGC that were excluded) so create an export table -- No we don't, I disabled this. If needed in the future maybe turn it back on
    #export = rbind(training_ns, validation_ns)
    #export = export %>% left_join(cohort)
    ## all the rest CGC are now coming in and their gene type will be NA
    #export$gene_type[is.na(export$gene_type)] = "CGC_discarded"
    
    
    message(paste0("Score file used: ", scores_fn))
    
    return(sys_cans)
    
}



## Copied scores in sysSVM/OAC
syscan = getSysCans(path = "~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/sysSVM/OAC/", cancer_type = "OAC", tissue = "Esophagus")
syscan = getSysCans(path = "~/rosalind_lustre/mourikisa/data/OAC/validation_cohort_OAC_TCGA/sysSVm_noPloidyCorrection/sysSVM/OAC/", cancer_type = "OAC", tissue = "Esophagus")


## Plot the scores
syscan = syscan%>%mutate(type=ifelse(gene_type=="cgc", "cgc", "other"))
ggplot(syscan, aes(x=type, y=score)) +
    geom_boxplot() +
    theme(
        strip.text.x = element_text(size = 5)
    ) +
    ylab("sysSVM score") +
    xlab("") +
    theme_boss() +
    theme(
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")
    ) +
    ggtitle("Comparison of sysSVM genes scores in 21 OACs")


wilcox.test(syscan$score[syscan$type=="cgc"], syscan$score[syscan$type=="other"])



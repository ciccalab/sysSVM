rm(list=ls())

## Define command line arguments
args = commandArgs(trailingOnly = TRUE)

## Parameters
CANCER_TYPE = args[1]
TISSUE_NCG = args[2]
TISSUE_GTEX = args[3] ## Code commented out so it is not used at the moment
GOF=T
GAINS=T
if (args[4]=="F" | args[4]=="FALSE"){
    GAINS=F
}
CANCER_DIR = paste0(args[5], "/", CANCER_TYPE)
CONFIG_PATH = args[6]

ISOMAP=T
if (args[7]=="F"  | args[7]=="FALSE"){
    ISOMAP=F
}

source(CONFIG_PATH)

## --------------------------------------
## Get the results from the training step
## --------------------------------------

cat("Pooling metrics from CV iterations...", "\n")
df = gatherMetrics(paste0(CANCER_DIR, "/Results"))
cv_stats_summary <- df %>% subset(!is.na(value)) %>%
    group_by(analysis, kernel, nu, gamma, degree, type, set) %>% summarize(iterations=n(),
                                                                           min=min(value),
                                                                           q1=quantile(value)[2],
                                                                           median=median(value),
                                                                           mean=mean(value),
                                                                           q3=quantile(value)[4],
                                                                           max=max(value),
                                                                           var=stats::var(value)) %>% ungroup
write.table(cv_stats_summary, file=paste0(CANCER_DIR,"/cv_stats_summary.tsv"), quote=F, row.names=F, sep="\t")

## -------------------
## Get the best models
## -------------------
best_models = NULL
cat("Getting best models for:", "\n")
for (k in unique(cv_stats_summary$kernel)){
    cat(k, "\n")
    d = cv_stats_summary %>% subset(kernel==k & type=="Sensitivity" & set=="test") %>% data.frame()
    bm = d %>% arrange(desc(mean)) %>% slice(1:5) %>% mutate(var_rank=row_number(var)) %>% subset(var_rank==1)
    best_models = rbind(best_models, bm)
}
write.table(best_models, file=paste0(CANCER_DIR,"/best_models.tsv"), quote=F, row.names=F, sep="\t")

## For testing purposes
#cv_stats_summary = read.table(file=paste0(CANCER_DIR, "/cv_stats_summary.tsv"), header = T)
#best_models = read.table(file=paste0(CANCER_DIR,"/best_models.tsv"), header = T)

## -----------------------------
## Predict using the best models
## -----------------------------
load("/mnt/lustre/users/k1469280/mourikisa/data/geneInfoNCG5.Rdata")
load("/mnt/lustre/users/k1469280/mourikisa/data/cancerGenesNCG5.Rdata")

## Check point for NCG tissues
ncg_tissues = cancerGenes %>% select(primary_site) %>% unique %>% .$primary_site
if (!(TISSUE_NCG%in%ncg_tissues)){
    stop("Please provide a valid tissue for NCG5!")
}


## Load training and validation sets
load(paste(CANCER_DIR, "/training_set.Rdata", sep=""))
load(paste(CANCER_DIR, "/validation_set.Rdata", sep=""))
load(paste(CANCER_DIR, "/validation_set_noScale.Rdata", sep=""))

## What features should be excluded
config = read.table(paste(CANCER_DIR, "/config.txt", sep=""), header = F)
config = config[10:nrow(config), 1]
s = c("cancer_type", "sample", "entrez", "no_ALL_muts","no_TRUNC_muts", "no_NTDam_muts",
      "no_GOF_muts", "Copy_number", "CNVLoss", "CNVGain")
s = s[!(s%in%config)]

validation_ns = validation_ns %>% tibble::rownames_to_column() %>% 
    separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>%
    mutate(entrez=as.numeric(entrez)) %>%
    select(one_of(s))

## Fix gene info table from NCG
geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
## Get a cancer gene with all the associated primary sites and cancer sites
cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>% 
    group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","), 
                                   cancer_site=paste(unique(cancer_site), collapse=",")) %>%
    ungroup
geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))

## Training also contains "N" class? (just in case to be sure between the trials)
training = training%>%subset(type=="C")

## Filter out columns that are NaN - this should be changed to replace the column with all 0s
training=Filter(function(x)!all(is.nan(x)), training)
validation=Filter(function(x)!all(is.nan(x)), validation)


cat("Predicting using kernel:", "\n")
grid_predictions_table = NULL
grid_predictions = list()
novel_drivers_genes_platt_all = NULL
novel_drivers_genes_platt_inNCG5_all = NULL
novel_drivers_genes_platt_rst_all = NULL
new_samples_covered_all = NULL
for (i in 1:nrow(best_models)){
    kern = best_models$kernel[i]
    analysis = best_models$analysis[i]
    mynu=best_models$nu[i]
    mygamma=best_models$gamma[i]
    mydegree=best_models$degree[i]
    
    ## Train on the whole set
    cat(paste0("Kernel: ", kern), "\n")

    if (kern=="linear"){
        svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, nu=mynu)
    }else if (kern=="radial"){
        svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, gamma=mygamma, nu=mynu)
    }else if (kern=="polynomial"){
        svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, gamma=mygamma, nu=mynu, degree=mydegree)
    }else if (kern=="sigmoid"){
        svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, gamma=mygamma, nu=mynu)
    }
    
    
    ## Asses training sensitivity
    trainset = training %>% select(type)
    trainset$fitted = svm.model$fitted
    trainset = trainset %>% select(type, fitted)
    trainset = trainset %>% mutate(f=ifelse(fitted==TRUE, "C", "N"))
    trainset$type = factor(as.character(trainset$type), levels = c("C", "N"))
    trainset$f = factor(as.character(trainset$f), levels = c("C", "N"))
    ## Confusion matrix
    cM <- confusionMatrix(trainset$f,trainset$type, positive = c("C"))
    training_sensitivity = cM$byClass["Sensitivity"]
    
    ## make prediction
    predictions = predict(svm.model, validation, decision.values = TRUE)
    ## add probability from Platt and Binning
    predictions_platt = plattScaling(svm.model, predictions)
    colnames(predictions_platt) = c("dv", "label", "prob_platt")
    preds = predictions_platt %>% tibble::rownames_to_column()
    
    ## Take information from gene info
    preds = preds %>%
        separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>%
        mutate(entrez=as.numeric(entrez)) %>%
        left_join(geneInfo%>%rename(gene_type=cancer_type), by=c("entrez"))
    
    ## Join also with systems level properties from validaton set
    preds = preds %>% left_join(validation_ns, by=c("cancer_type", "sample", "entrez"))
    

    save(svm.model, file=paste(CANCER_DIR,"/model_",analysis, ".Rdata", sep=""))
    save(preds, file=paste(CANCER_DIR,"/predictions_",analysis, ".Rdata", sep=""))

    
#     ## Save the predicted entrez to make the venn diagram
#     grid_predictions[[analysis]] = preds%>%subset(label==TRUE & prob_platt >0.95)%>%select(entrez)%>%unique%>%.$entrez
#     
#     training_samples = training%>%tibble::rownames_to_column()%>%separate(rowname, into=c("Cancer_type", "sample", "entrez"), sep="\\.")%>%select(sample)%>%unique%>%.$sample
#     prediction_samples = preds%>%subset(label==TRUE & prob_platt >0.95)%>%select(sample)%>%unique%>%.$sample
#     
#     ## Cutoffs
#     cutoffs = c(0.95, 0.85, 0.75)
#     new_samples = preds%>%subset(label==TRUE & prob_platt >0.95)%>%select(sample)%>%subset(!(sample%in%training_samples))%>%unique%>%.$sample
#     new_samples_no_0.95 = length(new_samples)
#     new_samples = preds%>%subset(label==TRUE & prob_platt >0.85)%>%select(sample)%>%subset(!(sample%in%training_samples))%>%unique%>%.$sample
#     new_samples_no_0.85 = length(new_samples)
#     new_samples = preds%>%subset(label==TRUE & prob_platt >0.75)%>%select(sample)%>%subset(!(sample%in%training_samples))%>%unique%>%.$sample
#     new_samples_no_0.75 = length(new_samples)
#     
#     df = data.frame(analysis=analysis, nu=mynu,
#                     gamma=mygamma, degree=mydegree,
#                     training_set_size = nrow(training),
#                     training_genes = training%>%tibble::rownames_to_column()%>%separate(rowname, into=c("Cancer_type", "sample", "entrez"), sep="\\.")%>%select(entrez)%>%unique%>%nrow,
#                     training_samples = training%>%tibble::rownames_to_column()%>%separate(rowname, into=c("Cancer_type", "sample", "entrez"), sep="\\.")%>%select(sample)%>%unique%>%nrow,
#                     training_sensitivity = training_sensitivity,
#                     prediction_set_size=nrow(preds),
#                     prediction_genes=preds%>%select(entrez)%>%unique%>%nrow,
#                     prediction_samples=preds%>%select(sample)%>%unique%>%nrow,
#                     novel_drivers_entries=preds%>%subset(label==TRUE)%>%nrow,
#                     novel_drivers_genes=preds%>%subset(label==TRUE)%>%select(entrez)%>%unique%>%nrow,
#                     novel_drivers_genes_inNCG5 = preds%>%subset(label==TRUE & (gene_type=="cgc" | gene_type=="can"))%>%select(entrez)%>%unique%>%nrow,
#                     novel_drivers_genes_rst = preds%>%subset(label==TRUE & gene_type=="rst")%>%select(entrez)%>%unique%>%nrow,
#                     novel_drivers_entries_platt_0.95=preds%>%subset(label==TRUE & prob_platt >0.95)%>%nrow,
#                     novel_drivers_genes_platt_0.95=preds%>%subset(label==TRUE & prob_platt >0.95)%>%select(entrez)%>%unique%>%nrow,
#                     novel_drivers_in_samples_mutated_0.95 = preds%>%subset(label==TRUE & prob_platt >0.95)%>%select(sample)%>%unique%>%nrow,
#                     novel_drivers_genes_platt_inNCG5_0.95 =preds%>%subset(label==TRUE & prob_platt >0.95 & (gene_type=="cgc" | gene_type=="can"))%>%select(entrez)%>%unique%>%nrow,
#                     novel_drivers_genes_platt_rst_0.95=preds%>%subset(label==TRUE & prob_platt >0.95 & gene_type=="rst")%>%select(entrez)%>%unique%>%nrow,
#                     novel_drivers_genes_platt_rst_in_samples_mutated_0.95=preds%>%subset(label==TRUE & prob_platt >0.95 & gene_type=="rst")%>%select(sample)%>%unique%>%nrow,
#                     new_samples_covered_0.95 = new_samples_no_0.95,
#                     novel_drivers_entries_platt_0.85=preds%>%subset(label==TRUE & prob_platt >0.85)%>%nrow,
#                     novel_drivers_genes_platt_0.85=preds%>%subset(label==TRUE & prob_platt >0.85)%>%select(entrez)%>%unique%>%nrow,
#                     novel_drivers_in_samples_mutated_0.85 = preds%>%subset(label==TRUE & prob_platt >0.85)%>%select(sample)%>%unique%>%nrow,
#                     novel_drivers_genes_platt_inNCG5_0.85 =preds%>%subset(label==TRUE & prob_platt >0.85 & (gene_type=="cgc" | gene_type=="can"))%>%select(entrez)%>%unique%>%nrow,
#                     novel_drivers_genes_platt_rst_0.85=preds%>%subset(label==TRUE & prob_platt >0.85 & gene_type=="rst")%>%select(entrez)%>%unique%>%nrow,
#                     novel_drivers_genes_platt_rst_in_samples_mutated_0.85=preds%>%subset(label==TRUE & prob_platt >0.85 & gene_type=="rst")%>%select(sample)%>%unique%>%nrow,
#                     new_samples_covered_0.85 = new_samples_no_0.85,
#                     novel_drivers_entries_platt_0.75=preds%>%subset(label==TRUE & prob_platt >0.75)%>%nrow,
#                     novel_drivers_genes_platt_0.75=preds%>%subset(label==TRUE & prob_platt >0.75)%>%select(entrez)%>%unique%>%nrow,
#                     novel_drivers_in_samples_mutated_0.75 = preds%>%subset(label==TRUE & prob_platt >0.75)%>%select(sample)%>%unique%>%nrow,
#                     novel_drivers_genes_platt_inNCG5_0.75 =preds%>%subset(label==TRUE & prob_platt >0.75 & (gene_type=="cgc" | gene_type=="can"))%>%select(entrez)%>%unique%>%nrow,
#                     novel_drivers_genes_platt_rst_0.75=preds%>%subset(label==TRUE & prob_platt >0.75 & gene_type=="rst")%>%select(entrez)%>%unique%>%nrow,
#                     novel_drivers_genes_platt_rst_in_samples_mutated_0.75=preds%>%subset(label==TRUE & prob_platt >0.75 & gene_type=="rst")%>%select(sample)%>%unique%>%nrow,
#                     new_samples_covered_0.75 = new_samples_no_0.75
#     )
#     grid_predictions_table = rbind(grid_predictions_table, df)
#     
#     
#     for (cutoff in cutoffs){
#         cat(cutoff, "\n")
#         positive_genes_platt_rst=preds%>%subset(label==TRUE & prob_platt >cutoff & gene_type=="rst")%>%select(entrez)%>%unique%>%nrow
#         if (positive_genes_platt_rst==0){
#             next
#         }else{
#             p = plotSysproperties(preds, cutoff, tissue=TISSUE_NCG)
#             mylegend<-g_legend(p[[1]])
#             for(i in 1:length(p)) p[[i]] = p[[i]] + theme(legend.position="none")
#             pdf(file=paste(CANCER_DIR, "/predictions_sysproperties_",analysis,"_",cutoff, ".pdf" , sep=""), h= 20, w=24)
#             grid.arrange(p[[1]],p[[3]],p[[6]],p[[7]],p[[8]],
#                          p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]], 
#                          mylegend, ncol=3)
#             dev.off()
#             p = plotMolProperties(preds, gains=GAINS)
#             mylegend<-g_legend(p[[1]])
#             for(i in 1:length(p)) p[[i]] = p[[i]] + theme(legend.position="none")
#             pdf(file=paste(CANCER_DIR, "/predictions_molproperties_", analysis,"_",cutoff, ".pdf", sep=""), h= 20, w=24)
#             if(GAINS & GOF) {
#                 grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]], p[[6]], mylegend, ncol=3)
#             }else{
#                 grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]], mylegend, ncol=3)
#             }
#             dev.off()
#         }
#     }
}
# 
# require(xlsx)
# wb<-createWorkbook(type="xlsx")
# 
# TEXT_STYLE           <- CellStyle(wb) + Font(wb,  heightInPoints=14)
# TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
#     Alignment(rotation=45, horizontal="ALIGN_CENTER") +
#     Border(color="black", position=c("TOP", "LEFT","BOTTOM", "RIGHT"), pen=rep("BORDER_THIN",4)) 
# 
# 
# write_text_in_cell = function(sheet, rowIndex, text, textStyle){
#     rows <-createRow( sheet, rowIndex=rowIndex )
#     sheetText <-createCell( rows, colIndex=1 )
#     setCellValue( sheetText[[1,1]], text )
#     setCellStyle( sheetText[[1,1]], textStyle )
# }
# 
# sheet <- createSheet(wb, sheetName = "predictions")
# 
# text = with( grid_predictions_table, paste0("TRAINING SET. Number_of_genes = ", training_genes[1], " ; Samples = ", training_samples[1], " ; Total_entries = ", training_set_size[1]))
# write_text_in_cell(sheet, 1, text, TEXT_STYLE)
# 
# text = with( grid_predictions_table, paste0("TEST SET.     Number_of_genes = ", prediction_genes[1], " ; Samples = ", prediction_samples[1], " ; Total_entries = ", prediction_set_size[1]))
# write_text_in_cell(sheet, 2, text, TEXT_STYLE)
# 
# to_remove = c("training_genes","training_samples",'training_set_size','prediction_genes','prediction_samples','prediction_set_size')
# rownames(grid_predictions_table) = NULL
# addDataFrame(grid_predictions_table[, !colnames(grid_predictions_table)%in%to_remove ], sheet, startRow=7, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, row.names = F)
# saveWorkbook(wb,file=paste(CANCER_DIR,"/predictions.xlsx", sep=""))


## -------------------
## Calculate the score
## -------------------
 
## Get the scores (based on the best models above)
cat("Calculating scores...", "\n")
LINEAR_PREDS = best_models$analysis[best_models$kernel=="linear"]
LINEAR_CVS = round(best_models$mean[best_models$kernel=="linear"], digits = 2)
LINEAR_VAR = round(best_models$var[best_models$kernel=="linear"], digits = 2)

POLYNOMIAL_PREDS = best_models$analysis[best_models$kernel=="polynomial"]
POLYNOMIAL_CVS = round(best_models$mean[best_models$kernel=="polynomial"], digits = 2)
POLYNOMIAL_VAR = round(best_models$var[best_models$kernel=="polynomial"], digits = 2)

RADIAL_PREDS = best_models$analysis[best_models$kernel=="radial"]
RADIAL_CVS = round(best_models$mean[best_models$kernel=="radial"], digits = 2)
RADIAL_VAR = round(best_models$var[best_models$kernel=="radial"], digits = 2)

SIGMOID_PREDS = best_models$analysis[best_models$kernel=="sigmoid"]
SIGMOID_CVS = round(best_models$mean[best_models$kernel=="sigmoid"], digits = 2)
SIGMOID_VAR = round(best_models$var[best_models$kernel=="sigmoid"], digits = 2)

scores = getScore(cancer_type=CANCER_TYPE, path=CANCER_DIR, 
                   linear_preds = paste0("/predictions_", LINEAR_PREDS, ".Rdata"), linearCVS = LINEAR_CVS, linearVar = LINEAR_VAR,
                   radial_preds = paste0("/predictions_", RADIAL_PREDS, ".Rdata"), radialCVS = RADIAL_CVS, radialVar = RADIAL_VAR,
                   polynomial_preds = paste0("/predictions_", POLYNOMIAL_PREDS, ".Rdata"), polynomialCVS = POLYNOMIAL_CVS, polynomialVar = POLYNOMIAL_VAR,
                   sigmoid_preds = paste0("/predictions_", SIGMOID_PREDS, ".Rdata"), sigmoidCVS = SIGMOID_CVS, sigmoidVar = SIGMOID_VAR)
 
save(scores, file=paste0(CANCER_DIR, "/scores.Rdata"))


if (ISOMAP){
    ## Get geneInfo
    annotation = getNCGannotation()
     
    ## Get features and scores (for dimensionality reduction)
    load(paste0(CANCER_DIR, "/training_set.Rdata"))
    load(paste0(CANCER_DIR, "/validation_set.Rdata"))

    cohort = rbind(training%>%tibble::rownames_to_column()%>%separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>% data.frame(), 
                   validation%>%tibble::rownames_to_column()%>%separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")%>%mutate(type="P")%>%data.frame())
    
    feature_names = names(validation)
    cohort[,feature_names] = sapply(cohort[,feature_names], function(x) as.numeric(as.character(x))) ## as.character is used to maintain factor's levels 
    
    cohort = cohort %>% mutate(entrez=as.numeric(entrez)) %>% left_join(scores[["scores"]])
    cohort$gene_type[is.na(cohort$gene_type)] = "cgc" ## because cgc is not in the scores file
    cohort = cohort %>% left_join(annotation[["gtex"]]%>%subset(tissue==TISSUE_GTEX)%>%select(entrez, exp_level)%>%rename(tissue_exp_level=exp_level))
    
    #cohort = cohort %>% select(-one_of(nearZeroVar(cohort[,feature_names], name=T)))
    
    ## Exclude features with near zero variance
    #cat("Features excluded from dimensionality reduction:", "\n")
    #print(nearZeroVar(cohort[,4:37], name=T))
     
    
    ## ------
    ## Isomap
    ## ------
    cat("Running Isomap...", "\n")
    iso_dim2 = Isomap(data=as.matrix(cohort[,feature_names]), dims=2)
    iso_dim2 = iso_dim2 %>% data.frame()
    colnames(iso_dim2) = c("X1", "X2")
    
    
    ## --------------------
    ## Make the final table
    ## --------------------
    cat("Exporting final table...", "\n")
    ## Get features and scores (for exporting - no scaling in the features)
    load(paste0(CANCER_DIR, "/training_set_noScale.Rdata"))
    load(paste0(CANCER_DIR, "/validation_set_noScale.Rdata"))
    cohort = rbind(training_ns%>%tibble::rownames_to_column()%>%separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>% data.frame(), 
                   validation_ns%>%tibble::rownames_to_column()%>%separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")%>%mutate(type="P")%>%data.frame())
    cohort = cohort %>% mutate(entrez=as.numeric(entrez)) %>% left_join(scores[["scores"]])
    cohort$gene_type[is.na(cohort$gene_type)] = "cgc"
    cohort = cohort %>% left_join(annotation[["gtex"]]%>%subset(tissue==TISSUE_GTEX)%>%select(entrez, exp_level)%>%rename(tissue_exp_level=exp_level))
    
    ## Molecular properties
    if (GAINS){
        mol_properties = c("no_ALL_muts", "no_TRUNC_muts", "no_NTDam_muts", "no_GOF_muts", "Copy_number", "CNVGain" ,"CNVLoss", "ExpT_ME", "ExpT_HE", "ExpT_LE", "ExpT_NE", "ExpT_NET")
    }else{
        mol_properties = c("no_ALL_muts", "no_TRUNC_muts", "no_NTDam_muts", "no_GOF_muts", "Copy_number", "CNVLoss", "ExpT_ME", "ExpT_HE", "ExpT_LE", "ExpT_NE", "ExpT_NET")
    }
    
    
    ## Annotate and export the table
    export = cohort %>% select(cancer_type, sample, entrez, type, gene_type, N, Nlog10, sum_of_scores, kernels_used, score, kernels_predicted, kernels_predicted_no ) %>% left_join(annotation[["geneInfo"]])
    export = export %>% left_join(cohort%>%select(cancer_type, sample, entrez, one_of(mol_properties))) %>% 
        left_join(cohort%>%select(-one_of(mol_properties), -type, -gene_type, -N, -Nlog10, -sum_of_scores, -kernels_used, -score, -kernels_predicted, -kernels_predicted_no))
    
    ## Add false positive annotation
    fp_dir="/home/mourikisa/data/NCG_false_positives.txt"
    false_positive_genes <- read.table(fp_dir, header = T, sep = "\t")
    false_positive_genes <- false_positive_genes %>% select(entrez) %>% .$entrez
    export = export %>% mutate(in_false_positives=ifelse(entrez%in%false_positive_genes, TRUE, FALSE))
    
    ## Add the coordinates
    export = cbind(export, iso_dim2)
    
    write.xlsx(export%>%data.frame(), file=paste0(CANCER_DIR, "/scores_plus_annotation.xlsx"), row.names = F)
    message("Isomap run finished")
}else{
    message("Finished without Isomap")
}




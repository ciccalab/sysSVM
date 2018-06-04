###############################################################################
##  Predicting cancer genes on the validation cohort using sysSVM best models
###############################################################################

detach(name = "package:dplyr", unload = T)
library(plyr)
library(dplyr)

ns = c("nonsynonymous","stopgain","frameshift deletion","splicing","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","stoploss","frameshift substitution")
dam = c("nonsynonymous","frameshift deletion","frameshift insertion","frameshift substitution","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","splicing","stopgain","stoploss")
trunc = c("frameshift deletion","frameshift insertion","frameshift substitution","stopgain","stoploss") ## Always damaging==TRUE
non_trunc = c("nonsynonymous","splicing")

## First load the data
## muts
load("~/rosalind_lustre/mourikisa/data/TCGA/01_03_2015/ESCA/Tumor/Somatic_Mutations/ESCA_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF_Annovar_dbNSFP_gAnn_SampleFilters_NonSilent_MutationFilters_OncodriveClust.Rdata")
muts = do.call(rbind, somatic_mutations)
rownames(muts) = NULL ## funny rownames from subsettings etc

## CNVs
load("~/rosalind_lustre/mourikisa/data/TCGA/01_03_2015/ESCA/Tumor/CNV/ESCA_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_GainLoss_BLAT_overlaps_segmentsUnique.Rdata")
cnvs = do.call(rbind, somatic_cnv)
rownames(cnvs) = NULL

## Find out which of them are OAC
clinical = read.table("~/rosalind_lustre/mourikisa/data/TCGA/01_03_2015/ESCA/Clinical/Biotab/nationwidechildrens.org_clinical_patient_esca.txt", header = T, sep = "\t", skip = 1)
clinical = clinical %>% slice(2:nrow(clinical))
samples = muts %>% select(Patient) %>% unique %>% left_join(clinical%>%select(bcr_patient_barcode, histological_type)%>%rename(Patient=bcr_patient_barcode))
samples_OAC = samples %>% subset(histological_type=="Esophagus Adenocarcinoma, NOS")

## From those I can only include samples for which I have ploidy
## Load ploidy data
load("~/rosalind_lustre/mourikisa/data/OAC/validation_cohort_OAC_TCGA/ploidy_TCGA_from_COSMIC.Rdata")
ploidy_TCGA = ploidy_TCGA %>% subset(Patient%in%samples_OAC$Patient)


## subset mutations and cnvs for OAC samples
muts = muts %>% subset(Patient%in%ploidy_TCGA$Patient)
cnvs = cnvs %>% subset(Patient%in%ploidy_TCGA$Patient) %>% left_join(ploidy_TCGA%>%select(Patient, Ploidy))

## Fix CNVs using ploidy
cnvs %>% mutate(CNV_type_corrected=ifelse(Copy_number>=2*Ploidy, "Gain", ifelse(Copy_number<2, "Loss", NA))) %>% subset(CNV_type_corrected=="Gain") %>% select(Sample)  %>% unique %>% nrow

## Create the iput for ML
fix_splicing=function(x){
    x[which(x$Func.refGene=="splicing"),"ExonicFunc.refGene"]="splicing"
    x
}

fix_exonic=function(x){
    x[which(x$Func.refGene=="exonic;splicing"),"Func.refGene"]="exonic"
    x
}

fix_exonicFunc=function(x){
    if(length(which(x$ExonicFunc.refGene=="nonsynonymous SNV"))>0) x[which(x$ExonicFunc.refGene=="nonsynonymous SNV"),"ExonicFunc.refGene"]="nonsynonymous"
    if(length(which(x$ExonicFunc.refGene=="synonymous SNV"))>0) x[which(x$ExonicFunc.refGene=="synonymous SNV"),"ExonicFunc.refGene"]="synonymous"
    if(length(which(x$ExonicFunc.refGene=="stopgain SNV"))>0) x[which(x$ExonicFunc.refGene=="stopgain SNV"),"ExonicFunc.refGene"]="stopgain"
    if(length(which(x$ExonicFunc.refGene=="stoploss SNV"))>0) x[which(x$ExonicFunc.refGene=="stoploss SNV"),"ExonicFunc.refGene"]="stoploss"
    x
}

is_nonsilent=function(x){
    x$nonsilent=x$ExonicFunc.refGene%in%ns
    x
}

set_IGV_code = function(x){
    x$IGV = paste(x$chr,":",x$end,sep="",coll="")
    x
}

## This function has been altered from the original version in functions.R because for TCGA samples we dont have any SV data
createTotalTable = function(muts=NULL, cnvs=NULL, exclude_samples=NULL){
    
    if(is.null(muts) | is.null(cnvs) ){
        stop("Missing data")
    }
    
    ## Make the lists
    message("Integrating SNVs...")
    df_mut = muts
    rm(muts)
    
    if(!is.null(exclude_samples)){
        df_mut = df_mut %>% subset(!sample%in%exclude_samples)
        message(paste0("Samples excluded in mutation data: ", paste0(exclude_samples, collapse=",")))
    }
    
    ## Fix nonsilent here
    df_mut=fix_splicing(df_mut)
    df_mut=fix_exonic(df_mut)
    df_mut=fix_exonicFunc(df_mut)
    df_mut=is_nonsilent(df_mut)
    df_mut=set_IGV_code(df_mut)
    df_mut=is_nonsilent(df_mut)
    
    ## In order to get the number of all mutations per gene and because
    ## I have WGS data, I exclude mutations that fall in the following categories
    df_mut = df_mut %>% subset(Func.refGene!="" &
                                   !grepl("downstream", df_mut$Func.refGene) &
                                   !grepl("upstream", df_mut$Func.refGene) &
                                   !grepl("intergenic", df_mut$Func.refGene) &
                                   !grepl("ncRNA", df_mut$Func.refGene) &
                                   !grepl("intronic", df_mut$Func.refGene) &    
                                   !grepl("UTR", df_mut$Func.refGene))
    
    ## Exclude genes that are not in 19014
    df_mut = df_mut %>% rename(sample=Sample)
    df_mut = df_mut %>% mutate(entrez_19014=as.numeric(entrez_19014)) ## in case it is not numeric
    df_mut = df_mut %>% subset(!is.na(entrez_19014))
    
    ## Create the total table
    total_muts = ddply(df_mut, .(sample, symbol_19014, entrez_19014), summarise,
                       no_ALL_muts=n(),
                       no_NSI_muts=sum(nonsilent),
                       no_TRUNC_muts = sum(ExonicFunc.refGene %in% trunc),
                       no_NTDam_muts = sum(ExonicFunc.refGene %in% non_trunc & damaging),
                       no_GOF_muts = sum(oncodriveClust), .progress = 'text'
    )
    
    ## Bring in the total_table the SVs
    message("Integrating SVs...")
    
    ## Simply add SVs here as 0s (this needs to be done automatically for missing features when the toll will be realeased on its own)
    total_table = total_muts %>% mutate(BND=0, INS=0, INV=0)
    total_table = total_table %>% mutate(key=paste(sample, entrez_19014, sep="."))
    
    
    
    message("Integrating CNVs...")
    
    ## Add also the genes that are in muts and SVs to get their ploidy and Copy number
    cnvs = cnvs %>% rename(entrez_19014=entrez, sample=Sample)
    cnvs = cnvs %>% mutate(key=paste(sample, entrez_19014, sep="."))
    df_cnv = cnvs %>% subset(key%in%total_table$key | !is.na(CNV_type)) ## Also get the real CNVs
    
    
    if(!is.null(exclude_samples)){
        df_cnv = df_cnv %>% subset(!sample%in%exclude_samples)
        message(paste0("Samples excluded in CNV data: ", paste0(exclude_samples, collapse=",")))
    }
    
    ## define Gains and Losses - this was done in previous step
    ## Deduplicate the CNV data, because some genes may fall into two regions (sometimes it can be gain and loss)
    message("Resolving duplicated entries in CNVs...")
    dups = df_cnv %>% group_by(key) %>% mutate(n=n(), types=paste(unique(CNV_type), collapse=","), ntypes=length(unique(CNV_type))) %>% subset(n>1)
    dups = dups %>% ungroup()
    
    if(nrow(dups)>0){
        ## In here you will find two kinds of duplications those that are duplicates but associated with one type of CNV (i.e Gain/Loss)
        ## And those that are associated with two types of CNVs
        ## I didn't use overlap function in the end
        overlap <- function(start1, end1, start2, end2){
            res = pmin(end1, end2) - pmax(start2, start1)
            if(res>=0){
                return(res) 
            }else if(res<0){
                res=0
                return(res) 
            } 
        }
        dups_refined = NULL
        for (t in unique(dups$types)){
            if (grepl(",NA|NA,", t)){ ## When arrange always on top will be Gain/Loss and those will be selected when deduplicate
                d = dups %>% subset(types==t) %>% arrange(key, desc(CNV_type)) %>% subset(!duplicated(key))
                dups_refined = rbind(dups_refined, d)
            }else if (t=="Gain" | t=="NA" | t=="Loss"){ ## Choose the one with the highest overlap
                d = dups %>% subset(types==t) %>% mutate(overlap=(overlap(start, end, Start, End)/(end-start))*100) %>% arrange(key, desc(overlap)) %>% subset(!duplicated(key)) %>% select(-overlap)
                dups_refined = rbind(dups_refined, d)
            }else if (grepl("Gain", t) & grepl("Loss", t)){ ## For those genes we have both gain and loss, I set CNV_type to NA because we cannot distinguish between the two
                d = dups %>% subset(types==t) %>% mutate(CNV_type=NA) %>% arrange(key) %>% subset(!duplicated(key))
                dups_refined = rbind(dups_refined, d)
            }
        }
        
        ## For now deduplicte them and keep as CNV_type_corrected the concatenation of both types to see how many they are in the drivers
        ## Take them out first from the df_cnv
        df_cnv = df_cnv %>% subset(!key%in%dups$key)
        df_cnv$n = 1
        ## Fix dups
        dups_refined = dups_refined %>% select(-types, -ntypes)
        ## Put the back in the df_cnv
        df_cnv = rbind(df_cnv, dups_refined)
        
        ## Create total table from mutations and CNVs
        total_table = total_table %>% subset(!is.na(entrez_19014)) %>% 
            full_join(df_cnv%>%select(sample, entrez_19014, Copy_number, CNV_type, n)%>%rename(CNV_entries=n)%>%subset(!is.na(entrez_19014)))
    }else{
        ## Create total table from mutations and CNVs
        total_table = total_table %>% subset(!is.na(entrez_19014)) %>% 
            full_join(df_cnv%>%select(sample, entrez_19014, Copy_number, CNV_type)%>%subset(!is.na(entrez_19014)))
    }
    
    total_table$na_19014 = apply(total_table[,c("symbol_19014", "entrez_19014")], 1, function(x) length(x[is.na(x)]))
    total_table = total_table %>% mutate(in_19014=ifelse(na_19014<2, TRUE, FALSE)) %>% select(-na_19014) %>% subset(in_19014==TRUE)
    
    return(total_table)
    
}
total_table = createTotalTable(muts=muts, cnvs = cnvs)

## Check the following line and then delete it
## df_mut %>% subset(symbol_19014=="TP53" & (damaging | oncodriveClust) & Sample=="TCGA-L5-A8NW-01") %>% data.frame()

save(total_table, file="~/rosalind_lustre/mourikisa/data/OAC/validation_cohort_OAC_TCGA/total_table.Rdata")
## Give to this function the output of the createTotalTable function above
getMLinput <- function(df, ct="OAC", geneProperties_dir="~/rosalind_lustre/mourikisa/data/geneProperties_final_mmImputed.Rdata"){
    
    df = df %>% mutate(Cancer_type=ct) %>% rename(Sample=sample, Entrez=entrez_19014) %>% select(-symbol_19014, -in_19014, -key)
    
    ## Replace numbers with names here
    ## We assume every gene with no mutation data that it's not mutated
    message("Fixing mutations...")
    df[,c("no_ALL_muts", "no_NSI_muts", "no_TRUNC_muts", 
          "no_NTDam_muts",
          "no_GOF_muts")][is.na(df[,c("no_ALL_muts", "no_NSI_muts", "no_TRUNC_muts",
                                      "no_NTDam_muts",
                                      "no_GOF_muts")])] <- 0
    
    message("Fixing CNVs...")
    ## Copy number (where copy number is NA, put copy number equal to 2)
    ## I integrated copy number data by selecting segment mean >|0.3| therefore I took only gains and losses
    ## But at the same time the unique number of genes in the CNV data is quite high, therefore whatever is left with NA is probably 2
    df$Copy_number[is.na(df$Copy_number)] <- 2
    
    message("Fixing SVs...") ## No structural variation data are available
    df$BND[is.na(df$BND)] = 0
    df$INS[is.na(df$INS)] = 0
    df$INV[is.na(df$INV)] = 0
    
    ## Join with gene Properties
    message("Joining table with systems-level properties...")
    load(geneProperties_dir)
    geneProperties = geneProperties_mmImputed
    geneProperties = geneProperties %>% select(-symbol, -cancer_type, -cancer_dom, -cancer_rec)
    
    df <- df %>% left_join(geneProperties, by=c("Entrez"))
    
    ## Convert categorical features to multiple factors
    message("Performing cleaning of categorical variables...")
    ## CNV type
    df <- df %>% 
        mutate(CNVGain=ifelse(is.na(CNV_type), 0, ifelse(CNV_type=="Gain",1, 0)), 
               CNVLoss=ifelse(Copy_number==0 | Copy_number==1, 1, 0)) %>%
        select(-CNV_type)
    
    ## age
    df <- df %>% 
        mutate(old=ifelse(is.na(age), NA, ifelse(age=="old",1, 0)),
               young=ifelse(is.na(age), NA, ifelse(age=="young",1, 0))) %>%
        select(-age)
    
    ## origin
    df <- df %>% 
        mutate(luca=ifelse(is.na(origin), NA, ifelse(origin=="LUCA",1, 0)), 
               eukaryotes=ifelse(is.na(origin), NA, ifelse(origin=="Eukaryotes",1, 0)),
               metazoans=ifelse(is.na(origin), NA, ifelse(origin=="Metazoans",1, 0)),
               vertebrates=ifelse(is.na(origin), NA, ifelse(origin=="Vertebrates",1, 0)),
               opisthokonts=ifelse(is.na(origin), NA, ifelse(origin=="Opisthokonts",1, 0)),
               mammals=ifelse(is.na(origin), NA, ifelse(origin=="Mammals",1, 0)),
               primates=ifelse(is.na(origin), NA, ifelse(origin=="Primates", 1, 0))) %>% 
        select(-origin)
    
    ## exp.breadth.class
    df <- df %>% 
        mutate(selective=ifelse(is.na(exp.breadth), NA, ifelse(exp.breadth=="Selective",1, 0)), 
               always.expressed=ifelse(is.na(exp.breadth), NA, ifelse(exp.breadth=="AlwaysExpressed",1, 0)),
               middle=ifelse(is.na(exp.breadth), NA, ifelse(exp.breadth=="Middle",1, 0)),
               one.tissue=ifelse(is.na(exp.breadth), NA, ifelse(exp.breadth=="OneTissue",1, 0)),
               never.expressed=ifelse(is.na(exp.breadth), NA, ifelse(exp.breadth=="Neverexpressed",1, 0))) %>%
        select(-exp.breadth)
    
    
    
    ## Before you add (change NAs in these columns to 0)
    #df[,c("High", "Low", "Medium", "NotExpressed")][is.na(df[,c("High", "Low", "Medium", "NotExpressed")])] <- 0
    #df <- df %>% ungroup %>% mutate(tot.tissues=High+Low+Medium)
    df <- data.frame(df)
    
    message("Converting features to factors...")
    fcols <- c("duplicated",
               "WGD", "hub", "central", "CNVGain", "CNVLoss",
               "ExpT_ME", "ExpT_HE", "ExpT_LE", "ExpT_NE",
               "ExpT_NET", "old", "young", "luca", "eukaryotes",
               "metazoans", "vertebrates", "opisthokonts",
               "mammals", "primates", "selective", "always.expressed",
               "middle", "one.tissue", "never.expressed")
    cols <- which(colnames(df) %in% fcols)
    for(i in cols){
        df[,i] = factor(df[,i], levels = c(0,1))
    }
    
    ## Reorder columns
    df = df[,c(12, 1:2, 3:11, 13:length(df))]
    
    return(df)
}
oac_TCGA_ML_input = getMLinput(df=total_table)
write.table(oac_TCGA_ML_input, file="~/rosalind_lustre/mourikisa/data/OAC/validation_cohort_OAC_TCGA/oac_TCGA_ML_input.tsv", row.names = F, sep = "\t", quote=F)


## And now prepare the table for the sysSVM
## OK - to speed up this step I run the first part of the pipeline which is reparation of training and prediction set.
## I submitted the job to rosalind using prepare.sh
## I will then merge them and do the prediction on everything
load("~/rosalind_lustre/mourikisa/data/OAC/validation_cohort_OAC_TCGA/sysSVM/OAC/training_set.Rdata")
load("~/rosalind_lustre/mourikisa/data/OAC/validation_cohort_OAC_TCGA/sysSVM/OAC/validation_set.Rdata")
prediction_set = rbind(training[,!names(training)%in%c("type")], validation)
save(prediction_set, file="~/rosalind_lustre/mourikisa/data/OAC/validation_cohort_OAC_TCGA/prediction_set.Rdata")

load("~/rosalind_lustre/mourikisa/data/OAC/validation_cohort_OAC_TCGA/sysSVM/OAC/training_set_noScale.Rdata")
load("~/rosalind_lustre/mourikisa/data/OAC/validation_cohort_OAC_TCGA/sysSVM/OAC/validation_set_noScale.Rdata")
prediction_set_ns = rbind(training_ns[,!names(training_ns)%in%c("type")], validation_ns)
save(prediction_set_ns, file="~/rosalind_lustre/mourikisa/data/OAC/validation_cohort_OAC_TCGA/prediction_set_noScale.Rdata")







#!/usr/bin/Rscript

## **************************************
##           mutation annotation
## **************************************

setwd("/mnt/lustre/users/k1469280/mourikisa/data/OAC")
cat("Sourcing relevant functions...")
source("functions.R")


## Command line arguments
SAMPLE_DIR = commandArgs(trailingOnly = TRUE)[1]
SNVs = commandArgs(trailingOnly = TRUE)[2]
INDELs = commandArgs(trailingOnly = TRUE)[3]
CNVs = commandArgs(trailingOnly = TRUE)[4]
CNV_stats = commandArgs(trailingOnly = TRUE)[5]
#SVs = commandArgs(trailingOnly = TRUE)[6]
GENOME = commandArgs(trailingOnly = TRUE)[6]

MUT_OUT = paste0(SAMPLE_DIR, "/parsing_and_annotation")
if(!dir.exists(MUT_OUT)){
    dir.create(MUT_OUT)
}

# Create output files if they dont exist
if(!dir.exists(paste0(MUT_OUT, "/Rdata"))){
    dir.create(paste0(MUT_OUT, "/Rdata"))
}

if(!dir.exists(paste0(MUT_OUT, "/annovar"))){
    dir.create(paste0(MUT_OUT, "/annovar"))
}


CNV_OUT = paste0(gsub("strelka", "ascat", SAMPLE_DIR), "/parsing_and_annotation")
if(!dir.exists(CNV_OUT)){
    dir.create(CNV_OUT)
}

# SV_OUT = paste0(gsub("strelka", "manta", SAMPLE_DIR), "/parsing_and_annotation")
# if(!dir.exists(SV_OUT)){
#     dir.create(SV_OUT)
# }

## Get SNVs
# cat("Parsing SNVs")
# snv = getSNVs(SNVs, GENOME)
# save(snv, file=paste0(MUT_OUT, "/Rdata/snvs.Rdata"))
# indel = getIndels(INDELs, GENOME)
# indel = unrowname(indel)
# save(indel, file=paste0(MUT_OUT, "/Rdata/indels.Rdata"))

## Run ANNOVAR to annotate SNVs and Indels
# muts = annotateMutations(snv=snv, indel=indel, g=GENOME, save_dir=paste0(MUT_OUT, "/annovar"))
# save(muts, file=paste0(MUT_OUT, "/annovar/muts_annovar_dbnsfp_19014.Rdata"))
# 
# ns = c("nonsynonymous","stopgain","frameshift deletion","splicing","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","stoploss","frameshift substitution")
# dam = c("nonsynonymous","frameshift deletion","frameshift insertion","frameshift substitution","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","splicing","stopgain","stoploss")
# trunc = c("frameshift deletion","frameshift insertion","frameshift substitution","stopgain","stoploss") ## Always damaging==TRUE
# non_trunc = c("nonsynonymous","splicing")
# ns_vep=c("missense_variant", "splice_region_variant", "splice_donor_variant", "stop_gained", "splice_acceptor_variant", "stop_lost")

## Get CNVs
cat("Parsing CNVs")
cnvs = getCNVs(cnv_fn=CNVs, stats_fn=CNV_stats, g=GENOME, save_dir=CNV_OUT)
## define Gains and Losses
cnvs[["df_cnvs_19014"]] = cnvs[["df_cnvs_19014"]] %>% mutate(CNV_type_corrected=ifelse(Total_CN>=2*ploidy, "Gain", ifelse(Total_CN<2, "Loss", NA)))
save(cnvs, file=paste0(CNV_OUT, "/cnvs.Rdata"))

## SVs
# svs = getSVs(sv_fn=SVs)
# save(svs, file=paste0(SV_OUT, "/svs.Rdata"))

# ## Plot overall
# toPlot = svs[,1:7] %>% gather(type, value, -Sample, -gene) %>% subset(value!=0) %>% group_by(Sample, type) %>% summarise(n=n()) %>% ungroup()
# p = ggplot(toPlot, aes(x=Sample, y=n)) + geom_bar(stat="identity") + 
#     geom_hline(data=toPlot %>% group_by(type) %>% summarise(m=median(n, na.rm=T)), aes(yintercept=m, color="red", linetype="dashed")) + 
#     facet_wrap(~type, nrow = 5, scales = "free") + 
#     theme(axis.text.x=element_blank()) + ylab("Genes (#)")
# pdf(file="129_OAC/Tables_and_Plots/SV_distributions.pdf", w=10,h=10)
# print(p)
# dev.off()
# 
# ## Get distributions
# svs[,c(1:7, 9)] %>% subset(!is.na(entrez_19014)) %>% select(-entrez_19014) %>% gather(type, value, -Sample, -gene) %>% subset(value!=0) %>% 
#     group_by(Sample, type) %>% summarise(n=n()) %>% ungroup() %>% group_by(type) %>% summarise(min=min(n), median=median(n), mean=mean(n), max=max(n), n=n())
# 
# svs[,c(1:7, 10)] %>% subset(cancer_type=="cgc" | cancer_type=="can") %>% select(-cancer_type) %>% gather(type, value, -Sample, -gene) %>% subset(value!=0) %>% 
#     group_by(Sample, type) %>% summarise(n=n()) %>% ungroup() %>% group_by(type) %>% summarise(min=min(n), median=median(n), mean=mean(n), max=max(n), n=n())






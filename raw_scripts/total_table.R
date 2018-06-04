source("config.R")
source("Thanos_TCGA/config.R")

load("Rdata/dataset_mutations_damaging.Rdata")
load("Rdata/dataset_cnvs.Rdata")
load("Rdata/dataset_rna.Rdata")

ns = c("nonsynonymous","stopgain","frameshift deletion","splicing","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","stoploss","frameshift substitution")
dam = c("nonsynonymous","frameshift deletion","frameshift insertion","frameshift substitution","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","splicing","stopgain","stoploss")
trunc = c("frameshift deletion","frameshift insertion","frameshift substitution","stopgain","stoploss") ## Always damaging==TRUE
non_trunc = c("nonsynonymous","splicing")

df_mut = subset(df_mut, !is.na(symbol_19014) & nonsilent)

df_cnv = df_cnv_19014
# GISTIC OPTIONS!!
df_cnv = subset(df_cnv_19014, gistic)

rm(df_cnv_19014, tpms, counts, counts_bygene)



total_muts = ddply(df_mut, .(sample, symbol_19014, entrez_19014), summarise,
                   no_NSI_muts=sum(nonsilent),
                   no_TRUNC_muts = sum(ExonicFunc.refGene %in% trunc),
                   no_NTDam_muts = sum(ExonicFunc.refGene %in% non_trunc & damaging),
                   no_NTDamFunction_muts = sum(damagingFunc),
                   no_NTDamCons_muts = sum(damagingCons),
                   no_NTDamSC_muts = sum(damagingSC),
                   no_GOF_muts = sum(oncodriveClust), .progress = 'text'
                   )

total_cnvs =  subset(df_cnv, Total_CN>=4 | Total_CN<2) [, c('sample', 'symbol_19014', 'entrez_19014', 'CNV_type' ,  'Total_CN')]
colnames(total_cnvs)[5] = 'Copy_number'

total_rna = df_rna[, c('sample', 'symbol_19014', 'expr_class','TPM')]

symbol = sort(unique(c(total_muts$symbol_19014, total_cnvs$symbol_19014, total_rna$symbol_19014)))
symbol = symbol[which(symbol!="-")]
sample = unique(df_cnv$sample)


total_muts$key = paster(total_muts$sample,".",total_muts$symbol_19014)
total_cnvs$key = paster(total_cnvs$sample,".",total_cnvs$symbol_19014)
total_rna$key = paster(total_rna$sample,".",total_rna$symbol_19014)

eac     = data.frame(sample=rep(sample, each=length(symbol)), symbol_19014 = symbol)
eac$key = paster(eac$sample,".",eac$symbol_19014)
eac = cbind(eac, 
            total_muts[ match(eac$key, total_muts$key), which(colnames(total_muts)%nin%c('sample', 'symbol_19014', 'entrez_19014', 'key'))],
            total_cnvs[ match(eac$key, total_cnvs$key), which(colnames(total_cnvs)%nin%c('sample', 'symbol_19014', 'entrez_19014', 'key'))], 
            total_rna[ match(eac$key, total_rna$key), which(colnames(total_rna)%nin%c('sample', 'symbol_19014', 'key'))])

eac$expr_class = gsub(" ","_", eac$expr_class)

# adding gene properties
load("Rdata/geneProperties.Rdata")


## Load gene information to annotate gene as dominant/recessive/candidate/rest

dev_mode()
# install.packages("~/Lavoro/ncglib")
library(ncglib)
ncg_connection = get_ncg_connection()
gene_info = get_geneinfo()
cancer_genes = subset(get_cancer_genes(), primary_site=="esophagus")
dev_mode()

write.xlsx(cancer_genes, file="EAC_cancer_genes.1202.xlsx", "Cancer genes", showNA = F)

driver = unique(cancer_genes$symbol)
eac = cbind( eac, gene_info[ match(eac$symbol_19014, gene_info$symbol), c("entrez","cancer_type","cancer_dom","cancer_rec")  ])
eac$cgc_esophagous = eac$symbol_19014%in%driver

vogel <- read.table("/Volumes/ceredam/Thanos/vogelstein_genes.csv", header = T, sep = " ")
row.names(vogel) <- NULL
vogel$genegroup = gsub("Vog.","",vogel$genegroup)
eac$vogel = vogel[ match(eac$entrez, vogel[,1]), 2]

eac = cbind(eac, geneProperties[match(eac$entrez, geneProperties$Entrez), 2:ncol(geneProperties)])
eac = eac[,c(1,2, 15:ncol(eac),3:14)]
colnames(eac)[3] = "Entrez"

eac.table=eac

# cat("Fixing copy number variations according to new cut-offs...", "\n")
# eac <- eac %>% mutate(CNV_type=ifelse(Copy_number>=4, "Gain", ifelse(Copy_number<2, "Loss", NA)))

## Replace numbers with names here
eac.table[,c("no_NSI_muts", "no_TRUNC_muts", "no_NTDam_muts", "no_NTDamFunction_muts",
      "no_NTDamCons_muts", "no_NTDamSC_muts","no_GOF_muts")][
        is.na(eac.table[,c("no_NSI_muts", "no_TRUNC_muts", "no_NTDam_muts", "no_NTDamFunction_muts", "no_NTDamCons_muts", "no_NTDamSC_muts","no_GOF_muts")])] <- 0


## Convert categorical features to multiple factors
## CNV type
eac.table <- eac.table %>% 
  mutate(CNVGain=ifelse(is.na(CNV_type), 0, ifelse(CNV_type=="Gain",1, 0)), 
         CNVLoss=ifelse(is.na(CNV_type), 0, ifelse(CNV_type=="Loss",1, 0))) 
# %>%    select(-CNV_type)

## Copy number (where copy number is NA, put copy number equal to 2)
eac.table$Copy_number[which(is.na(eac.table$Copy_number))] <- 2

## Convert the length to bp instead of aa
eac.table <- eac.table %>% mutate(Length.fullrefseq=Length.fullrefseq*3)

## expr_class
eac.table <- eac.table %>% 
  mutate(ExpT_ME=ifelse(is.na(expr_class), 0, ifelse(expr_class=="ME",1, 0)), 
         ExpT_HE=ifelse(is.na(expr_class), 0, ifelse(expr_class=="HE",1, 0)),
         ExpT_LE=ifelse(is.na(expr_class), 0, ifelse(expr_class=="LE",1, 0)),
         ExpT_NE=ifelse(is.na(expr_class), 0, ifelse(expr_class=="not_expr",1, 0)),
         ExpT_NET=ifelse(is.na(expr_class), 0, ifelse(expr_class=="not_expr_NT",1, 0)))
#%>%  select(-expr_class)

## age
eac.table <- eac.table %>% 
  mutate(old=ifelse(is.na(age), 0, ifelse(age=="old",1, 0)),
         young=ifelse(is.na(age), 0, ifelse(age=="young",1, 0)))
# %>%  select(-age)

## origin
eac.table <- eac.table %>% 
  mutate(luca=ifelse(is.na(origin), 0, ifelse(origin=="LUCA",1, 0)), 
         eukaryotes=ifelse(is.na(origin), 0, ifelse(origin=="Eukaryotes",1, 0)),
         metazoans=ifelse(is.na(origin), 0, ifelse(origin=="Metazoans",1, 0)),
         vertebrates=ifelse(is.na(origin), 0, ifelse(origin=="Vertebrates",1, 0)),
         opisthokonts=ifelse(is.na(origin), 0, ifelse(origin=="Opisthokonts",1, 0)),
         mammals=ifelse(is.na(origin), 0, ifelse(origin=="Mammals",1, 0))) 
# %>%  select(-origin)

## exp.breadth.class
eac.table <- eac.table %>% 
  mutate(selective=ifelse(is.na(exp.breadth), 0, ifelse(exp.breadth=="Selective",1, 0)), 
         always.expressed=ifelse(is.na(exp.breadth), 0, ifelse(exp.breadth=="AlwaysExpressed",1, 0)),
         middle=ifelse(is.na(exp.breadth), 0, ifelse(exp.breadth=="Middle",1, 0)),
         one.tissue=ifelse(is.na(exp.breadth), 0, ifelse(exp.breadth=="OneTissue",1, 0)),
         never.expressed=ifelse(is.na(exp.breadth), 0, ifelse(exp.breadth=="Neverexpressed",1, 0))) 
# %>%  select(-exp.breadth)

## Degree, Betweenness, Hub, Central have NAs
eac.table <- eac.table %>% mutate(degree=ifelse(is.na(degree), 0, degree), betweenness=ifelse(is.na(betweenness), 0, betweenness),
                    hub=ifelse(is.na(hub), 0, hub), central=ifelse(is.na(central), 0, central))

## Before you add (chnage NAs in these columns to 0)
eac.table[,c("High", "Low", "Medium", "NotExpressed")][is.na(eac.table[,c("High", "Low", "Medium", "NotExpressed")])] <- 0
eac.table <- eac.table %>% ungroup %>% mutate(tot.tissues=High+Low+Medium)
eac.table <- data.frame(eac.table)

fcols <- c("Genic", "memberofcomplex",
           "WGD", "hub", "central", "CNVGain", "CNVLoss",
           "ExpT_ME", "ExpT_HE", "ExpT_LE", "ExpT_NE",
           "ExpT_NET", "old", "young", "luca", "eukaryotes",
           "metazoans", "vertebrates", "opisthokonts",
           "mammals", "selective", "always.expressed",
           "middle", "one.tissue", "never.expressed")
cols <- which(colnames(eac.table) %in% fcols)
for(i in cols){
  eac.table[,i] = factor(eac.table[,i])
}

eac.table = eac.table[, which(!colnames(eac.table)%in%c("CNV_type","expr_class","age","origin","exp.breadth"))]

# eac.table <- eac.table %>% mutate(key=paste(Cancer_type, Sample, Entrez, sep=".")) %>% 
#   select(-Cancer_type, -Sample, -Entrez)
# rnames <- eac.table$key
# eac.table <- eac.table %>% select(-key)
# eac.table <- data.frame(eac.table)
# row.names(eac.table) <- rnames

save(eac, eac.table, file="Rdata/EAC.12_02_16.Rdata")

save(eac, eac.table, file="Rdata/EAC.gistic.16_02_16.Rdata")



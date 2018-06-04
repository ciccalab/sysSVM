########################################
##      Check expression data
########################################

rm(list=ls())

#library('biomaRt')
library(gplots)

setwd("~/rosalind_lustre/mourikisa/data/OAC/")

## Get the data for the new samples
fns = list.files("rnaseq_data/FPKM_Juliane_03_12_17/counts/", pattern = "txt")


mappings_16 = read.table("rnaseq_data/old/16_samples/wgs_rnaseq_mapping.txt", header = F, sep="\t")
colnames(mappings_16) = c("id1", "sample", "id")
mappings_16 = mappings_16 %>% dplyr::select(sample, id)
mappings_19 = read.table("rnaseq_data/old/19_samples/rna_seq_mapping.txt", header = F, sep="\t")
colnames(mappings_19) = c("id1", "sample", "id")
mappings_19 = mappings_19 %>% mutate(id2=gsub("(\\.[^\\.]*)$", "", id)) %>% 
    mutate(id2=gsub("(\\.[^\\.]*)$", "", id2)) %>% 
    dplyr::select(-id) %>% rename(id=id2) %>% 
    dplyr::select(sample, id)
mappings = rbind(mappings_16, mappings_19)

fpkms = NULL
for(fn in fns){
    d = read.table(paste0("rnaseq_data/FPKM_Juliane_03_12_17/counts/", fn), header = T, sep="\t")
    d = d %>% tibble::rownames_to_column() %>% mutate(fn=fn) %>% rename(gene=rowname)
    s = unlist(strsplit(fn, "_primaryReads"))[1]
    d = d %>% mutate(id=s)
    d = d %>% left_join(mappings)
    fpkms = rbind(fpkms, d)
}

## Exclude sample LP6005690-DNA_F02 because this was sequenced twice - the one excluded sequenced in the pilot of ICGC
fpkms = fpkms %>% subset(id!="SLX-10113.D705_D502")

load("~/rosalind_lustre/mourikisa/data/geneInfoNCG5_2.Rdata")
geneInfo = geneInfo %>% subset(duplicability!=-1) %>% dplyr::select(entrez, symbol, refdna) %>% 
    mutate(refdna=strsplit(refdna, ";")) %>% unnest()

fpkms = fpkms %>% rename(refdna=gene) %>% left_join(geneInfo)
save(fpkms, file="rnaseq_data/FPKM_Juliane_03_12_17/fpkms.Rdata")

# ## To check if I can reduce the lost NMs via biomart
# ## convert NMs ids to entrez ids
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# genes <- unique(fpkms$refdna)
# glist <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna", "entrezgene", "hgnc_symbol"),values=genes,mart= mart)
# glist = glist %>% rename(refdna=refseq_mrna, entrez_biomart=entrezgene, hgnc_symbol_biomart=hgnc_symbol)
# detach("package:biomaRt", unload=TRUE)
# ## Sometime I run the script from this bit that's why the next couple of lines are repeated
# ## I will intersect the list from biomart with NCG
# load("~/rosalind_lustre/mourikisa/data/geneInfoNCG5_2.Rdata")
# geneInfo = geneInfo %>% subset(duplicability!=-1) %>% dplyr::select(entrez, symbol, refdna) %>% 
#     mutate(refdna=strsplit(refdna, ";")) %>% unnest()
# glist = glist %>% subset(entrez_biomart%in%geneInfo$entrez)
# 
# fpkms = fpkms %>% left_join(glist)
# save(fpkms, file="~/rosalind_lustre/mourikisa/data/OAC/rnaseq_data/fpkms_plus_biomart_annotation.Rdata")

## I collapse different mRNAs and take the max fpkm
expression = fpkms %>% subset(!is.na(entrez)) %>% group_by(sample, entrez, symbol) %>% summarise(R.FPKM_union=max(R.FPKM_union), R.FPKM_intersectionNotEmpt=max(R.FPKM_intersectionNotEmpt)) %>% ungroup()


load("~/rosalind_lustre/mourikisa/data/geneInfoNCG5_2.Rdata")
load("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/training_set_noScale.Rdata")
load("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/validation_set_noScale.Rdata")
cohort = rbind(training_ns%>% tibble::rownames_to_column(), validation_ns%>% tibble::rownames_to_column()%>%mutate(type="P"))
cohort = cohort %>% separate(rowname, into=c("ct", "sample", "entrez"), sep="\\.") %>% 
    separate(sample, into=c("sample", "normal"), sep="_vs_") %>% 
    mutate(entrez=as.numeric(entrez))
cohort = cohort %>% left_join(geneInfo%>%subset(duplicability!=-1)%>%dplyr::select(entrez, symbol, cancer_type)%>%rename(gene_type=cancer_type))

## Check how many samples in our cohort
venn(list(rnaseq=unique(fpkms$sample), samples_161=unique(cohort$sample)))
## Select only samples in our cohort
expression = expression %>% subset(sample%in%cohort$sample)

## Amplified helpers
load("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/patient_stratification_onlyTop10.Rdata")
dt = d
dt = dt %>% separate(sample, into=c("sample", "normal"), sep="_vs_")
a = dt %>% subset(CNVGain==1) %>% dplyr::select(sample, entrez, symbol) %>% unique %>% mutate(type="Amplification") %>% mutate(key=paste0(sample, ".", entrez))
## Plots for amplifications (I decided to go with entrez ids) & get rid of some duplicated entrez due to the conversion
toPlot = expression %>% subset(entrez%in%a$entrez) %>% mutate(key=paste0(sample, ".", entrez))
genesWithAtLeastOneAmp = toPlot %>% subset(key%in%a$key) %>% .$entrez %>% unique ## we do that because we want only those that are amplified in at least one sample
toPlot = toPlot %>% subset(entrez%in%genesWithAtLeastOneAmp) %>% mutate(type=ifelse(key %in% a$key, "Amplification", "WT"))

toPlot = toPlot %>% left_join(toPlot %>% group_by(type) %>% summarise(n=n())) %>% mutate(label=paste0(type, "(n=", n, ")"))
ggplot(toPlot, aes(x=label, y=R.FPKM_union)) + geom_violin() + geom_boxplot(width=0.1) +
    scale_y_log10(breaks=c(-10, 0, 1, 3, 5, 10, 25, 50, seq(min(toPlot$R.FPKM_union), max(toPlot$R.FPKM_union), 100)))+ 
    #theme_boss() + 
    theme(
        axis.title.x = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")
    ) +
    ylab("log10(FPKM)") + xlab("") + ggtitle("Expression comparison for 202 amplified helper genes in 32 samples")
write.table(toPlot, file="~/Desktop/amplified_vs_WT_helpers.tsv", row.names = F, quote = F, sep = "\t")
print(wilcox.test(toPlot$FPKM[toPlot$type=="Amplification"], toPlot$FPKM[toPlot$type=="WT"]))
print(summary(toPlot$FPKM[toPlot$type=="Amplification"]))
print(summary(toPlot$FPKM[toPlot$type=="WT"]))

print(wilcox.test(toPlot$R.FPKM_union[toPlot$type=="Amplification"], toPlot$R.FPKM_union[toPlot$type=="WT"]))
print(summary(toPlot$R.FPKM_union[toPlot$type=="Amplification"]))
print(summary(toPlot$R.FPKM_union[toPlot$type=="WT"]))


print(wilcox.test(toPlot$R.FPKM_intersectionNotEmpt[toPlot$type=="Amplification"], toPlot$R.FPKM_intersectionNotEmpt[toPlot$type=="WT"]))
print(summary(toPlot$R.FPKM_intersectionNotEmpt[toPlot$type=="Amplification"]))
print(summary(toPlot$R.FPKM_intersectionNotEmpt[toPlot$type=="WT"]))


## Amplified drivers
load("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/patient_stratification_onlyCGC.Rdata")
dc = d
dc = dc %>% separate(sample, into=c("sample", "normal"), sep="_vs_")
a = dc %>% subset(CNVGain==1) %>% dplyr::select(sample, entrez, symbol) %>% unique %>% mutate(type="Amplification") %>% mutate(key=paste0(sample, ".", entrez))
toPlot = expression %>% subset(entrez%in%a$entrez) %>% mutate(key=paste0(sample, ".", entrez))
genesWithAtLeastOneAmp = toPlot %>% subset(key%in%a$key) %>% .$entrez %>% unique ## we do that because we want only those that are amplified in at least one sample
toPlot = toPlot %>% subset(entrez%in%genesWithAtLeastOneAmp) %>% mutate(type=ifelse(key %in% a$key, "Amplification", "WT"))

toPlot = toPlot %>% left_join(toPlot %>% group_by(type) %>% summarise(n=n())) %>% mutate(label=paste0(type, "(n=", n, ")"))

print(wilcox.test(toPlot$R.FPKM_union[toPlot$type=="Amplification"], toPlot$R.FPKM_union[toPlot$type=="WT"]))
print(summary(toPlot$R.FPKM_union[toPlot$type=="Amplification"]))
print(summary(toPlot$R.FPKM_union[toPlot$type=="WT"]))


print(wilcox.test(toPlot$R.FPKM_intersectionNotEmpt[toPlot$type=="Amplification"], toPlot$R.FPKM_intersectionNotEmpt[toPlot$type=="WT"]))
print(summary(toPlot$R.FPKM_intersectionNotEmpt[toPlot$type=="Amplification"]))
print(summary(toPlot$R.FPKM_intersectionNotEmpt[toPlot$type=="WT"]))


## Amplified genes overall
t = cohort
a = t %>% subset(CNVGain==1) %>% dplyr::select(sample, entrez, symbol) %>% unique %>% mutate(type="Amplification") %>% mutate(key=paste0(sample, ".", entrez))
toPlot = expression %>% subset(entrez%in%a$entrez) %>% mutate(key=paste0(sample, ".", entrez))
genesWithAtLeastOneAmp = toPlot %>% subset(key%in%a$key) %>% .$entrez %>% unique ## we do that because we want only those that are amplified in at least one sample
toPlot = toPlot %>% subset(entrez%in%genesWithAtLeastOneAmp) %>% mutate(type=ifelse(key %in% a$key, "Amplification", "WT"))

toPlot = toPlot %>% left_join(toPlot %>% group_by(type) %>% summarise(n=n())) %>% mutate(label=paste0(type, "(n=", n, ")"))

print(wilcox.test(toPlot$R.FPKM_union[toPlot$type=="Amplification"], toPlot$R.FPKM_union[toPlot$type=="WT"]))
print(summary(toPlot$R.FPKM_union[toPlot$type=="Amplification"]))
print(summary(toPlot$R.FPKM_union[toPlot$type=="WT"]))


print(wilcox.test(toPlot$R.FPKM_intersectionNotEmpt[toPlot$type=="Amplification"], toPlot$R.FPKM_intersectionNotEmpt[toPlot$type=="WT"]))
print(summary(toPlot$R.FPKM_intersectionNotEmpt[toPlot$type=="Amplification"]))
print(summary(toPlot$R.FPKM_intersectionNotEmpt[toPlot$type=="WT"]))



##############
## Old code
##############

fns = list.files("rnaseq_data/16_samples/")
fns = fns[!grepl("mapping", fns)]
rpkms_16 = data.frame()
for(fn in fns){
    d = read.table(paste0("rnaseq_data/16_samples/", fn), header = T, sep="\t")
    d = d %>% tibble::rownames_to_column() %>% mutate(id=gsub("_rpkm.txt", "", fn)) %>% rename(gene=rowname)
    rpkms_16 = rbind(rpkms_16, d)
}

## Get the mappings
mappings_16 = read.table("rnaseq_data/16_samples/wgs_rnaseq_mapping.txt", header = F, sep="\t")
colnames(mappings_16) = c("id1", "sample", "id")
mappings_16 = mappings_16 %>% dplyr::select(sample, id)
rpkms_16 = rpkms_16 %>% left_join(mappings_16)

###########################################################
##      C = Number of reads mapped to a gene
##      N = Total mapped reads in the experiment
##      L = exon length in base-pairs for a gene
##      Equation: RPKM = (10^9 * C)/(N * L)
###########################################################
## Get RPKM
samples2N = rpkms_16 %>% group_by(id) %>% summarise(N=sum(readcounts))
rpkms_16 = rpkms_16 %>% left_join(samples2N)
rpkms_16 = rpkms_16 %>% mutate(rpkm=exp( log(readcounts) + log(1e9) - log(genelength) - log(N) ))

## Get the data for the old samples
fns = list.files("rnaseq_data/19_samples/")
fns = fns[!grepl("mapping", fns)]
rpkms_19 = data.frame()
for(fn in fns){
    d = read.table(paste0("rnaseq_data/19_samples/", fn), header = T, sep="\t")
    d = d %>% tibble::rownames_to_column() %>% mutate(id=gsub("_rpkm.txt", "", fn)) %>% rename(gene=rowname)
    rpkms_19 = rbind(rpkms_19, d)
}

## Get the mappings
mappings_19 = read.table("rnaseq_data/19_samples/rna_seq_mapping.txt", header = F, sep="\t")
colnames(mappings_19) = c("id1", "sample", "id")
mappings_19 = mappings_19 %>% mutate(id2=gsub("(\\.[^\\.]*)$", "", id)) %>% 
    mutate(id2=gsub("(\\.[^\\.]*)$", "", id2)) %>% 
    dplyr::select(-id) %>% rename(id=id2) %>% 
    dplyr::select(sample, id)
rpkms_19 = rpkms_19 %>% left_join(mappings_19)

## Exclude sample LP6005690-DNA_F02 because this was sequenced twice
rpkms_19 = rpkms_19 %>% subset(sample!="LP6005690-DNA_F02")

## Get RPKM also for this one
samples2N = rpkms_19 %>% group_by(id) %>% summarise(N=sum(readcounts))
rpkms_19 = rpkms_19 %>% left_join(samples2N)
rpkms_19 = rpkms_19 %>% mutate(rpkm=exp( log(readcounts) + log(1e9) - log(genelength) - log(N) ))


## Combine the two tables
rpkms = rbind(rpkms_16, rpkms_19 %>% dplyr::select(-R.FPKM))

## convert ensembl ids to entrez ids
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rpkms$gene
genes = gsub("\\..*","",genes)
genes = unique(genes)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),values=genes,mart= mart)
G_list = G_list %>% rename(gene=ensembl_gene_id)
G_list = G_list %>% subset(!duplicated(gene)) ## some ensembl ids return more than one symbols (they are only 10 anyways)

rpkms$gene = gsub("\\..*","",rpkms$gene)
rpkms = rpkms %>% left_join(G_list)
save(rpkms, file="~/athena/data/OAC/rnaseq_data/rpkm.Rdata")

load("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/training_set_noScale.Rdata")
load("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/validation_set_noScale.Rdata")
cohort = rbind(training_ns%>% tibble::rownames_to_column(), validation_ns%>% tibble::rownames_to_column()%>%mutate(type="P"))
cohort = cohort %>% separate(rowname, into=c("ct", "sample", "entrez"), sep="\\.") %>% 
    separate(sample, into=c("sample", "normal"), sep="_vs_") %>% 
    mutate(entrez=as.numeric(entrez))
load("~/rosalind_lustre/mourikisa/data/geneInfoNCG5_2.Rdata")
cohort = cohort %>% left_join(geneInfo%>%subset(duplicability!=-1)%>%dplyr::select(entrez, symbol, cancer_type)%>%rename(gene_type=cancer_type))

## Check how many samples in our cohort
venn(list(rnaseq=unique(rpkms$sample), samples_161=unique(cohort$sample)))

rpkms = rpkms %>% subset(sample%in%cohort$sample)
## Exclude genes that cannot be mapped to entrez/symbol
rpkms = rpkms %>% subset(!(is.na(entrezgene) & is.na(hgnc_symbol)))

## Amplified helpers
load("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.top10.plusCGC.containers.corrected.Rdata")
t = genes2paths %>% subset(gene_type!="cgc")
t = t %>% separate(sample, into=c("sample", "normal"), sep="_vs_")
a = t %>% subset(CNVGain==1) %>% dplyr::select(sample, entrez, symbol) %>% unique %>% mutate(type="Amplification") %>% mutate(key=paste0(sample, ".", entrez))
## Plots for amplifications (I decided to go with entrez ids) & get rid of some duplicated entrez due to the conversion
toPlot = rpkms %>% subset(entrezgene%in%a$entrez) %>% mutate(key=paste0(sample, ".", entrezgene)) %>% arrange(sample, desc(rpkm)) %>% subset(!duplicated(key))
genesWithAtLeastOneAmp = toPlot %>% subset(key%in%a$key) %>% .$entrezgene %>% unique ## we do that because we want only those that are amplified in at least one sample
toPlot = toPlot %>% subset(entrezgene%in%genesWithAtLeastOneAmp) %>% mutate(type=ifelse(key %in% a$key, "Amplification", "WT"))

toPlot = toPlot %>% left_join(toPlot %>% group_by(type) %>% summarise(n=n())) %>% mutate(label=paste0(type, "(n=", n, ")"))
ggplot(toPlot, aes(x=label, y=rpkm)) + geom_violin() + geom_boxplot(width=0.1) + 
    scale_y_log10(breaks=c(0, 1, 3, 5, 10, 25, 50, seq(min(toPlot$rpkm), max(toPlot$rpkm), 100)))+ 
    theme_boss() + 
    theme(
        axis.title.x = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")
    ) +
    ylab("Expression (RPKM)") + xlab("") + ggtitle("Expression comparison for 209 amplified helper genes in 32 samples \n (top 10 plus CGC)")
write.table(toPlot, file="~/Desktop/amplified_helpers.tsv", row.names = F, quote = F, sep = "\t")
print(wilcox.test(toPlot$rpkm[toPlot$type=="Amplification"], toPlot$rpkm[toPlot$type=="WT"]))
print(wilcox.test(toPlot$rpkm[toPlot$type=="Amplification"], toPlot$rpkm[toPlot$type=="WT"], alternative = "greater"))
print(summary(toPlot$rpkm[toPlot$type=="Amplification"]))
print(summary(toPlot$rpkm[toPlot$type=="WT"]))


## Amplified drivers
load("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.top10.plusCGC.containers.corrected.Rdata")
t = genes2paths %>% subset(gene_type=="cgc")
t = t %>% separate(sample, into=c("sample", "normal"), sep="_vs_")
a = t %>% subset(CNVGain==1) %>% dplyr::select(sample, entrez, symbol) %>% unique %>% mutate(type="Amplification") %>% mutate(key=paste0(sample, ".", entrez))
toPlot = rpkms %>% subset(entrezgene%in%a$entrez) %>% mutate(key=paste0(sample, ".", entrezgene)) %>% arrange(sample, desc(rpkm)) %>% subset(!duplicated(key))
genesWithAtLeastOneAmp = toPlot %>% subset(key%in%a$key) %>% .$entrezgene %>% unique ## we do that because we want only those that are amplified in at least one sample
toPlot = toPlot %>% subset(entrezgene%in%genesWithAtLeastOneAmp) %>% mutate(type=ifelse(key %in% a$key, "Amplification", "WT"))

toPlot = toPlot %>% left_join(toPlot %>% group_by(type) %>% summarise(n=n())) %>% mutate(label=paste0(type, "(n=", n, ")"))
ggplot(toPlot, aes(x=label, y=rpkm)) + geom_violin() + geom_boxplot(width=0.1) + 
    scale_y_log10(breaks=c(0, 1, 3, 5, 10, 25, 50, seq(min(toPlot$rpkm), max(toPlot$rpkm), 100)))+ 
    theme_boss() + 
    theme(
        axis.title.x = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")
    ) +
    ylab("Expression (RPKM)") + xlab("") + ggtitle("Expression comparison for 52 amplified driver genes (top 10 plus CGC)")
write.table(toPlot, file="~/Desktop/amplified_drivers.tsv", row.names = F, quote = F, sep = "\t")
print(wilcox.test(toPlot$rpkm[toPlot$type=="Amplification"], toPlot$rpkm[toPlot$type=="WT"]))
print(wilcox.test(toPlot$rpkm[toPlot$type=="Amplification"], toPlot$rpkm[toPlot$type=="WT"], alternative = "greater"))
print(summary(toPlot$rpkm[toPlot$type=="Amplification"]))
print(summary(toPlot$rpkm[toPlot$type=="WT"]))


## Amplified genes in prediction set
t = cohort %>% subset(gene_type!="cgc")
a = t %>% subset(CNVGain==1) %>% dplyr::select(sample, entrez, symbol) %>% unique %>% mutate(type="Amplification") %>% mutate(key=paste0(sample, ".", entrez))
toPlot = rpkms %>% subset(entrezgene%in%a$entrez) %>% mutate(key=paste0(sample, ".", entrezgene)) %>% arrange(sample, desc(rpkm)) %>% subset(!duplicated(key))
genesWithAtLeastOneAmp = toPlot %>% subset(key%in%a$key) %>% .$entrezgene %>% unique ## we do that because we want only those that are amplified in at least one sample
toPlot = toPlot %>% subset(entrezgene%in%genesWithAtLeastOneAmp) %>% mutate(type=ifelse(key %in% a$key, "Amplification", "WT"))

toPlot = toPlot %>% left_join(toPlot %>% group_by(type) %>% summarise(n=n())) %>% mutate(label=paste0(type, "(n=", n, ")"))
ggplot(toPlot, aes(x=label, y=rpkm)) + geom_violin() + geom_boxplot(width=0.1) + 
    scale_y_log10(breaks=c(0, 1, 3, 5, 10, 25, 50, 100, 500, 1000, 2000, seq(min(toPlot$rpkm), max(toPlot$rpkm), 3000)))+ 
    theme_boss() + 
    theme(
        axis.title.x = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")
    ) +
    ylab("Expression (RPKM)") + xlab("") + ggtitle("Expression comparison for 5,039 amplified genes in prediction set")
write.table(toPlot, file="~/Desktop/amplified_prediction_set.tsv", row.names = F, quote = F, sep = "\t")
print(wilcox.test(toPlot$rpkm[toPlot$type=="Amplification"], toPlot$rpkm[toPlot$type=="WT"]))
print(wilcox.test(toPlot$rpkm[toPlot$type=="Amplification"], toPlot$rpkm[toPlot$type=="WT"], alternative = "greater"))
print(summary(toPlot$rpkm[toPlot$type=="Amplification"]))
print(summary(toPlot$rpkm[toPlot$type=="WT"]))


## Check individual genes

## E2F
load("~/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.top10.plusCGC.containers.corrected.Rdata")
t = genes2paths %>% subset(grepl("E2F", symbol))
t = t %>% separate(sample, into=c("sample", "normal"), sep="_vs_")
rpkms %>% subset(sample%in%t$sample) %>% select(sample) %>% unique
rpkms %>% subset(entrezgene%in%t$entrez) %>% select(entrezgene) %>% unique %>% nrow

genes = c("E2F1", "E2F2", "E2F3", "E2F4", "E2F5" )
plots = list()
for(g in genes){
  s_amp = t %>% subset(symbol==g) %>% .$sample %>% unique
  length(s_amp)
  expr = rpkms %>% subset(entrezgene==unique(t$entrez[t$symbol==g]))
  expr = expr %>% mutate(type=ifelse(sample%in%s_amp, "Amplification", "WT"))
  
  expr = expr %>% left_join(expr %>% count(type)) %>% mutate(label=paste0(type, "(n=", n, ")"))
  
  p = ggplot(expr, aes(x=type, y=rpkm)) + geom_boxplot(width=0.1) +
    geom_point() +
    theme_boss() + 
    theme(
      axis.title.x = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")
    ) +
    ylab(paste0(g, " expression (RPKM)")) + xlab("")
  
  plots[[g]] = p
  
}

grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], ncol=3)


## Bring Copy number in to select the gene
toPlot = toPlot %>% left_join(syscans%>%dplyr::select(sample, symbol, CNVGain, Copy_number)%>%unique)
## Bring also whether they have 
sigma = read_xlsx("~/Desktop/amplified_genes_no_cgc_sigma.xlsx", 1) %>% rename(symbol=Gene)

toPlot = toPlot %>% left_join(sigma)
write.table(toPlot, file="~/Desktop/amplifications_vs_WT_syscans.tsv", quote = F, row.names = F, sep="\t")


## Plots for Homozygous deletions
homs = syscans %>% subset(Copy_number==0 & symbol%in%h) %>% dplyr::select(sample, symbol) %>% unique %>% mutate(type="HomDel")
toPlot = rpkms %>% subset(hgnc_symbol%in%h) %>% rename(symbol=hgnc_symbol) %>% left_join(homs) %>% mutate(type=ifelse(is.na(type), "WT", type))

ggplot(toPlot, aes(x=type, y=rpkm)) + geom_boxplot() + 
    scale_y_continuous(limits = c(0,10)) + 
    theme_boss() + 
    theme(
        axis.title.x = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")
    )

print(wilcox.test(toPlot$rpkm[toPlot$type=="HomDel"], toPlot$rpkm[toPlot$type=="WT"]))
print(summary(toPlot$rpkm[toPlot$type=="HomDel"]))
print(summary(toPlot$rpkm[toPlot$type=="WT"]))

print(length(toPlot$rpkm[toPlot$type=="HomDel"]))
print(length(toPlot$rpkm[toPlot$type=="WT"]))

print(length(unique(toPlot$entrez[toPlot$type=="HomDel"])))
print(length(unique(toPlot$sample[toPlot$type=="HomDel"])))
print(length(unique(toPlot$entrez[toPlot$type=="WT"])))
print(length(unique(toPlot$sample[toPlot$type=="WT"])))

write.table(toPlot, file="~/Desktop/homdels_vs_WT_syscans.tsv", quote = F, row.names = F, sep="\t")

## Get the genes we selected
cans = c("KAT2A", "KAT2B", "ROCK1", "ABI2", "PTK2", "IQGAP1", "COL4A1", "EIF4G1", "NBEA") 

plots = list()
for (c in cans){
    print(c)
    amplifs = syscans %>% subset(CNVGain==1 & symbol==c) %>% dplyr::select(sample, symbol) %>% unique %>% mutate(type="A")
    if(sum(amplifs$sample%in%unique(rpkms$sample))==0){
        stop("None of the samples has expression")
    }
    toPlot = rpkms %>% subset(hgnc_symbol==c) %>% rename(symbol=hgnc_symbol) %>% left_join(amplifs) %>% mutate(type=ifelse(is.na(type), "WT", type))
    
    p = ggplot(toPlot, aes(x=type, y=rpkm)) + geom_boxplot() + 
        scale_y_continuous(limits = c(0,100)) + 
        theme_boss() + 
        theme(
            axis.title.x = element_blank(),
            axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour = "black")
        )
    plots[[c]] = p
    
    print(wilcox.test(toPlot$rpkm[toPlot$type=="A"], toPlot$rpkm[toPlot$type=="WT"]))
    print(summary(toPlot$rpkm[toPlot$type=="A"]))
    print(summary(toPlot$rpkm[toPlot$type=="WT"]))
    
    print(length(toPlot$rpkm[toPlot$type=="A"]))
    print(length(toPlot$rpkm[toPlot$type=="WT"]))
    
    print(length(unique(toPlot$entrez[toPlot$type=="A"])))
    print(length(unique(toPlot$sample[toPlot$type=="A"])))
    print(length(unique(toPlot$entrez[toPlot$type=="WT"])))
    print(length(unique(toPlot$sample[toPlot$type=="WT"])))
    
}







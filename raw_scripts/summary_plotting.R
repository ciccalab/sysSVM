## Script for summarising and plotting data/mutations/annotations etc

## Get first the directories (samples may be in different directories)

library(RColorBrewer)


es = c("LP6005690-DNA_E02_vs_LP6005689-DNA_E02", 
       "LP6008280-DNA_F02_vs_LP6008264-DNA_F02",
       "LP6008202-DNA_F01_vs_LP6008201-DNA_F01",
       "LP6005935-DNA_C01_vs_LP6005934-DNA_C01",
       "LP6008031-DNA_E03_vs_LP6008032-DNA_A04")


## For the mutations I just need to load them - because runOncodriveClust script concatenated the data from all samples
load("~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/Rdata/muts_129_66_71_OACs_annovar_dbnsfp_oncodriveClust.Rdata")

## Save the data for 19,014
#muts = muts %>% subset(!is.na(entrez_19014))
#save(muts, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/Rdata/muts_129_66_71_OACs_annovar_dbnsfp_oncodriveClust_19014.Rdata")

## Plot the number of all mutations
samples2muts = muts %>% group_by(sample) %>% summarise(all_muts=n())
samples2muts$sample = factor(as.character(samples2muts$sample), levels = samples2muts$sample[order(samples2muts$all_muts, decreasing = F)])

sm = data.frame(type= names(summary(samples2muts$all_muts)),
           value=unname(c(summary(samples2muts$all_muts))))

p = ggplot(samples2muts, aes(x=sample, y=all_muts)) +
    geom_bar(stat = "identity") + xlab("samples") + ylab("Mutations (#)") +
    theme_boss() +
    theme(
        axis.text.x=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.ticks.x = element_blank()
    ) + 
    scale_y_continuous(breaks = seq(0, max(samples2muts$all_muts), max(samples2muts$all_muts)/10)) + 
    annotation_custom(tableGrob(sm, cols = NULL, rows = NULL), xmin=100, xmax=150, ymin=150000, ymax=200000) +
    ggtitle("All mutations (this sample order forced to all other plots)")

ns = c("nonsynonymous","stopgain","frameshift deletion","splicing","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","stoploss","frameshift substitution")
dam = c("nonsynonymous","frameshift deletion","frameshift insertion","frameshift substitution","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","splicing","stopgain","stoploss")
trunc = c("frameshift deletion","frameshift insertion","frameshift substitution","stopgain","stoploss") ## Always damaging==TRUE
non_trunc = c("nonsynonymous","splicing")
ns_vep=c("missense_variant", "splice_region_variant", "splice_donor_variant", "stop_gained", "splice_acceptor_variant", "stop_lost")

d = muts %>% count(sample, Func.refGene) %>% data.frame()
d = d %>% left_join(samples2muts) %>% mutate(perc=n/all_muts)
d$sample = factor(as.character(d$sample), levels = samples2muts$sample[order(samples2muts$all_muts, decreasing = F)])
## First plot check the fraction of exonic overall

n <- 15
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=col_vector[1:n])
cols = c(col_vector[1:2], "red", col_vector[4:n])

p1 = ggplot(d,
    aes(x=sample, y=perc, fill=Func.refGene)) +
    geom_bar(stat = "identity") + xlab("samples") + ylab("Mutations (fraction)") +
    scale_fill_manual(values = cols) +
    theme_boss() +
    theme(
        axis.text.x=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.ticks.x = element_blank()
    ) + labs(fill='Effect') + 
    ggtitle("All mutations")

## Get distribution of the percentage of exonic
tb1 = rbind(muts %>% group_by(sample) %>% count(Func.refGene) %>% mutate(Func.refGene=ifelse(Func.refGene=="exonic", "exonic", "other")) %>%
    group_by(sample, Func.refGene) %>% summarise(n=sum(n)) %>%left_join(samples2muts) %>% mutate(perc=n/all_muts) %>% subset(Func.refGene=="exonic") %>% .$perc %>% summary())

## Now check the categories of the exonic
samples2exonic = muts %>% subset(Func.refGene=="exonic") %>% group_by(sample) %>% summarise(exonic_muts=n())
samples2exonic$sample = factor(as.character(samples2exonic$sample), levels = samples2muts$sample[order(samples2muts$all_muts, decreasing = F)])

sm = data.frame(type= names(summary(samples2exonic$exonic_muts)),
                value=unname(c(summary(samples2exonic$exonic_muts))))

p2 = ggplot(samples2exonic, aes(x=sample, y=exonic_muts)) +
    geom_bar(stat = "identity") + xlab("samples") + ylab("Mutations (#)") +
    theme_boss() +
    theme(
        axis.text.x=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.ticks.x = element_blank()
    ) + 
    scale_y_continuous(breaks = seq(0, max(samples2exonic$exonic_muts), max(samples2exonic$exonic_muts)/10)) + 
    annotation_custom(tableGrob(sm, cols = NULL, rows = NULL), xmin=100, xmax=150, ymin=700, ymax=900) +
    ggtitle("Exonic mutations")



d = muts %>% subset(Func.refGene=="exonic") %>%count(sample, ExonicFunc.refGene) %>% data.frame()
d = d %>% left_join(samples2exonic) %>% mutate(perc=n/exonic_muts)
d$sample = factor(as.character(d$sample), levels = samples2muts$sample[order(samples2muts$all_muts, decreasing = F)])

p3 = ggplot(d ,
            aes(x=sample, y=perc, fill=ExonicFunc.refGene)) +
    geom_bar(stat = "identity") + xlab("samples") + ylab("Mutations (fraction)") +
    scale_fill_manual(values = cols) +
    theme_boss() +
    theme(
        axis.text.x=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.ticks.x = element_blank()
    ) + labs(fill='Effect') +
    ggtitle("Exonic mutations")

## Now check the damaging
samples2damaging = muts %>% subset(damaging) %>% group_by(sample) %>% summarise(damaging_muts=n())
samples2damaging = samples2damaging %>% full_join(samples2muts%>%select(sample)) ## Not all samples have damaging mutations
samples2damaging$damaging_muts[is.na(samples2damaging$damaging_muts)] = 0
samples2damaging$sample = factor(as.character(samples2damaging$sample), levels = samples2muts$sample[order(samples2muts$all_muts, decreasing = F)])


sm = data.frame(type= names(summary(samples2damaging$damaging_muts)),
                value=unname(c(summary(samples2damaging$damaging_muts))))

p4 = ggplot(samples2damaging, aes(x=sample, y=damaging_muts)) +
    geom_bar(stat = "identity") + xlab("samples") + ylab("Mutations (#)") +
    theme_boss() +
    theme(
        axis.text.x=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.ticks.x = element_blank()
    ) + 
    scale_y_continuous(breaks = seq(0, max(samples2damaging$damaging_muts), max(samples2damaging$damaging_muts)/10)) + 
    annotation_custom(tableGrob(sm, cols = NULL, rows = NULL), xmin=100, xmax=150, ymin=250, ymax=300) +
    ggtitle("Damaging mutations")

## Now the gain of function mutations
samples2gof = muts %>% subset(oncodriveClust) %>% group_by(sample) %>% summarise(gof_muts=n())
samples2gof = samples2gof %>% full_join(samples2muts%>%select(sample)) ## Not all samples have damaging mutations
samples2gof$gof_muts[is.na(samples2gof$gof_muts)] = 0
samples2gof$sample = factor(as.character(samples2gof$sample), levels = samples2muts$sample[order(samples2muts$all_muts, decreasing = F)])


sm = data.frame(type= names(summary(samples2gof$gof_muts)),
                value=unname(c(summary(samples2gof$gof_muts))))

p5 = ggplot(samples2gof, aes(x=sample, y=gof_muts)) +
    geom_bar(stat = "identity") + xlab("samples") + ylab("Mutations (#)") +
    theme_boss() +
    theme(
        axis.text.x=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.ticks.x = element_blank()
    ) + 
    scale_y_continuous(breaks = seq(0, max(samples2gof$gof_muts), 5)) + 
    annotation_custom(tableGrob(sm, cols = NULL, rows = NULL), xmin=10, xmax=60, ymin=10, ymax=15) +
    ggtitle("Gain-of-function mutations")


grid.arrange(p, p1, p2, p3, p4, ncol=2)
grid.arrange(p, p1, p2, p3, p4, p5, ncol=2)

## CNVs

## Gather all CNVs
# mainDirs = c("~/data/OAC/71_OAC/ascat/",
#              "~/data/OAC/87_OAC/66_ICGC/ascat/",
#              "~/data/OAC/129_OAC/ascat/")

mainDirs = c("~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/ascat/")


message("Getting CNVs...")
all_cnvs = data.frame()
count = 0
for(b in mainDirs){
    samples = list.dirs(b, recursive = F)
    for(s in samples){
        cat(s, "\n")
        fn = paste0(s, "/parsing_and_annotation/cnvs.Rdata")
        load(fn)
        sname = unlist(strsplit(s, "/"))
        sname = sname[length(sname)]
        d = cnvs[["df_cnvs_19014"]]
        all_cnvs = rbind(all_cnvs, d %>% mutate(sample=sname))
        count = count +1
    }
}
cat(paste0("Samples: ", count))

## Save raw data
cnvs = all_cnvs
save(cnvs, file="~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/Rdata/cnvs_21_literature_OACs.Rdata")


samples2cnvs = all_cnvs %>% subset(!is.na(entrez_19014)) %>% select(sample, entrez_19014, CNV_type_corrected) %>% unique %>% count(sample, CNV_type_corrected)
all_samples = rbind(samples2cnvs%>%select(sample)%>%unique%>%mutate(CNV_type_corrected="Gain"), samples2cnvs%>%select(sample)%>%unique%>%mutate(CNV_type_corrected="Loss"))
samples2cnvs = samples2cnvs %>% full_join(all_samples)
samples2cnvs$n[is.na(samples2cnvs$n)] = 0
samples2cnvs = samples2cnvs %>% subset(!is.na(CNV_type_corrected))

sm1 = data.frame(type= names(summary(samples2cnvs$n[samples2cnvs$CNV_type_corrected=="Loss"])),
                value=unname(c(summary(samples2cnvs$n[samples2cnvs$CNV_type_corrected=="Loss"]))))

sm2 = data.frame(type= names(summary(samples2cnvs$n[samples2cnvs$CNV_type_corrected=="Gain"])),
                 value=unname(c(summary(samples2cnvs$n[samples2cnvs$CNV_type_corrected=="Gain"]))))

p = ggplot(samples2cnvs %>% subset(n>0), aes(x=sample, y=n, fill=CNV_type_corrected)) +
    geom_bar(stat = "identity", position = "dodge") +
    ylab("Genes (#)") +
    xlab("Samples") + theme_boss() +
    theme(
        axis.text.x=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.ticks.x = element_blank()
    ) +
    ggtitle("Gains>=2*ploidy; Losses=CN<2")

sm1 = tableGrob(sm1, cols = NULL, rows = NULL)
sm2 = tableGrob(sm2, cols = NULL, rows = NULL)

grid.arrange(arrangeGrob(sm1, sm2, ncol=2), 
             arrangeGrob(p, nrow=1, ncol=1), heights=c(0.2, 0.8))


## SVs
mainDirs = c("~/athena/data/OAC/71_OAC/manta/",
             "~/athena/data/OAC/87_OAC/66_ICGC/manta/",
             "~/athena/data/OAC/129_OAC/manta/")


message("Getting SVs...")
all_svs = data.frame()
count = 0
ss = NULL
for(b in mainDirs){
    samples = list.dirs(b, recursive = F)
    for(s in samples){
        cat(s, "\n")
        fn = paste0(s, "/parsing_and_annotation/svs.Rdata")
        if(file.exists(fn)){
            load(fn)
            sname = unlist(strsplit(s, "/"))
            sname = sname[length(sname)]
            ss = c(ss, sname)
            all_svs = rbind.fill(all_svs, svs %>% mutate(sample=sname))
            count = count +1 
        }else{
            next
        }

    }
}
cat(paste0("Samples: ", count))

## Save raw data
svs = all_svs
svs[,2:6][is.na(svs[,2:6])] = 0
save(svs, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/Rdata/svs.Rdata")
## Save the data on the 19,014
#svs = svs %>% subset(!is.na(entrez_19014)) %>% data.frame()
#save(svs, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/Rdata/svs_19014.Rdata")

load("~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/Rdata/svs.Rdata")
samples2svs2type = svs %>% select(sample, gene, DEL, DUP, INV, BND, INS) %>% gather(type, value, -sample, -gene) %>% subset(value!=0) %>% group_by(sample, type) %>% summarise(n=sum(value))
samples2svs = samples2svs2type %>% group_by(sample) %>% summarise(svs=sum(n))
samples2svs$sample = factor(as.character(samples2svs$sample), levels = unique(samples2svs$sample[order(samples2svs$svs, decreasing = F)]))

samples2svs2type = samples2svs2type %>% left_join(samples2svs) %>% mutate(perc=(n/svs)*100)
samples2svs2type$sample = factor(as.character(samples2svs2type$sample), levels = unique(samples2svs$sample[order(samples2svs$svs, decreasing = F)]))


sm = data.frame(type= names(summary(samples2svs$svs)),
                value=unname(c(summary(samples2svs$svs))))

p1 = ggplot(samples2svs, aes(x=sample, y=svs)) +
    geom_bar(stat = "identity") + xlab("samples") + ylab("Genes (#)") +
    theme_boss() +
    theme(
        axis.text.x=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.ticks.x = element_blank()
    ) + 
    scale_y_continuous(breaks = seq(0, max(samples2svs$svs), 100)) + 
    annotation_custom(tableGrob(sm, cols = NULL, rows = NULL), xmin=50, xmax=100, ymin=1000, ymax=1200) +
    ggtitle("All SVs (this sample order forced to all other plots)")



p2 = ggplot(samples2svs2type, aes(x=sample, y=perc, fill=type)) +
    geom_bar(stat = "identity") + xlab("samples") + ylab("Genes (%)") +
    theme_boss() +
    theme(
        axis.text.x=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.ticks.x = element_blank()
    ) + 
    ggtitle("SV types")

grid.arrange(p1, p2, nrow=1)


## ---------------------------------------
##      Plots for the drivers
## ---------------------------------------

load("~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/training_set_noScale.Rdata")
load("~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/validation_set_noScale.Rdata")
training_ns = training_ns %>% tibble::rownames_to_column() %>% separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")
validation_ns = validation_ns %>% tibble::rownames_to_column() %>% separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")
cohort = rbind.fill(training_ns, validation_ns)
## Get count for basic alterations
train_toPlot = training_ns %>% select(sample, no_TRUNC_muts, no_NTDam_muts, no_GOF_muts, BND, INS, INV, CNVGain) %>% 
    mutate(CNVGain=as.numeric(as.character(CNVGain))) %>% gather(type, value, -sample) %>% 
    subset(value==1) %>% group_by(type) %>% summarise(n=n()) %>% mutate(all=4091, perc=(as.numeric(n)/all)*100) %>% 
    rbind(c("Homo_dels", training_ns %>% subset(Copy_number==0) %>% nrow, 4091, (training_ns %>% subset(Copy_number==0) %>% nrow)/4091)) %>% 
    rbind(c("Multiple_hits", training_ns %>% subset(Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)) %>% nrow, 4091, (training_ns %>% subset(Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)) %>% nrow)/4091)) %>%
    mutate(label="training set (4091 obs)")


validation_toPlot = validation_ns %>% select(sample, no_TRUNC_muts, no_NTDam_muts, no_GOF_muts, BND, INS, INV, CNVGain) %>% 
    mutate(CNVGain=as.numeric(as.character(CNVGain))) %>% gather(type, value, -sample) %>% 
    subset(value==1) %>% group_by(type) %>% summarise(n=n()) %>% mutate(all=112898, perc=(as.numeric(n)/all)*100) %>% 
    rbind(c("Homo_dels", training_ns %>% subset(Copy_number==0) %>% nrow, 112898, (training_ns %>% subset(Copy_number==0) %>% nrow)/112898)) %>% 
    rbind(c("Multiple_hits", training_ns %>% subset(Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)) %>% nrow, 112898, (training_ns %>% subset(Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)) %>% nrow)/112898)) %>% 
    mutate(label="prediction set (112,898 obs)")

toPlot = rbind(train_toPlot, validation_toPlot)
ggplot(toPlot, aes(x=type, y=as.numeric(perc))) + 
    geom_bar(stat = "identity", position="dodge", color="black", fill="grey50") +
    facet_wrap(~label) +
    xlab("") + ylab("Drivers (%)") +
    theme(
        axis.text.x=element_text(angle=90),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")

    ) + 
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10))

write.table(toPlot, file="~/Desktop/drivers_261.tsv", quote = F, sep = "\t", row.names = F)


## And the same for the sys-candidates
load("~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.noAmp.top10.plusCGC.Rdata")
toPlot_noAmp = gsea.noAmp.top10.plusCGC[["genes"]] %>% subset(gene_type!="cgc") %>% select(sample, no_TRUNC_muts, no_NTDam_muts, no_GOF_muts, BND, INS, INV, CNVGain) %>% 
    mutate(CNVGain=as.numeric(as.character(CNVGain))) %>% gather(type, value, -sample) %>% 
    subset(value==1) %>% group_by(type) %>% summarise(n=n()) %>% mutate(all=2598, perc=(as.numeric(n)/all)*100) %>% 
    rbind(c("Homo_dels", training_ns %>% subset(Copy_number==0) %>% nrow, 2598, (training_ns %>% subset(Copy_number==0) %>% nrow)/2598)) %>% 
    rbind(c("Multiple_hits", training_ns %>% subset(Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)) %>% nrow, 2598, (training_ns %>% subset(Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)) %>% nrow)/2598)) %>% 
    mutate(label="Sys-candidates without Amplification (2,598 obs)")


load("~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/gsea.withAmp.top10.plusCGC.Rdata")
toPlot_withAmp = gsea.withAmp.top10.plusCGC[["genes"]] %>% subset(gene_type!="cgc") %>% select(sample, no_TRUNC_muts, no_NTDam_muts, no_GOF_muts, BND, INS, INV, CNVGain) %>% 
    mutate(CNVGain=as.numeric(as.character(CNVGain))) %>% gather(type, value, -sample) %>% 
    subset(value==1) %>% group_by(type) %>% summarise(n=n()) %>% mutate(all=2608, perc=(as.numeric(n)/all)*100) %>% 
    rbind(c("Homo_dels", training_ns %>% subset(Copy_number==0) %>% nrow, 2608, (training_ns %>% subset(Copy_number==0) %>% nrow)/2608)) %>% 
    rbind(c("Multiple_hits", training_ns %>% subset(Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)) %>% nrow, 2608, (training_ns %>% subset(Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)) %>% nrow)/2608)) %>% 
    mutate(label="Sys-candidates with Amplifications (2,608 obs)")

toPlot = rbind(toPlot_noAmp, toPlot_withAmp)

ggplot(toPlot, aes(x=type, y=as.numeric(perc))) + 
    geom_bar(stat = "identity", position="dodge", color="black", fill="grey50") +
    facet_wrap(~label) +
    xlab("") + ylab("Drivers (%)") +
    theme(
        axis.text.x=element_text(angle=90),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")
        
    ) + 
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10))

write.table(toPlot, file="~/Desktop/drivers_syscans.tsv", quote = F, sep = "\t", row.names = F)







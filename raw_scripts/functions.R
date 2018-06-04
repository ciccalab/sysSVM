## *****************************************************************************
##           R functions for parsing of ICGC OAC samples 
## *****************************************************************************


## NOTE: this makes the code less portable
options(stringsAsFactors = F)
## ---------------------------------------
source("http://bioconductor.org/biocLite.R")

pkgs  =  c("VariantAnnotation", "GenomicRanges")
for (pkg in pkgs){
    if (!pkg %in% rownames(installed.packages())) { 
        biocLite(pkg)
        library(pkg, character.only = TRUE, quietly = TRUE)
    } else { 
        library(pkg, character.only = TRUE, quietly = TRUE)
    }
}

## make sure you load dplyr last not to mask select etc
pkgsR  =  c('plyr',
            'tidyr',
            'xlsx',
            'data.table',
            'doMC',
            'foreach',
            'ggplot2',
            'grid',
            'gridExtra',
            'scales',
            'snow',
            'dplyr')
for (pkgR in pkgsR){
    if (!pkgR %in% rownames(installed.packages())) { 
        #install.packages(pkgR, dependencies=TRUE, INSTALL_opts = c('--no-lock'))
        library(pkgR, character.only = TRUE, quietly = TRUE)
    } else { 
        library(pkgR, character.only = TRUE, quietly = TRUE)
    }
}

## Load dplyr always last
detach("package:dplyr", unload=TRUE)
library(dplyr)

## Other configurations
cl = makeCluster(10)

ns = c("nonsynonymous","stopgain","frameshift deletion","splicing","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","stoploss","frameshift substitution")
dam = c("nonsynonymous","frameshift deletion","frameshift insertion","frameshift substitution","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","splicing","stopgain","stoploss")
trunc = c("frameshift deletion","frameshift insertion","frameshift substitution","stopgain","stoploss") ## Always damaging==TRUE
non_trunc = c("nonsynonymous","splicing")

magicNumber <- function (x) {
    samples = length(unique(x$sample))
    c.samples = length( unique( subset(x, cgc_esophagous)$sample ) )
    r.sample = length(unique(subset(x, cancer_type=='rst')$sample))
    
    u.genes = with( unique(x[,c("symbol_19014","cancer_type")]), c('tot'=length(symbol_19014), as.array(table(cancer_type)[c('cancer','rst')])))
    u.c.genes = with( unique(subset(x, cgc_esophagous)[,c("symbol_19014","cancer_type")]), c('tot'=length(symbol_19014), as.array(table(cancer_type)[c('cancer','rst')])))
    genes = paste(u.genes, " (", u.c.genes, ")", sep=""); names(genes) = names(u.genes)
    genes['rst'] = u.genes['rst']
    
    r.genes = with( (x[,c("symbol_19014","cancer_type")]),  c('tot.entry'=length(symbol_19014), as.array(table(cancer_type)[c('cancer','rst')])))
    r.c.genes = with( (subset(x, cgc_esophagous)[,c("symbol_19014","cancer_type")]),  c('tot.entry'=length(symbol_19014), as.array(table(cancer_type)[c('cancer','rst')])))
    entry = paste(r.genes, " (", r.c.genes, ")", sep=""); names(entry) = names(r.genes)
    entry['rst'] = r.genes['rst']
    
    y = c('samples' = paste0(samples, " (", c.samples, ")"), genes, 'samples.rst'=r.sample,entry)
    names(y)[7:8] =c('cancer.entry','rst.entry')
    y
}

base_breaks_x = function(x, xend, br, la){
    d <- data.frame(y=-Inf, yend=-Inf, x=x, xend=xend)
    list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE), scale_x_reverse(breaks=br, labels=la))
}

base_breaks_y = function(y, yend, br, la ){
    d <- data.frame(x=-Inf, xend=-Inf, y=y, yend=yend)
    list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE), scale_y_discrete(breaks=br, labels=la))
}

get_boxplot = function(t, main="", log10=F){
    p = wilcox.test(t$value[t[,1]=="cancer"], t$value[t[,1]=="rst"] )$p.value
    r = NULL
    if(log10!=T){
        r = ggplot(t, aes_string(fill=colnames(t)[1], y=colnames(t)[3], x=colnames(t)[2] )) +geom_boxplot(col="black")+theme_boss()+
            scale_fill_manual(values=c("white","grey"))+ggtitle(paste0(main," P=",round(p,3)))
    }else{
        r = ggplot(t, aes_string(fill=colnames(t)[1], y=colnames(t)[3], x=colnames(t)[2] )) +geom_boxplot(col="black")+theme_boss()+
            scale_fill_manual(values=c("white","grey"))+ggtitle(paste0(main," P=",round(p,3)))+scale_y_log10()
    }
    r
}

roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

get_barplot = function(t, main=""){
    t = as.matrix(t)
    s = apply(t,1,sum)
    p = c()
    for(i in 1:ncol(t)) p[i] = fisher.test( matrix(c(t[,i],s),nc=2))$p.value
    names(p) = colnames(t)  
    t = t/apply(t,1,sum)
    m = melt(t)
    m$p = p[match(m[,2], names(p))]
    p = unique(ddply(m, colnames(m)[2], mutate, value = value[which.max(value)], p=round(p[which.max(value)],3), pr = p[which.max(value)]))
    fdr = p.adjust(unique(p[,c(2,5)])$pr,'BH')
    names(fdr) = unique(p[,2])
    p$fdr = fdr[p[,2]]
    p$star = ""
    p$star[p$fdr<=0.1] = "*"
    p$print = paste0(p$p,p$star)
    ymax= roundUpNice(max(m$value))
    mid = ymax/2
    r = 
        ggplot(m, aes_string(fill=colnames(m)[1], x=colnames(m)[2], y=colnames(m)[3] )) +
        geom_bar(stat="identity", position=position_dodge(), col="black")+
        geom_text( data=p, aes_string(x = colnames(p)[2], y = colnames(p)[3], label=colnames(p)[8]), nudge_y = 0.025)+
        theme_boss()+scale_fill_manual(values=c("white","grey"))+ggtitle(main)#+theme(rect             = element_blank())+
    # base_breaks_y(0, ymax, c(0,mid, ymax), as.character( c(0,mid, ymax) ))
    # base_breaks_y(0, 1, c(0, 1), as.character( c(0,1) ))
    r
}

get.dataset.barplot <- function (d, varid) {
    m = melt(d, id.vars=c('Cancer_type'), measure.vars = varid)
    
    tmp = subset(m, Cancer_type!="EAC")
    
    tmp$Cancer_type='ALL'
    m = rbind(m, tmp)
    t = table(m$Cancer_type)
    med = ddply(m, .(Cancer_type), summarise, m = round(median(value),2))
    
    m$code = paste(m[,1],' (', t[m$Cancer_type], ',', med[match(m$Cancer_type,med[,1]),2], ')',sep="", coll="")
    
    code = unique(m$code);
    names(code)=unique(m$Cancer_type)
    
    tmp = subset(med, Cancer_type!="ALL")
    ix = order(tmp[,2], decreasing = F)
    
    m$code = factor(as.character(m$code), levels = c(code[tmp[ix,1]],code["ALL"]), ordered = T)
    m$Cancer_type = factor(as.character(m$Cancer_type), levels = names(c(code[tmp[ix,1]],code["ALL"])), ordered = T)
    m
}

get.dataset.barplot.double.sorting <- function (d, varid) {
    m = melt(d, id.vars=c('Cancer_type'), measure.vars = varid)
    
    tmp = subset(m, Cancer_type!="EAC")
    
    tmp$Cancer_type='ALL'
    m = rbind(m, tmp)
    t = table(m$Cancer_type)
    med = ddply(m, .(Cancer_type), summarise, m = round(median(value),2), mn =round(mean(value),2))
    
    m$code = paste(m[,1],' (', t[m$Cancer_type], ',', med[match(m$Cancer_type,med[,1]),2],',', med[match(m$Cancer_type,med[,1]),3], ')',sep="", coll="")
    
    code = unique(m$code);
    names(code)=unique(m$Cancer_type)
    
    tmp = subset(med, Cancer_type!="ALL")
    ix = order(tmp[,2],tmp[,3], decreasing = F)
    
    m$code = factor(as.character(m$code), levels = c(code[tmp[ix,1]],code["ALL"]), ordered = T)
    m$Cancer_type = factor(as.character(m$Cancer_type), levels = names(c(code[tmp[ix,1]],code["ALL"])), ordered = T)
    m
}

theme_boss <- function(base_size = 12, base_family = "sans"){
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(axis.text        = element_text(size=rel(0.85)),
              axis.title.x     = element_text(size=rel(1)),
              axis.title.y     = element_text(size=rel(1)),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),  
              panel.background = element_blank(), 
              panel.border=element_blank()
        )
}

lapply_pb = function(X, FUN, ...){
    env <- environment()
    pb_Total <- length(X)
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   
    
    # wrapper around FUN
    wrapper <- function(...){
        curVal <- get("counter", envir = env)
        assign("counter", curVal +1 ,envir=env)
        setTxtProgressBar(get("pb", envir=env), curVal +1)
        FUN(...)
    }
    res <- lapply(X, wrapper, ...)
    close(pb)
    res
}

mutation_frequency = function(x, target=51.18932)  round(nrow(x)/target,1)

get_gname=function(data){
    gname=sapply(strsplit(sapply(strsplit(sapply(strsplit(data$Gene.refGene,";"),function(x) x[1]),split="\\("),function(x) x[1]),split="\\,"),function(x) x[1])
    gname
}

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

get_19014 = function(x, geneSymbols){
    require(dplyr)
    
    ix = which(colnames(x)=="Gene.refGene")
    x$symbol_19014 = sapply(x[,ix], function(x) unlist(strsplit(x, ","))[which(unlist(strsplit(x, ",")) %in% geneSymbols$symbol)][1] )
    x = x %>% left_join(geneSymbols%>%select(symbol, Entrez)%>%rename(symbol_19014=symbol, entrez_19014=Entrez))
    x
}

set_mutation_code = function(x){
    ns         = c("NS",'SG','SL','SP','FD','FI','NFD','NFI')
    names(ns)  = c("nonsynonymous" , "stopgain", "stoploss",  "splicing" , "frameshift deletion", "frameshift insertion" , "nonframeshift deletion", "nonframeshift insertion")
    x$mtype = ns[match(x$ExonicFunc.refGene,names(ns))];
    x$mtype[which(x$mtype=="NS" & x$damaging)]='NSD'
    x$assignedCNV[which(is.na(x$assignedCNV))]=''
    x
}

set_IGV_code = function(x){
    x$IGV = paste(x$chr,":",x$end,sep="",coll="")
    x
}

correct_dam = function(x){
    y=x
    y[which(y$ExonicFunc.refGene%in%c("stopgain","stoploss")),'damaging']=T
    return(y)
}

MutationTaster_convert <- function(df) {
    df <- as.data.frame(df,stringsAsFactors = F)
    colnames(df) <- c("score", "pred")
    df$score <- as.numeric(df$score)
    df$converted_score <- NA
    for (i in 1:nrow(df)){
        if (df$pred[i]=="D" | df$pred[i]=="A"){
            df$converted_score[i] <- df$score[i]
        }else if (df$pred[i]=="N" | df$pred[i]=="P"){
            df$converted_score[i] <- 1-df$score[i]
        }else if (df$pred[i]=="."){
            df$converted_score[i] <- "."
        }
    }
    
    return(paste(df$converted_score,collapse = ";"))
}

MutationAssessor_convert <- function(df) {
    df <- as.data.frame(df,stringsAsFactors = F)
    colnames(df) <- c("score")
    df$score <- as.numeric(df$score)
    df$converted_score <- NA
    ma_range <- 5.545+5.975
    ma_min <- -5.545
    for (i in 1:nrow(df)){
        if (!is.na(df$score[i])){
            df$converted_score[i] <- (df$score[i]-ma_min)/ma_range
        }else if (is.na(df$score[i])){
            df$converted_score[i] <- "."
        }
    }
    return(paste(df$converted_score,collapse = ";"))
}

FATHMM_convert <- function(df) {
    df <- as.data.frame(df,stringsAsFactors = F)
    colnames(df) <- c("score")
    df$score <- as.numeric(df$score)
    df$converted_score <- NA
    ftm_range <- 16.13+10.64
    ftm_min <- -16.13
    for (i in 1:nrow(df)){
        if (!is.na(df$score[i])){
            df$converted_score[i] <- 1-(df$score[i]-ftm_min)/ftm_range
        }else if (is.na(df$score[i])){
            df$converted_score[i] <- "."
        }
    }
    
    return(paste(df$converted_score,collapse = ";"))
}

# === RNA ======

fpkmToTpm <- function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToFpkm <- function(counts, tlen){
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(tlen) - log(N) )
}

countToTpm <- function(counts, effLen)
{
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}

getClinicalInfo = function(clinic_dir = "Clinical_data/old_clinical/icgc-dataset-1473066214193"){
    
    ## Get the clinical data to check sample names
    message("Reading clinical data...")
    donor = read.delim(paste0(clinic_dir, "/donor.tsv"), header = T)
    donor_exposure = read.delim(paste0(clinic_dir, "/donor_exposure.tsv"), header = T)
    donor_family = read.delim(paste0(clinic_dir, "/donor_family.tsv"), header = T)
    donor_therapy = read.delim(paste0(clinic_dir, "/donor_therapy.tsv"), header = T)
    sample = read.delim(paste0(clinic_dir, "/sample.tsv"), header = T)
    specimen = read.delim(paste0(clinic_dir, "/specimen.tsv"), header = T)
    ## Combine them
    clinic = sample %>% left_join(specimen) %>% left_join(donor) %>% left_join(donor_exposure) %>% left_join(donor_family) %>% left_join(donor_therapy)
    return(clinic)
}


getCohortInfo = function(mut_dir = "129_OAC/129_raw_data/strelka/",
                         icgc_clinic = T,
                         clinic_dir = "Clinical_data/release_24/",
                         cnv_dir = "129_OAC/129_raw_data/ascat/",
                         sample_annot_fn = NULL){
    
    ## Get the annotation of the 18 OAC samples from the pilot data
    # load("18_OAC/Rdata/dataset_mutations.Rdata")
    # samples_18_OAC = names(muts)
    # rm(muts)
    
    ## Get the clinical data to check sample names
    if(icgc_clinic){
        message("Reading clinical data...")
        donor = read.delim(paste0(clinic_dir, "/donor.tsv"), header = T)
        donor_exposure = read.delim(paste0(clinic_dir, "/donor_exposure.tsv"), header = T)
        donor_family = read.delim(paste0(clinic_dir, "/donor_family.tsv"), header = T)
        donor_therapy = read.delim(paste0(clinic_dir, "/donor_therapy.tsv"), header = T)
        sample = read.delim(paste0(clinic_dir, "/sample.tsv"), header = T)
        specimen = read.delim(paste0(clinic_dir, "/specimen.tsv"), header = T)
        ## Combine them
        clinic = sample %>% left_join(specimen) %>% left_join(donor) %>% left_join(donor_exposure) %>% left_join(donor_family) %>% left_join(donor_therapy)
    }else if(!is.null(sample_annot_fn)){
        clinic = read.table(sample_annot_fn, header = T)
    }

    
    # === Mutation data (SNV + indel as this stage, I am just checking sample names) ======
    message("Reading mutation data...")
    mut_dirs = list.files(path=mut_dir)
    mut_dirs = as.list(mut_dirs)
    
    ## Check sample names from SNV data
    muts = data.frame(directory=unlist(mut_dirs)) %>% mutate(submitted_sample_id=directory) %>% mutate(submitted_sample_id=strsplit(submitted_sample_id, "_vs_")) %>% unnest() 
    if(icgc_clinic){
        muts = muts %>% left_join(clinic%>%select(submitted_sample_id, specimen_type))
    }else if(!is.null(sample_annot_fn)){
        muts = muts %>% left_join(clinic%>%select(submitted_sample_id, specimen_type))
    }
    muts = muts %>% mutate(strelka_data=TRUE)
    
    
    # === Copy number data (SNV + indel as this stage, I am just checking sample) ======
    message("Reading copy number variation data...")
    cnv_dirs = list.files(path=cnv_dir)
    cnv_dirs = as.list(cnv_dirs)
    
    cnvs = data.frame(directory=unlist(cnv_dirs)) %>% mutate(submitted_sample_id=directory) %>% mutate(submitted_sample_id=strsplit(submitted_sample_id, "_vs_")) %>% unnest()
    if(icgc_clinic){
        cnvs = muts %>% left_join(clinic%>%select(submitted_sample_id, specimen_type))
    }else if(!is.null(sample_annot_fn)){
        cnvs = muts %>% left_join(clinic%>%select(submitted_sample_id, specimen_type))
    }
    cnvs = cnvs %>% mutate(ascat_data=TRUE)
    
    ## Also, get the cellularity from ascat for each sample
    # message("Getting ASCAT cellularity...")
    # cellular = NULL
    # for (i in cnv_dirs){
    #     stats = read.csv( paste0(cnv_dir,i,"/1.4/results/",i,".samplestatistics.csv"), sep=" ", header = F)
    #     colnames(stats) = c("type", "value")
    #     stats = stats$value[stats$type=="rho"]
    #     cellular = rbind(cellular, data.frame(directory=i, ascat_rho=stats))
    # }
    # cnvs = cnvs %>% left_join(cellular)
    
    ## Combine two tables
    res = muts %>% left_join(cnvs)
    
    ## Annotate 18 samples
    #found_18 = res %>% mutate(in_18=ifelse(Primary_tumour_solid_tissue%in%samples_18_OAC, TRUE, FALSE)) %>% subset(in_18) %>% .$Primary_tumour_solid_tissue
    #setdiff(samples_18_OAC, found_18)
    ## only 16/18 found in the 129 - I think this is because they removed samples that were stage 3 and both of them were stage 3 at the time of diagnosis
    ## Produce the table
    #res = res %>% mutate(in_18=ifelse(Primary_tumour_solid_tissue%in%samples_18_OAC, TRUE, FALSE))
    ## Join clinical data
    #res = res %>% left_join(clinic, by=c("Primary_tumour_solid_tissue"="submitted_sample_id"))
    
    ## Add mutational signature grouping
    # if(!is.null(group_fn)){
    #     groups = read.table(group_fn, header = T, sep = "\t")
    #     res = res %>% left_join(groups, by=c("directory"="Sample"))
    # }
    
    return(res) 
}

getSNVs = function(mut_fn, g){
    
    ## Get sample names etc
    if(is.null(mut_fn)){
        stop("SNV file does not exist")
    }
    
    x = readVcf(mut_fn,genome=g )
    
    ## Get information from the VCF object

    snv_df = data.frame(
        chr= as.character(seqnames(rowRanges(x))),
        start=start(ranges(rowRanges(x))),
        end=end(ranges(rowRanges(x))),
        ref=as.character(rowRanges(x)$REF),
        alt=sapply(strsplit(names(rowRanges(x)[,1]), "\\/"), function(x) x[2]),
        ReadCount = info(x)$ReadCount,
        ReadCountControl = info(x)$ReadCountControl,
        VariantAlleleCount = info(x)$VariantAlleleCount,
        VariantAlleleCountControl = info(x)$VariantAlleleCountControl,
        freq = info(x)$VariantAlleleFrequency,
        VariantStrandBias = info(x)$VariantStrandBias
    )
    return(snv_df)
}


getIndels = function(indel_fn, g){
    
    ## Get sample names etc
    if(is.null(indel_fn)){
        stop("Indel file does not exist")
    }
    
    x = readVcf( indel_fn, genome="hg19" )
    
    ## Get information from the VCF object

    indel_df = data.frame(
        chr= as.character(seqnames(rowRanges(x))),
        start=start(ranges(rowRanges(x))),
        end=end(ranges(rowRanges(x))),
        ref=as.character(rowRanges(x)$REF),
        alt=sapply(strsplit(names(rowRanges(x)[,1]), "\\/"), function(x) x[2]),
        ReadCount = geno(x)$DP[,"TUMOR"],
        ReadCountControl = geno(x)$DP[,"NORMAL"],
        VariantAlleleCount = geno(x)$TAR[,2,1],
        VariantAlleleCountControl = NA,
        freq = (geno(x)$TAR[,2,1])/(geno(x)$DP[,"TUMOR"]),
        VariantStrandBias = NA
    )

    return(indel_df)
}

## Careful here with the version of human genome (cmds below use hg19 add a parameter fi you want that changed or change it directly to the commands)
annotateMutations = function(snv=NULL, indel=NULL, g=NULL,
                             annovar_path="/mnt/lustre/users/k1469280/FC/Software/annovar_2015_12_14",
                             dbnsfp_path="/mnt/lustre/users/k1469280/FC/DB/dbNSFPv3.0b2a/",
                             dbNSFP2_ver = "search_dbNSFP30b2a",
                             gene_symbols_fn="/mnt/lustre/users/k1469280/mourikisa/data/19014_GeneSymbols.Rdata",
                             save_dir=NULL){
    
    if(is.null(snv) | is.null(indel) | is.null(g) | is.null(save_dir)){
        stop("Please provide all parameters")
    }
    
    if(!file.exists(save_dir)){
        dir.create(save_dir)
    }
    
    message("Creating ANNOVAR's input")
    ## SNVs
    ann = unique(snv[,1:5])
    ann = unrowname(ann)
    write.table(ann, file=paste0(save_dir, "/annovar.snv"), col.names=F, row.names=F, quote=F, sep="\t")
    ## Indels
    ann = unique(indel[,1:5])
    ann = unrowname(ann)
    write.table(ann, file=paste0(save_dir, "/annovar.indel"), col.names=F, row.names=F, quote=F, sep="\t")

    snv_fn = paste0(save_dir, "/annovar.snv")
    snv_out = paste0(save_dir, "/annovar.snv.out")

    indel_fn = paste0(save_dir, "/annovar.indel")
    indel_out = paste0(save_dir, "/annovar.indel.out")

    ## Run ANNOVAR for SNVs
    oldwd = getwd()
    setwd(annovar_path)
    cmd = paste0("perl table_annovar.pl ", snv_fn, " humandb/ -buildver ", g , " -out ", snv_out, " -remove -protocol refGene,esp6500siv2_all,esp6500siv2_ea,1000g2015aug_all,1000g2015aug_eas,snp138  -operation g,f,f,f,f,f")
    system(noquote(cmd))
    cmd = paste0("perl table_annovar.pl ", indel_fn, " humandb/ -buildver ", g, " -out ", indel_out, " -remove -protocol refGene,esp6500siv2_all,esp6500siv2_ea,1000g2015aug_all,1000g2015aug_eas,snp138  -operation g,f,f,f,f,f")
    system(noquote(cmd))
    setwd(oldwd)

    ## Parse ANNOVAR's output
    message("Parsing ANNOVAR's results...")
    as = read.delim(file=paste0(snv_out, ".hg19_multianno.txt"))
    ai = read.delim(file=paste0(indel_out, ".hg19_multianno.txt"))

    ia = with(as,paste0( Chr,'.',Start,'.',End,'.',Ref,'.',Alt))
    is = with(snv,paste0( chr,'.',start,'.',end,'.',ref,'.',alt))
    asnv= cbind(snv,as[match(is,ia),c("Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","esp6500siv2_all","esp6500siv2_ea","X1000g2015aug_all","X1000g2015aug_eas","snp138")])


    ia = with(ai,paste0( Chr,'.',Start,'.',End,'.',Ref,'.',Alt))
    is = with(indel,paste0( chr,'.',start,'.',end,'.',ref,'.',alt))
    aindel = cbind(indel,ai[match(is,ia),c("Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","esp6500siv2_all","esp6500siv2_ea","X1000g2015aug_all","X1000g2015aug_eas","snp138")])


    muts = rbind(asnv, aindel)

    muts=fix_splicing(muts)
    muts=fix_exonic(muts)
    muts=fix_exonicFunc(muts)
    muts=is_nonsilent(muts)
    muts=set_IGV_code(muts)

    save(muts, file=paste0(save_dir, "/muts_annovar.Rdata"))

    # # === DBNSFP ====
    message("Running DBNSP...")
    ns = subset(muts, ExonicFunc.refGene%in%c("nonsynonymous SNV","nonsynonymous"))
    df_ns <- unique(ns[,c("chr","end","ref","alt")])
    ns_infile = paste0(save_dir, "/ns.oac")
    write.table(df_ns, file=ns_infile, col.names=F, row.names=F, quote=F, sep="\t")

    sp = subset(muts, Func.refGene%in%c("splicing"))
    df_sp <- unique(sp[,c("chr","end","ref","alt")])
    sp_infile = paste0(save_dir, "/sp.oac")
    write.table(df_sp, file=sp_infile, col.names=F, row.names=F, quote=F, sep="\t")

    oldwd = getwd()
    ## Non-synonymous
    ns_outfile = paste0(save_dir, "/ns.oac.dbnsfp")
    setwd(dbnsfp_path)
    if(nrow(ns)>0){
        cmd = paste("java -Xmx4g ",dbNSFP2_ver," -i ",ns_infile," -o ",ns_outfile, " -v ", g," ",coll="",sep="")
        system(noquote(cmd))
    }

    ## splicing
    outfile <- paste0(save_dir, "/sp.oac.dbnsfp")
    scoutfile <- paste0(save_dir, "/sp.oac.dbnsfp.sc")
    if(nrow(sp)>0){
        cmd <- paste("java -Xmx4g ",dbNSFP2_ver," -i ",sp_infile," -o ",outfile, " -s ", scoutfile , " -v ", g, " ",coll="",sep="")
        system(noquote(cmd))
    }
    setwd(oldwd)

    ## **********************************
    ## Get dbnsfp results for non-silent
    ## **********************************
    message("Parsing DBNSFP results...")
    if(nrow(ns)>0){
        ns = read.delim(ns_outfile,stringsAsFactor=F,quote = "")
        if(nrow(ns)>0){
            # Add MutationTaster_converted_score column
            mt_convert <- lapply(1:nrow(ns), function(x) cbind( unlist(strsplit(as.character(ns$MutationTaster_score[x]), ";")), unlist(strsplit(ns$MutationTaster_pred[x], ";"))) )
            ns$MutationTaster_converted_score <- unlist(lapply(mt_convert, MutationTaster_convert))
            ## Add MutationAssessor_converted_score column
            ma_convert <- lapply(1:nrow(ns), function(x) unlist(strsplit(as.character(ns$MutationAssessor_score[x]), ";")))
            ns$MutationAssessor_converted_score <- unlist(lapply(ma_convert, MutationAssessor_convert))
            ## Add FATHMM_converted_score column
            ftm_convert <- lapply(1:nrow(ns),function(x) unlist(strsplit(as.character(ns$FATHMM_score[x]), ";")))
            ns$FATHMM_converted_score <- unlist(lapply(ftm_convert, FATHMM_convert))
            ns <- ns %>% mutate(key=paste("chr",hg19_chr, ".", hg19_pos.1.based., ".", ref, ".", alt , sep=""))
            ## damaging based on functional impact for ns
            damFunc <- cbind(
                key <- ns$key,
                sift <- data.frame(sapply(strsplit(as.character(ns$SIFT_pred), ";"), function(x) "D" %in% x)),
                pp2.div <- data.frame(sapply(strsplit(as.character(ns$Polyphen2_HDIV_score), ";"), function(x) length(which(as.numeric(x)>0.5))>0)),
                pp2.var <- data.frame(sapply(strsplit(as.character(ns$Polyphen2_HVAR_score), ";"), function(x) length(which(as.numeric(x)>0.5))>0)),
                lrt <- data.frame(sapply(strsplit(as.character(ns$LRT_pred), ";"), function(x) "D" %in% x)),
                mt <- data.frame(sapply(strsplit(as.character(ns$MutationTaster_converted_score), ";"), function(x) length(which(as.numeric(x)>0.5))>0)), ## "." checked for consistency length(which(as.numeric(".")>0.5))>0 == length(which(NA>0.5))>0 == FALSE
                ma  <- data.frame(sapply(strsplit(as.character(ns$MutationAssessor_converted_score), ";"), function(x) length(which(as.numeric(x)>0.65))>0)),
                ftm <- data.frame(sapply(strsplit(as.character(ns$FATHMM_converted_score), ";"), function(x) length(which(as.numeric(x)>=0.45))>0))
            )
            colnames(damFunc) <- c("key", "sift", "pp2.div", "pp2.var", "lrt", "mt", "ma", "ftm")
            ## I consider as damaging all the mutation that pass at least 5/7 function scores
            damFunc$damagingFunc <- apply(damFunc[,2:length(damFunc)],1, function(x) sum(x) >= 5)
            ## I order it in descending order (because FALSE==0 and TRUE==1) and deduplicate the data frame so only the TRUE (first occurence) for mutations with multiple rows is kept
            damFunc <- damFunc %>% arrange(key, desc(damagingFunc)) %>% subset(!duplicated(key))
            ## Step 2: Infer damaging effect for the nonsynonymous mutations based on conservation
            ## damaging based on functional impact for ns
            ns$phyloP7way_vertebrate <- as.character(ns$phyloP7way_vertebrate)
            ns$GERP.._RS <- as.character(ns$GERP.._RS)
            ns$SiPhy_29way_logOdds <- as.character(ns$SiPhy_29way_logOdds)
            damCons <- cbind(
                key <- ns$key,
                phylop <- data.frame(sapply(strsplit(ns$phyloP7way_vertebrate, ";"), function(x) length(which(as.numeric(x)>1.6))>0)),
                gerp.rs <- data.frame(sapply(strsplit(ns$GERP.._RS, ";"), function(x) length(which(as.numeric(x)>4.4))>0)),
                siphy <- data.frame(sapply(strsplit(ns$SiPhy_29way_logOdds, ";"), function(x) length(which(as.numeric(x)>12.17))>0))
            )
            colnames(damCons) <- c("key", "phylop", "gerp.rs", "siphy")
            ## I consider as damaging all the mutation that pass at least 2/3 conservation scores
            damCons$damagingCons <- apply(damCons[,2:length(damCons)],1, function(x) sum(x) >= 2)
            ## I order it in descending order (because FALSE==0 and TRUE==1) and deduplicate the data frame so only the TRUE (first occurence) for mutations with multiple rows is kept
            damCons <- damCons %>% arrange(key, desc(damagingCons)) %>% subset(!duplicated(key))
            ## Combine the Functional and Conservation scores
            damNS <- damFunc %>% left_join(damCons, by=c("key"))
            ## Include overall damaging prediction based on either 5/7 Functional scores or 2/3 Conservation scores
            damNS$damagingNS <- apply(damNS[,c(9,13)],1, function(x) sum(x) >= 1)
            damNS$key <- as.character(damNS$key)
        }
    }
    ## ********************************
    ## Get dbnsfp results for splicing
    ## ********************************
    ## Step 3: Infer damaging mutations for the splicing mutations
    if(nrow(sp)>0){ ## Sometimes either there are no splicing or no splicing returned from dbNSFP
        sc <- read.delim(scoutfile,stringsAsFactor=F,quote = "")
        if(nrow(sc)>0){
            sc <- sc %>% mutate(key=paste("chr",chr, ".", pos, ".", ref, ".", alt , sep=""))
            sc$ada_score <- as.character(sc$ada_score)
            sc$rf_score <- as.character(sc$rf_score)
            ## damaging based on functional impact for ns
            damSC <- cbind(
                key <- sc$key,
                ada <- data.frame(sapply(strsplit(sc$ada_score, ";"), function(x) length(which(as.numeric(x) >= 0.6))>0)),
                rf <- data.frame(sapply(strsplit(sc$rf_score, ";"), function(x) length(which(as.numeric(x) >= 0.6))>0))
            )
            colnames(damSC) <- c("key", "ada_score", "rf_score")
            ## Consider the mutation as damaging if at least one of the scores predict it as damaging
            damSC$damagingSC <- apply(damSC[,2:length(damSC)],1, function(x) sum(x) >= 1)
            ## I order it in descending order (because FALSE==0 and TRUE==1) and deduplicate the data frame so only the TRUE (first occurence) for mutations with multiple rows is kept
            damSC <- damSC %>% arrange(key, desc(damagingSC)) %>% subset(!duplicated(key))
            damSC$key <- as.character(damSC$key)
        }
    }

    muts <- muts %>% mutate(key=paste(chr, ".", end, ".",ref, ".", alt , sep=""))
    if(nrow(ns)>0){
        muts <- muts %>% left_join(damNS, by=c("key"))
    }else{ ## In case ns is 0 put NAs
        muts = muts %>% mutate(damagingNS=NA)
    }

    if(nrow(sp)>0 & exists("sc")){ ## Sometime either there are no splicing or very few and no splicing return from dbNSFP
      if(nrow(sc)>0){
        muts <- muts %>% left_join(damSC, by=c("key"))
      }else{
        muts = muts %>% mutate(damagingSC=NA)
      }
    }else{ ## Put NAs if there are no splicing mutations
        muts = muts %>% mutate(damagingSC=NA)
    }

    ## Finally define damaging as TRUE or FALSE regardless if it is nonsynonymous or splicing
    muts$damaging <- apply(muts[,c("damagingNS", "damagingSC")],1, function(x) sum(x[!is.na(x)]) >= 1)
    ix <- which(muts$ExonicFunc.refGene %in% c("stopgain","stoploss","stopgain","stoploss", "frameshift deletion", "frameshift insertion", "frameshift substitution"))
    if(length(ix)>0){
        muts$damaging[ix] <- TRUE
    }

    save(muts, file=paste0(save_dir, "/muts_annovar_dbnsfp.Rdata"))
    
    # adding 19014
    load(gene_symbols_fn) ## object name geneSymbols
    ## Check which of the gene symbol are in the symbols of our 19014
    ## First get index of Gene.refGene
    message("Annotating 19,014...")
    muts = get_19014(muts, geneSymbols=geneSymbols)
    return(muts)
    
}

## Feed the mutations from above to the oncodrive clust function
## Here I need to run all the samples together because OncodriveClust does the prediction using recurrent mutations
runOncodriveClust = function(muts=NULL,
                             annovar_path="/mnt/lustre/users/k1469280/FC/Software/annovar_2015_12_14",
                             ucsc_refgene_path=NULL,
                             onco_path = "/mnt/lustre/users/k1469280/FC/Software/OncodriveClust_0.4.1",
                             cgc_phen_path = "/mnt/lustre/users/k1469280/mourikisa/data/CGC_phen.tsv",
                             save_dir=NULL){
    
    if(is.null(muts)){
        stop("No mutations provided")
    }
    
    # === ONCOCLUSTER =====
    ## For OncodriveClust I need 1: Synonymous mutations, 2: Non-Synonymous mutations, 3: Transcript lengths
    message("Running OncodriveClust...")
    # Prepare INPUT files (if no path provided)
    if(is.null(ucsc_refgene_path)){
        ucsc_refgene_path = paste0(save_dir, "/transcript_length.tsv")
        ucsc_refgene <- read.table(paste0(annovar_path, "/humandb/hg19_refGene.txt"), header = F, stringsAsFactors = F)
        colnames(ucsc_refgene) <- c("bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames")
        ## Replace the exon start with the cdsStart to exclude UTRs 5&3 and same for the cdsEnd
        pb <- txtProgressBar(min = 0, max = nrow(ucsc_refgene), style = 3)
        for (row in 1:nrow(ucsc_refgene)){
            ucsc_refgene$exonStarts[row] <- paste(paste(c(ucsc_refgene$cdsStart[row],unlist(strsplit(ucsc_refgene$exonStarts[row], "\\,"))[which(unlist(strsplit(ucsc_refgene$exonStarts[row], "\\,")) > ucsc_refgene$cdsStart[row] & unlist(strsplit(ucsc_refgene$exonStarts[row], "\\,")) < ucsc_refgene$cdsEnd[row])]), collapse=","), ",", sep = "")
            ucsc_refgene$exonEnds[row] <- paste(paste(c(unlist(strsplit(ucsc_refgene$exonEnds[row], "\\,"))[which(unlist(strsplit(ucsc_refgene$exonEnds[row], "\\,")) < ucsc_refgene$cdsEnd[row] & unlist(strsplit(ucsc_refgene$exonEnds[row], "\\,")) > ucsc_refgene$cdsStart[row])], ucsc_refgene$cdsEnd[row]), collapse=","), ",", sep = "")
            setTxtProgressBar(pb, row)
        }
        close(pb)
        ucsc_refgene <- ucsc_refgene %>% select(name2, name, chrom, exonStarts, exonEnds)
        ucsc_refgene <- do.call(rbind, lapply(split(ucsc_refgene, rownames(ucsc_refgene)), function(x) cbind(x[,1],x[,2],x[,3],  unlist(strsplit(x[,4],"\\,")), unlist(strsplit(x[,5],"\\,"))))) %>%
            data.frame(stringsAsFactors=F)
        colnames(ucsc_refgene) <- c("Symbol", "Transcript.id", "chrom", "exon_start", "exon_end")
        ucsc_refgene$exon_start <- as.numeric(ucsc_refgene$exon_start)
        ucsc_refgene$exon_end <- as.numeric(ucsc_refgene$exon_end)
        #ucsc_refgene = ddply(ucsc_refgene, .(Symbol, Transcript.id), mutate, n=1:length(exon_start), .progress = 'text')
        ucsc_refgene <- ucsc_refgene %>% mutate(length=exon_end-exon_start) %>% group_by(Symbol, Transcript.id) %>% summarize(CDS.length=sum(length)) %>% ungroup %>% subset(CDS.length!=0)
        write.table(ucsc_refgene, file=ucsc_refgene_path, quote = F, row.names = F, sep = "\t")
    }
    
    df_mut = as.data.frame(muts, stringsAsFactors=F)
    df_mut = unrowname(df_mut)
    nsyn = subset(df_mut, ExonicFunc.refGene=="nonsynonymous")
    syn = subset(df_mut, ExonicFunc.refGene=="synonymous")
    
    tmp = do.call('rbind', strsplit(sapply(strsplit(nsyn$AAChange.refGene,"\\,"), function(x) x[1] ),"\\:"))
    tmp = data.frame(
        symbol =  tmp[,1],
        Transcript.id =  tmp[,2],
        exon =  tmp[,3],
        nChange =  tmp[,4],
        aaChange =  tmp[,5],
        key = with(nsyn, paste0(chr,'.',start,'.',end,'.',ref,'.',alt, '.', sample)),
        aa.position = gsub("[^0-9]", "", tmp[,5])
    )
    nsyn_fn = paste0(save_dir, "/nsyn.onco")
    write.table(tmp, file=nsyn_fn, row.names = F, quote = F, sep = '\t')
    
    tmp = do.call('rbind', strsplit(sapply(strsplit(syn$AAChange.refGene,"\\,"), function(x) x[1] ),"\\:"))
    tmp = data.frame(
        symbol =  tmp[,1],
        Transcript.id =  tmp[,2],
        exon =  tmp[,3],
        nChange =  tmp[,4],
        aaChange =  tmp[,5],
        key = with(syn, paste0(chr,'.',start,'.',end,'.',ref,'.',alt,'.',sample)),
        aa.position = gsub("[^0-9]", "", tmp[,5])
    )
    syn_fn = paste0(save_dir, "/syn.onco")
    write.table(tmp, file=syn_fn, col.names=T, row.names = F, quote = F, sep = '\t')
    
    # cd /home/FC/Software/OncodriveClust_0.4.1
    # source env/bin/activate
    # module load python
    # oncodriveclust -m 3 -c --cgc /home/mourikisa/data/CGC_phen.tsv ~/OAC/nsyn.onco ~/OAC/syn.onco /home/mourikisa/data/transcript_length.tsv -o ~/OAC/oncodriveclust-results.tsv
    onco_out_fn = paste0(save_dir, "/oncodriveclust-results.tsv")
    message(nsyn_fn)
    message(syn_fn)
    message(onco_out_fn)
    cmd <- paste0("cd ", onco_path, ";source env/bin/activate;oncodriveclust -c -m 5 --cgc ", cgc_phen_path ,
                  " -o ", onco_out_fn, " ", nsyn_fn, " ", syn_fn, " ", ucsc_refgene_path)
    message(noquote(cmd))
    system(noquote(cmd))
    
    message("Parsing the output of OncodriverClust...")
    onco_out <- read.table(onco_out_fn, header = T, sep = "\t")
    onco_out <- onco_out %>% subset(QVALUE<=0.1)
    onco_in <- read.table(nsyn_fn, header = T, sep = "\t")
    onco_in$Sample <- apply(onco_in, 1, function(x) unlist(strsplit(x[6], "\\."))[6])
    
    ## Each cluster a separate row - it will be easier to gather patients and check if clusters are costant across cancer types
    onco_out <- onco_out %>% mutate(CLUST_COORDS=strsplit(CLUST_COORDS, "\\,\\[")) %>%
        unnest(CLUST_COORDS) %>% mutate(CLUST_COORDS=gsub("\\[|\\]", "", CLUST_COORDS)) %>%
        separate(CLUST_COORDS, into=c("CLUST_COORDS_START", "CLUST_COORDS_END"), sep="\\,") %>%
        separate(CLUST_COORDS_END, into=c("CLUST_COORDS_END", "NUMBER_OF_MUTS_IN_CLUST"), sep="\\:")
    
    ## Find which samples have mutations in the clusters
    find_samples <- function(gene, start, end, onco_in){
        keys <- onco_in %>% subset(symbol==gene & as.numeric(aa.position) >= start & as.numeric(aa.position) <= end) %>% select(key)
        keys <- paste(keys$key, collapse=",")
        return(keys)
    }
    onco_out$KEY <- apply(onco_out, 1, function(x) find_samples(x[1], as.numeric(x["CLUST_COORDS_START"]), as.numeric(x["CLUST_COORDS_END"]), onco_in))
    
    ## Expand the key components in the table
    onco_out <- onco_out %>% mutate(KEY=strsplit(KEY, "\\,")) %>% unnest(KEY)
    onco_out$SAMPLE <- apply(onco_out, 1, function(x) unlist(strsplit(x[14], "\\."))[6])
    
    df_mut$key = with(df_mut, paste0(chr,'.',start,'.',end,'.',ref,'.',alt,'.', sample))
    df_mut$oncodriveClust = df_mut$key%in%onco_out$KEY
    return(df_mut)
}


## Get VEP annotation and check the numbers with VEP
checkANNOVAR = function(annovar = muts,
                        muts_dir="/Volumes/mourikisa/data/OAC/129_OAC/129_raw_data/strelka/"){
    
    muts = do.call(rbind, annovar)
    ## Get coding(exonic) mutations from ANNOVAR
    samples2exonic = muts %>% group_by(sample) %>% count(Func.refGene) %>% mutate(Func.refGene=ifelse(Func.refGene=="exonic", "exonic", "other")) %>% 
            group_by(sample, Func.refGene) %>% summarise(annovar=sum(n)) %>% ungroup() %>% subset(Func.refGene=="exonic")
    ## Get non silent
    samples2ns = muts %>% group_by(sample) %>% count(ExonicFunc.refGene) %>% mutate(ExonicFunc.refGene=ifelse(ExonicFunc.refGene%in%ns, "nonsilent", "other")) %>%
        group_by(sample, ExonicFunc.refGene) %>% summarise(annovar=sum(n)) %>% ungroup() %>% subset(ExonicFunc.refGene=="nonsilent")
    
    vep_numbers = NULL
    for (s in unique(muts$sample)){
        cat(s, "\n")
        ## Get the same number from VEP
        dirs = list.files(muts_dir)
        dirs = dirs[grepl(s, dirs)]
        vep_fn=paste0(muts_dir, "/", dirs, "/1.3/strelka/vep/", dirs, ".snp.pass.coding.supplemented.vep")
        vep = read.delim(vep_fn, header = T, sep="\t")
        vep.nonsil <- vep[which(grepl("missense",vep$Consequence) |
                                            grepl("nonsense",vep$Consequence) |
                                            grepl("stop_gained", vep$Consequence) |
                                            grepl("stop_lost", vep$Consequence) |
                                            grepl("splice_donor_variant", vep$Consequence) |
                                            grepl("splice_acceptor_variant", vep$Consequence) |
                                            grepl("splice_region_variant", vep$Consequence) |
                                            grepl("initiator_codon_variant", vep$Consequence)),]
        #vep = vep %>% mutate(Consequence=strsplit(Consequence, ",")) %>% unnest
        ## VEP's output is one transcript per line
        #vep = vep %>% mutate(key=paste0(sample_id, ".", Uploaded_variation))
        
        no.vep.nonsil = vep.nonsil %>% select(Uploaded_variation) %>% unique %>% nrow
        d = data.frame(sample=s, vep=no.vep.nonsil)
        vep_numbers = rbind(vep_numbers, d)
    }
    
    #mut_numbers = samples2exonic %>% left_join(vep_numbers)
    mut_numbers = samples2ns %>% left_join(vep_numbers)
    return(mut_numbers)
}

## Muts argument corresponds to the annotated SNVs and INDELs all together
summariseMuts = function(snvs=snv, indel=indel,muts=muts){
    
    ## This function was created to get the numbers of mutations
    
    ## Get SNVs per patient
    l = lapply(snv, function(x) c(unique(x$sample), nrow(x)))
    summary.table = t(as.data.frame(l)) %>% unrowname() %>% data.frame()
    colnames(summary.table) = c("sample", "SNVs")
    summary.table$SNVs = as.numeric(summary.table$SNVs)
    
    ## Get INDELs per patient
    l = lapply(indel, function(x) c(unique(x$sample), nrow(x)))
    st = t(as.data.frame(l)) %>% unrowname() %>% data.frame()
    colnames(st) = c("sample", "INDELs")
    st$INDELs = as.numeric(st$INDELs)
    summary.table = summary.table %>% left_join(st)
    
    ## Get mutations annotated forchecking purposes
    l = lapply(muts, function(x) c(unique(x$sample), nrow(x)))
    st = t(as.data.frame(l)) %>% unrowname() %>% data.frame()
    colnames(st) = c("sample", "Mutations_annotated")
    st$Mutations_annotated = as.numeric(st$Mutations_annotated)
    summary.table = summary.table %>% left_join(st)
    
    ## Get silent/nonsilent mutations
    l = lapply(muts, function(x){
        x=x%>%subset(!ExonicFunc.refGene%in%ns)
        c(unique(x$sample), nrow(x)) 
    })
    st = t(as.data.frame(l)) %>% unrowname() %>% data.frame()
    colnames(st) = c("sample", "Silent_mutations")
    st$Silent_mutations = as.numeric(st$Silent_mutations)
    summary.table = summary.table %>% left_join(st)
    
    l = lapply(muts, function(x){
        x=x%>%subset(ExonicFunc.refGene%in%ns)
        c(unique(x$sample), nrow(x)) 
    })
    st = t(as.data.frame(l)) %>% unrowname() %>% data.frame()
    colnames(st) = c("sample", "Nonsilent_mutations")
    st$Nonsilent_mutations = as.numeric(st$Nonsilent_mutations)
    summary.table = summary.table %>% left_join(st)
    
    ## Get number of genes in the nonsilent
    l = lapply(muts, function(x){
        x=x%>%subset(ExonicFunc.refGene%in%ns)
        c(unique(x$sample), x%>%select(Gene.refGene)%>%unique%>%nrow) 
    })
    st = t(as.data.frame(l)) %>% unrowname() %>% data.frame()
    colnames(st) = c("sample", "Nonsilent_genes")
    st$Nonsilent_genes = as.numeric(st$Nonsilent_genes)
    summary.table = summary.table %>% left_join(st)
    
    ## Get number of 19014 in the nonsilent
    l = lapply(muts, function(x){
        x=x%>%subset(ExonicFunc.refGene%in%ns)
        c(unique(x$sample), x%>%select(entrez_19014)%>%unique%>%nrow) 
    })
    st = t(as.data.frame(l)) %>% unrowname() %>% data.frame()
    colnames(st) = c("sample", "Nonsilent_genes_19014")
    st$Nonsilent_genes_19014 = as.numeric(st$Nonsilent_genes_19014)
    summary.table = summary.table %>% left_join(st)
    
    ## Get damaging/nondamaging mutations
    l = lapply(muts, function(x){
        x=x%>%subset(ExonicFunc.refGene%in%ns & damaging==FALSE)
        c(unique(x$sample), nrow(x)) 
    })
    st = t(as.data.frame(l)) %>% unrowname() %>% data.frame()
    colnames(st) = c("sample", "Nondamaging_mutations")
    st$Nondamaging_mutations = as.numeric(st$Nondamaging_mutations)
    summary.table = summary.table %>% left_join(st)
    
    l = lapply(muts, function(x){
        x=x%>%subset(damaging==TRUE)
        c(unique(x$sample), nrow(x)) 
    })
    st = t(as.data.frame(l)) %>% unrowname() %>% data.frame()
    colnames(st) = c("sample", "Damaging_mutations")
    st$Damaging_mutations = as.numeric(st$Damaging_mutations)
    summary.table = summary.table %>% left_join(st)
    
    ## Get genes in damaging mutations
    l = lapply(muts, function(x){
        x=x%>%subset(damaging==TRUE)
        c(unique(x$sample), x%>%select(Gene.refGene)%>%unique%>%nrow) 
    })
    st = t(as.data.frame(l)) %>% unrowname() %>% data.frame()
    colnames(st) = c("sample", "Damaging_genes")
    st$Damaging_genes = as.numeric(st$Damaging_genes)
    summary.table = summary.table %>% left_join(st)
    ## And the 19014
    l = lapply(muts, function(x){
        x=x%>%subset(damaging==TRUE)
        c(unique(x$sample), x%>%select(entrez_19014)%>%unique%>%nrow) 
    })
    st = t(as.data.frame(l)) %>% unrowname() %>% data.frame()
    colnames(st) = c("sample", "Damaging_genes_19014")
    st$Damaging_genes_19014 = as.numeric(st$Damaging_genes_19014)
    summary.table = summary.table %>% left_join(st)
    
    stats = sapply(summary.table[,2:length(summary.table)], function(x) summary(x))
    
    return(list(summary.table=summary.table, stats=stats))
}

summariseCNVs = function(cnv=cnvs){
    
    ## This function was created to get the numbers of CNV data
    cnv  = cnv[["df_cnv_19014"]]
    
    ## Annotate Gains and Losses with our previous definition 
    cnv = cnv %>% mutate(Gain_old=ifelse(Total_CN>=4, 1, 0), Loss_old=ifelse(Total_CN==1 | Total_CN==0, 1, 0))
    ## Annotate Gains and Losses using ploidy
    cnv = cnv %>% mutate(Gain_ploidy=ifelse(Total_CN>=(2*ploidy), 1, 0), Loss_ploidy=ifelse(Total_CN<=(0.5*ploidy), 1, 0))
    ## Define Hetero/Hom Losses    
    cnv = cnv %>% mutate(Hetero_Loss=ifelse(Total_CN==1, 1, 0), Homo_Loss=ifelse(Total_CN==0, 1, 0))
    
    summary.table = cnv %>% subset(Gain_old==1) %>% group_by(sample) %>% summarise(Gains_old_genes=length(unique(entrez_19014)))
    st = cnv %>% subset(Loss_old==1) %>% group_by(sample) %>% summarise(Loss_old_genes=length(unique(entrez_19014)))
    summary.table = summary.table %>% left_join(st)    
    
    st = cnv %>% subset(Gain_ploidy==1) %>% group_by(sample) %>% summarise(Gain_ploidy_genes=length(unique(entrez_19014)))
    summary.table = summary.table %>% left_join(st)
    st = cnv %>% subset(Loss_ploidy==1 & (Hetero_Loss==1 | Homo_Loss==1)) %>% group_by(sample) %>% summarise(Loss_ploidy_genes=length(unique(entrez_19014)))
    summary.table = summary.table %>% left_join(st)
    
    st = cnv %>% subset(Loss_ploidy==1 & Hetero_Loss==1) %>% group_by(sample) %>% summarise(Loss_ploidy_genes_Hetero=length(unique(entrez_19014)))
    summary.table = summary.table %>% left_join(st)
    st = cnv %>% subset(Loss_ploidy==1 & Homo_Loss==1) %>% group_by(sample) %>% summarise(Loss_ploidy_genes_Homo=length(unique(entrez_19014)))
    summary.table = summary.table %>% left_join(st)
    
    summary.table[is.na(summary.table)] = 0
    stats = sapply(summary.table[,2:length(summary.table)], function(x) summary(x))
    
    return(list(summary.table=summary.table, stats=stats))   
}


summariseSVs = function(svs=svs){
    
    ## This function was created to get the numbers of SV data
    
    summary.table = svs[,1:7] %>% gather(type, value, -Sample, -gene) %>% subset(value!=0) %>% group_by(Sample, type) %>% summarise(n=n()) %>% ungroup() %>% spread(type, n)
    
    svs[,c(1:7, 9)] %>% subset(!is.na(entrez_19014)) %>% select(-entrez_19014) %>% gather(type, value, -Sample, -gene) %>% subset(value!=0) %>% 
        group_by(Sample, type) %>% summarise(n=n()) %>% ungroup() %>% group_by(type) %>% summarise(min=min(n), median=median(n), mean=mean(n), max=max(n), n=n())
    
    svs[,c(1:7, 10)] %>% subset(cancer_type=="cgc" | cancer_type=="can") %>% select(-cancer_type) %>% gather(type, value, -Sample, -gene) %>% subset(value!=0) %>% 
        group_by(Sample, type) %>% summarise(n=n()) %>% ungroup() %>% group_by(type) %>% summarise(min=min(n), median=median(n), mean=mean(n), max=max(n), n=n())
    
}

## Get those in 19014, refine indels as we did for TCGA
applyMutationFilters = function(df=NULL){
    
    ## Part of the 19014 code here
    
    
    if(is.null(df)){
        return(NULL)
    }else{
        ## Apply the same filters as for nonsynonymous
        df <- subset(df, MeanMutFreq >= 0.10 | is.na(MeanMutFreq) | is.nan(MeanMutFreq))
        df <- subset(df, !(Center == "hgsc.bcm.edu" & MeanMutFreq == 1) | is.na(MeanMutFreq) | is.nan(MeanMutFreq))
        df <- subset(df, !(t_alt_count == 2 & t_ref_count == 4) | is.na(t_alt_count) | is.nan(t_alt_count) | is.na(t_ref_count) | is.nan(t_ref_count))
        return(df)
    }
}

annotate_genes = function(x, gene_coord_fn=NULL, save_dir=NULL){
    
    if(is.null(gene_coord_fn) | is.null(save_dir)){
        stop("annotate_genes: parameters missing")
    }
    
    load(gene_coord_fn)
    ## Use only genes with mainEntrez==TRUE
    coord <- subset(coord, mainEntrez==TRUE)
    coord_file <- paste0(save_dir, "/tmp.blat.genes.bed")
    
    write.table(unique(coord[,c("chrom","start","end","symbol","Entrez")]),file=coord_file,row.names=F,col.names = F,quote = F,sep="\t")
    write.table(x[,c("Chromosome","Start","End")], file=paste0(save_dir, "/tmp.bed"), row.names=F,col.names = F,quote = F,sep="\t")
    
    #system("intersectBed -a tmp.blat.genes.bed -b tmp.bed -f 0.25 -wa -wb > tmp.genes.bed", wait = TRUE)
    system(paste0("intersectBed -a ", coord_file, " -b ", paste0(save_dir, "/tmp.bed"), " -f 0.25 -wa -wb > ", paste0(save_dir, "/tmp.genes.bed")), wait = TRUE)
    tmp = read.table(paste0(save_dir, "/tmp.genes.bed"), h=F)
    tmp$key=paste0(tmp$V6,".",tmp$V7,'.',tmp$V8,'.',tmp$V9)
    key=paste0(x$Chromosome,".",x$Start,'.',x$End,'.',x$sample)
    y = cbind(tmp[,1:5], x[match(tmp$key, key),])
    colnames(y)[1:5] = c("chrom","start","end","symbol_19014","entrez_19014")
    y
}

getCNVs = function(cnv_fn=NULL, stats_fn=NULL,
                   g=NULL,
                   gene_coord_fn="/mnt/lustre/users/k1469280/mourikisa/data/blat_exon_coordinates.Rdata",
                   save_dir=NULL){
    
    
    
    if(is.null(cnv_fn) | is.null(stats_fn) | is.null(g) | is.null(save_dir)){
        stop("Please provide file name")
    }
    
    chr_sizes_fn=paste0("/mnt/lustre/users/k1469280/FC/DB/", g, "/", g, ".chrom.sizes")
    message("Reading ASCAT files...")

    x =  read.csv(cnv_fn)
    x$Chromosome = paste0("chr",x$Chromosome)
    x$Major_CN = with(x, Total_CN-Minor_CN)
    x = x[,c('Chromosome','Start','End', 'Total_CN', 'Major_CN', 'Minor_CN' )]
    x$CNV_type = NA
    x$CNV_type[ which(x$Total_CN>2) ] =  "Gain"
    x$CNV_type[ which(x$Total_CN<2) ] =  "Loss"
    x$CNV_type[ which(x$Total_CN==2 & x$Minor_CN==0) ] =  "Loss"
    x$real_CNV_type = NA
    x$real_CNV_type[ which(x$Total_CN>2) ] =  "Gain"
    x$real_CNV_type[ which(x$Total_CN<2) ] =  "Loss"
    x$real_CNV_type[ which(x$Total_CN>=2 & x$Minor_CN==0) ] =  "LOH"
    ploidy = read.csv( stats_fn, sep=" ", header = F)
    colnames(ploidy) = c("type", "value")
    ploidy = round(ploidy$value[ploidy$type=="Ploidy"], digits = 2)
    x$ploidy = ploidy
        
    
    chr = read.delim(chr_sizes_fn, h=F)
    colnames(chr) = c("Chromosome", "chrSize")
    x$region_len = with(x, End-Start)
    print(nrow(x)) ## to check joining
    x = x %>% left_join(chr) %>% data.frame()
    print(nrow(x)) ## to check joining
    
    message("Running gene intersection...")
    x_19014 = annotate_genes(x, gene_coord_fn=gene_coord_fn, save_dir=save_dir)
    
    return(list(df_cnvs=x, df_cnvs_19014=x_19014))
    
}

getSVs = function(sv_fn=NULL,
                  geneSymbols_fn="/home/mourikisa/data/19014_GeneSymbols.Rdata",
                  geneInfo_fn="/home/mourikisa/data/geneInfoNCG5.Rdata",
                  cancerGenes_fn="/home/mourikisa/data/cancerGenesNCG5.Rdata"){
    
    if(is.null(sv_fn)){
        stop("Please provide file name")
    }
    message("Reading Manta calls...")
    svs = read.delim(sv_fn, header = T)
    
    ## Get table with SVs per patient and gene
    res = NULL
    for (t in unique(svs$Type)){
        d = svs %>% subset(Type==t) %>% select(Type, GeneName.node1, GeneName.node2) %>% gather(node, gene, -Type) %>% subset(gene!="-") %>% select(-node, -Type)
        if(nrow(d)==0){ ## in case a certain type of SV has 0 hits in a sample
            res[,t] = NA
        }else{
            d = d %>% unique
            d[,t] = 1
            if(is.null(res)){
                res = rbind(res, d)
            }else if(!is.null(res)){
                res = res %>% full_join(d)
            }
        }
    }
    
    ## Add Fusions - Fusions is a subcategory and can be present in any of the 5 main SV types
    ## For now I am not sure how to parse fusions (i.e fusions+ versus fusions- and which alteration affects what?)
    ##svs %>% subset(grepl("FUSION", FusionPrediction.node1) | grepl("FUSION", FusionPrediction.node2)) %>% head
    
    
    res[is.na(res)] = 0
    
    ## get the annotation of the 19014
    load(geneSymbols_fn)
    res = res %>% mutate(symbol_19014=ifelse(gene %in% geneSymbols$symbol, gene, NA)) %>% 
        left_join(geneSymbols%>%select(symbol, Entrez)%>%rename(symbol_19014=symbol, entrez_19014=Entrez))
    
    
    ## Add annotation for the 19014 and for NCG5
    load(geneInfo_fn)
    load(cancerGenes_fn)
    ## Fix gene info table from NCG
    geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
    ## Get a cancer gene with all the associated primary sites and cancer sites
    cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>% 
        group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","), 
                                       cancer_site=paste(unique(cancer_site), collapse=",")) %>%
        ungroup
    geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))
    
    ## Mark genes that are in NCG
    res = res %>% left_join(geneInfo%>%rename(entrez_19014=entrez))
    return(res)
}

## This function combines all types of data to single table 
## and prepares them for the extraction of drivers and prediction
createTotalTable = function(muts=NULL, cnvs=NULL, svs=NULL, exclude_samples=NULL){
    
    ns = c("nonsynonymous","stopgain","frameshift deletion","splicing","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","stoploss","frameshift substitution")
    dam = c("nonsynonymous","frameshift deletion","frameshift insertion","frameshift substitution","splicing","stopgain","stoploss")
    trunc = c("frameshift deletion","frameshift insertion","frameshift substitution","stopgain","stoploss") ## Always damaging==TRUE
    non_trunc = c("nonsynonymous","splicing")
    
    ## Make the lists
    message("Integrating SNVs...")
    df_mut = muts
    rm(muts)
    
    if(!is.null(exclude_samples)){
        df_mut = df_mut %>% subset(!sample%in%exclude_samples)
        message(paste0("Samples excluded in mutation data: ", paste0(exclude_samples, collapse=",")))
    }
    
    ## Fix nonsilent here
    df_mut = df_mut %>% select(-nonsilent)
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
    df_mut = df_mut %>% subset(!is.na(entrez_19014))
    
    ## Create the total table
    total_muts = ddply(df_mut, .(sample, symbol_19014, entrez_19014), summarise,
                       no_ALL_muts=n(),
                       no_NSI_muts=sum(nonsilent),
                       no_TRUNC_muts = sum(ExonicFunc.refGene %in% trunc),
                       no_NTDam_muts = sum(ExonicFunc.refGene %in% non_trunc & damaging),
                       no_GOF_muts = sum(oncodriveClust), .progress = 'text'
    )
    
    ## Add protein position if needed
    #aa = df_mut %>% select(sample, entrez_19014, symbol_19014, AAChange.refGene) %>% group_by(sample, entrez_19014, symbol_19014) %>% summarise(AAChange=paste(unique(AAChange.refGene), collapse=","))
    #total_muts = total_muts %>% left_join(aa)
    
    # ## Check that you see a difference in the number of total mutations
    # geneInfo_fn="/Volumes/mourikisa/data/geneInfoNCG5.Rdata"
    # cancerGenes_fn="/Volumes/mourikisa/data/cancerGenesNCG5.Rdata"
    # load(geneInfo_fn)
    # load(cancerGenes_fn)
    # ## Fix gene info table from NCG
    # geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
    # ## Get a cancer gene with all the associated primary sites and cancer sites
    # cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
    #     group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
    #                                    cancer_site=paste(unique(cancer_site), collapse=",")) %>%
    #     ungroup
    # geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))
    # test = total_muts %>% left_join(geneInfo, by=c("entrez_19014"="entrez"))
    # 
    # test = test %>% mutate(sumDrivers = rowSums(.[6:8]))
    # test = test %>% mutate(dVa=sumDrivers/no_ALL_muts)
    # test %>% subset(dVa==1) %>% head
    # wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$dVa, test%>%subset(cancer_type=="can")%>%.$dVa)
    # summary(test%>%subset(cancer_type=="cgc")%>%.$dVa)
    # summary(test%>%subset(cancer_type=="can")%>%.$dVa)
    
    
    ## Bring in the total_table the SVs
    message("Integrating SVs...")
    
    if(!is.null(exclude_samples)){
        svs = svs %>% subset(!sample%in%exclude_samples)
        message(paste0("Samples excluded in SV data: ", paste0(exclude_samples, collapse=",")))
    }
    
    if(!is.null(svs)){
        ## There are 2 genes duplicated in 4 samples due to aliases
        svs = svs %>% subset(!is.na(entrez_19014)) %>% mutate(key=paste(sample, entrez_19014, sep=".")) %>% subset(!duplicated(key)) %>% select(-key)
        svs = svs %>% subset(BND>0 | INS>0 | INV>0) %>% select(-gene, -cancer_type, -primary_site, -cancer_site)
        ## And put them in the total table as well
        total_table = total_muts %>% full_join(svs%>%select(sample, BND, INS, INV, entrez_19014)%>%subset(!is.na(entrez_19014)))
        total_table = total_table %>% mutate(key=paste(sample, entrez_19014, sep="."))
    }else{
        total_table = total_muts %>% mutate(BND=0, INS=0, INV=0)
        total_table = total_table %>% mutate(key=paste(sample, entrez_19014, sep="."))
    }

    
    
    
    message("Integrating CNVs...")
    
    ## Add also the genes that are in muts and SVs to get their ploidy and Copy number
    cnvs = cnvs %>% mutate(key=paste(sample, entrez_19014, sep="."))
    df_cnv = cnvs %>% subset(key%in%total_table$key | !is.na(CNV_type_corrected)) ## Also get the real CNVs

    
    if(!is.null(exclude_samples)){
        df_cnv = df_cnv %>% subset(!sample%in%exclude_samples)
        message(paste0("Samples excluded in CNV data: ", paste0(exclude_samples, collapse=",")))
    }
    
    ## define Gains and Losses - this was done in previous step
    ## Deduplicate the CNV data, because some genes may fall into two regions (sometimes it can be gain and loss)
    message("Resolving duplicated entries in CNVs...")
    dups = df_cnv %>% group_by(key) %>% mutate(n=n(), types=paste(unique(CNV_type_corrected), collapse=","), ntypes=length(unique(CNV_type_corrected))) %>% subset(n>1)
    dups = dups %>% ungroup()
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
            d = dups %>% subset(types==t) %>% arrange(key, desc(CNV_type_corrected)) %>% subset(!duplicated(key))
            dups_refined = rbind(dups_refined, d)
        }else if (t=="Gain" | t=="NA" | t=="Loss"){ ## Choose the one with the highest overlap
            d = dups %>% subset(types==t) %>% mutate(overlap=(overlap(start, end, Start, End)/(end-start))*100) %>% arrange(key, desc(overlap)) %>% subset(!duplicated(key)) %>% select(-overlap)
            dups_refined = rbind(dups_refined, d)
        }else if (grepl("Gain", t) & grepl("Loss", t)){ ## For those genes we have both gain and loss, I set CNV_type to NA because we cannot distinguish between the two
            d = dups %>% subset(types==t) %>% mutate(CNV_type_corrected=NA) %>% arrange(key) %>% subset(!duplicated(key))
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
        full_join(df_cnv%>%select(sample, entrez_19014, Total_CN, CNV_type_corrected, ploidy, n)%>%rename(CNV_entries=n)%>%subset(!is.na(entrez_19014)))
    
    
    total_table$na_19014 = apply(total_table[,c("symbol_19014", "entrez_19014")], 1, function(x) length(x[is.na(x)]))
    total_table = total_table %>% mutate(in_19014=ifelse(na_19014<2, TRUE, FALSE)) %>% select(-na_19014) %>% subset(in_19014==TRUE)
    
    return(total_table)
    
}

## Give to this function the output of the createTotalTable function above
getMLinput <- function(df, geneProperties_dir="Mountpoints/rosalind_lustre/mourikisa/data/geneProperties_final_mmImputed.Rdata"){
    
    df = df %>% mutate(Cancer_type="OAC") %>% rename(Sample=sample, Entrez=entrez_19014, Copy_number=Total_CN, CNV_type=CNV_type_corrected) %>% select(-symbol_19014, -in_19014, -ploidy, -CNV_entries, -key)
    
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
    
    message("Fixing SVs...")
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


getMutationDetails = function(mainDir = "~/rosalind_lustre/mourikisa/data/OAC/87_OAC/66_ICGC/"){
    sample_dirs = list.dirs(paste0(mainDir, "/strelka"), recursive = F)
    total_ann = data.frame()
    for(s in sample_dirs){
        cat(s, "\n")
        ann_fn = paste0(s, "/parsing_and_annotation/annovar/annovar.snv.out.hg19_multianno.txt")
        ann = read.table(ann_fn, sep="\t", header = T)
        total_ann = rbind(total_ann, ann)
    }
    return(total_ann)
}








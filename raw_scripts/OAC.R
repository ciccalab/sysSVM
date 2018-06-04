source("config.R")

#***********************************
#           MUTATIONS
#***********************************


# === SNV ======

a = list.files(path="strelka/")
a = as.list(a)

r = read.table("cufflinks/rna_seq_mapping.txt",sep = "\t")
mapping=as.data.frame(do.call("rbind", strsplit(unlist(a),"_vs_")))
colnames(mapping) = c('Tumor','Control')
mapping$RNAseq=mapping[,1]%in%r[,2]
save(mapping, file="Rdata/env.Rdata")

l = lapply(a, function(x) readVcf( paster("strelka/",x,"/1.3/strelka/filtered_results/",x,".snp.pass.vcf"),genome="hg19" ))
names(l) = paster(mapping[,1], " (", sapply(l,nrow),")")

snv = list()
for(i in 1:length(l)){
  print(i)
  x = l[[i]]
  snv[[i]] = data.frame(
              chr= as.character(seqnames(rowData(x))),
              start=start(ranges(rowData(x))),
              end=end(ranges(rowData(x))),
              ref=as.character(rowData(x)$REF),
              alt=sapply(strsplit(names(rowData(x)[,1]), "\\/"), function(x) x[2]),
              freq = info(x)$VariantAlleleFrequency,
              ReadCount = info(x)$ReadCount,
              VariantAlleleCount = info(x)$VariantAlleleCount,
              ReadCountControl = info(x)$ReadCountControl,
              VariantAlleleCountControl = info(x)$VariantAlleleCountControl,
              VariantStrandBias = info(x)$VariantStrandBias
  )
}
names(snv) = names(l)
save(snv,l, file="snvs.Rdata")

ann = unique(do.call("rbind", snv)[,1:5])
ann = unrowname(ann)
write.table(ann, file="annovar.snv", col.names=F, row.names=F, quote=F, sep="\t")

# === INDEL ====

a = list.files(path="strelka/")
a = as.list(a)

s = lapply(a, function(x) readVcf( paster("strelka/",x,"/1.3/strelka/results/passed.somatic.indels.vcf"),genome="hg19" ))
names(s) = mapping[,1]


indel = list()
for(i in 1:length(s)){
  print(i)
  x = s[[i]]
  indel[[i]] = data.frame(
    chr= as.character(seqnames(rowData(x))),
    start=start(ranges(rowData(x))),
    end=end(ranges(rowData(x))),
    ref=as.character(rowData(x)$REF),
    alt=sapply(strsplit(names(rowData(x)[,1]), "\\/"), function(x) x[2]),
    freq = (geno(x)$TAR[,2,1])/(geno(x)$DP[,"TUMOR"]),
    ReadCount = geno(x)$DP[,"TUMOR"],
    VariantAlleleCount = geno(x)$TAR[,2,1],
    ReadCountControl = geno(x)$DP[,"NORMAL"],
    VariantAlleleCountControl = NA,
    VariantStrandBias = NA
  )
}

indel = lapply(indel, function(x) subset(x, freq>=0.10))
indel = lapply(indel, function(x){ x$freq[which(x$freq>=1)]=1; return(x)})

ann = unique(do.call("rbind", indel)[,1:5])
ann = unrowname(ann)
write.table(ann, file="annovar.indel", col.names=F, row.names=F, quote=F, sep="\t")


ai(file="allele_freq_snv.pdf", w=10,h=10)
par(mfrow=c(5,4))
for(i in 1:18) hist(info(l[[i]])$VariantAlleleFrequency, main=names(l)[i], xlab = "Allele Freq")
dev.off()

ai(file="allele_freq_indels.pdf", w=10,h=10)
par(mfrow=c(5,4))
for(i in 1:19) hist(indel[[i]]$freq, main=names(s)[i], xlab = "Allele Freq")
dev.off()

# === ANNOVAR ====
# cd /home/FC/Software/annovar_2015_12_14
# perl table_annovar.pl -buildver hg19 -protocol refGene,
# perl table_annovar.pl ~/OAC/annovar.snv humandb/ -buildver hg19 -out ~/OAC/annovar.snv.out -remove -protocol refGene,esp6500siv2_all,esp6500siv2_ea,1000g2015aug_all,1000g2015aug_eas,snp138  -operation g,f,f,f,f,f
# perl table_annovar.pl ~/OAC/annovar.indel humandb/ -buildver hg19 -out ~/OAC/annovar.indel.out -remove -protocol refGene,esp6500siv2_all,esp6500siv2_ea,1000g2015aug_all,1000g2015aug_eas,snp138  -operation g,f,f,f,f,f

as = read.delim(file="annovar.snv.out.hg19_multianno.txt")
ai = read.delim(file="annovar.indel.out.hg19_multianno.txt")

ia = with(as,paster( Chr,'.',Start,'.',End,'.',Ref,'.',Alt))
asnv=list()
for(i in 1:length(snv)){
  is = with(snv[[i]],paster( chr,'.',start,'.',end,'.',ref,'.',alt))
  asnv[[i]] = cbind(snv[[i]],as[match(is,ia),c("Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","esp6500siv2_all","esp6500siv2_ea","X1000g2015aug_all","X1000g2015aug_eas","snp138")])
}

ia = with(ai,paster( Chr,'.',Start,'.',End,'.',Ref,'.',Alt))
aindel=list()
for(i in 1:length(indel)){
  is = with(indel[[i]],paster( chr,'.',start,'.',end,'.',ref,'.',alt))
  aindel[[i]] = cbind(indel[[i]],ai[match(is,ia),c("Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","esp6500siv2_all","esp6500siv2_ea","X1000g2015aug_all","X1000g2015aug_eas","snp138")])
}


lapply(asnv, function(x) table(x$ExonicFunc.refGene))
lapply(asnv, function(x) table(x$Func.refGene))

names(asnv)=names(aindel)=mapping[,1]
muts = mapply(rbind, asnv, aindel, SIMPLIFY = F)

muts=lapply(muts, fix_splicing)
muts=lapply(muts, fix_exonic)
muts=lapply(muts, fix_exonicFunc)
muts=lapply(muts, is_nonsilent)
muts=lapply(muts, set_IGV_code)

names(muts) = mapping[,1]
muts =muts[subset(mapping, RNAseq)[,1]]
save(muts, file="Rdata/dataset_mutations.Rdata")


x = lapply(muts, function(x) table(x$Func))
n = unique(unlist(lapply(x, names)))
m = matrix(NA, nr=length(n), nc=18, dimnames = list(n,names(x)))
for(i in names(x)){
  m[names(x[[i]]),i] = x[[i]]
}

x = lapply(muts, function(x) table(x$Exonic))
n = unique(unlist(lapply(x, names)))
m2 = matrix(NA, nr=length(n), nc=18, dimnames = list(n,names(x)))
for(i in names(x)){
  m2[,i] = x[[i]][n]
}

write.xlsx(m, file="Results/Variant_analysis.xlsx","stats",row.names=T, showNA = F)
write.xlsx(m2, file="Results/Variant_analysis.xlsx","stats ex",row.names=T, showNA = F, append=T)

x = cbind(as.matrix(s(sapply(muts, function(x) mutation_frequency(subset(x,nonsilent))))),as.matrix(s(sapply(muts, mutation_frequency, target=3000))))
colnames(x) =c("nonsilent","WG")
write.xlsx(x, file="Results/Variant_analysis.xlsx","mutation freq",row.names=T, showNA = F, append=T)

# === DBNSFP ====

ns = lapply(muts, function(x) subset(x, ExonicFunc.refGene%in%c("nonsynonymous SNV","nonsynonymous")))
df <- unique(do.call("rbind",ns)[,c("chr","end","ref","alt")])
write.table(df, file="ns.oac", col.names=F, row.names=F, quote=F, sep="\t")

sp = lapply(muts, function(x) subset(x, Func.refGene%in%c("splicing")))
df <- unique(do.call("rbind",sp)[,c("chr","end","ref","alt")])
write.table(df, file="sp.oac", col.names=F, row.names=F, quote=F, sep="\t")


# cd /home/FC/DB/dbNSFPv3.0b2a/
# java -Xmx4g search_dbNSFP30b2a  -i ~/OAC/ns.oac -o ~/OAC/ns.oac.dbnsfp -v hg19
# java search_dbNSFP30b2a  -i ~/OAC/sp.oac -o ~/OAC/sp.oac.dbnsfp -s ~/OAC/sp.oac.dbnsfp.sc -v hg19

ns = read.delim("ns.oac.dbnsfp")

## Add MutationTaster_converted_score column
mt_convert <- sapply(1:nrow(ns), function(x) cbind( unlist(strsplit(ns$MutationTaster_score[x], ";")), unlist(strsplit(ns$MutationTaster_pred[x], ";"))) )
ns$MutationTaster_converted_score <- unlist(lapply(mt_convert, MutationTaster_convert))

## Add MutationAssessor_converted_score column
ma_convert <- sapply(1:nrow(ns), function(x) unlist(strsplit(ns$MutationAssessor_score[x], ";")))
ns$MutationAssessor_converted_score <- unlist(lapply(ma_convert, MutationAssessor_convert))

## Add FATHMM_converted_score column
ftm_convert <- sapply(1:nrow(ns),function(x) unlist(strsplit(ns$FATHMM_score[x], ";")))
ns$FATHMM_converted_score <- unlist(lapply(ftm_convert, FATHMM_convert))

ns <- ns %>% mutate(key=paste("chr",hg19_chr, ".", hg19_pos.1.based., ".", ref, ".", alt , sep=""))

## damaging based on functional impact for ns
damFunc <- cbind(
  key <- ns$key,
  sift <- data.frame(sapply(strsplit(ns$SIFT_pred, ";"), function(x) "D" %in% x)),
  pp2.div <- data.frame(sapply(strsplit(ns$Polyphen2_HDIV_score, ";"), function(x) length(which(as.numeric(x)>0.5))>0)),
  pp2.var <- data.frame(sapply(strsplit(ns$Polyphen2_HVAR_score, ";"), function(x) length(which(as.numeric(x)>0.5))>0)),
  lrt <- data.frame(sapply(strsplit(ns$LRT_pred, ";"), function(x) "D" %in% x)),
  mt <- data.frame(sapply(strsplit(ns$MutationTaster_converted_score, ";"), function(x) length(which(as.numeric(x)>0.5))>0)), ## "." checked for consistency length(which(as.numeric(".")>0.5))>0 == length(which(NA>0.5))>0 == FALSE
  ma  <- data.frame(sapply(strsplit(ns$MutationAssessor_converted_score, ";"), function(x) length(which(as.numeric(x)>0.65))>0)),
  ftm <- data.frame(sapply(strsplit(ns$FATHMM_converted_score, ";"), function(x) length(which(as.numeric(x)>=0.45))>0))
)
colnames(damFunc) <- c("key", "sift", "pp2.div", "pp2.var", "lrt", "mt", "ma", "ftm")

## I consider as damaging all the mutation that pass at least 5/7 function scores
#damFunc$damagingFunc <- apply(damFunc[,2:length(damFunc)],1, function(x) floor((sum(x)/length(x))*100) >= 70) ## We consider at least 5/7 (71%)
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



## Step 3: Infer damaging mutations for the splicing mutations

sc <- read.delim("sp.oac.dbnsfp.sc",stringsAsFactor=F,quote = "")

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

## Integrate both nonsynonymous and splicing predictions to the somatic mutation table and save it
damNS$key <- as.character(damNS$key)
damSC$key <- as.character(damSC$key)

for(i in 1:length(muts)){
  s=muts[[i]]
  s <- s %>% mutate(key=paste(chr, ".", end, ".",ref, ".", alt , sep="")) 
  s <- s %>% left_join(damNS, by=c("key")) %>% left_join(damSC, by=c("key"))
  ## Finally define damaging as TRUE or FALSE regardless if it is nonsynonymous or splicing
  s$damaging <- apply(s[,c("damagingNS", "damagingSC")],1, function(x) sum(x[!is.na(x)]) >= 1)
  ix <- s$ExonicFunc.refGene %in% c("stopgain","stoploss","stopgain","stoploss", "frameshift deletion", "frameshift insertion", "frameshift substitution")
  s[ix,]$damaging <- TRUE
  muts[[i]]=s
}

sapply(muts, function(x) sum(x$damaging))
save(muts, file="Rdata/dataset_mutations_damaging.Rdata")

# === ONCOCLUSTER =====

# Prepare INPUT files
# ucsc_refgene <- read.table("/Volumes/FC/Software/annovar_2015_12_14/humandb/hg19_refGene.txt", header = F)
# colnames(ucsc_refgene) <- c("bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames")
# 
# ## Replace the exon start with the cdsStart to exclude UTRs 5&3 and same for the cdsEnd
# pb <- txtProgressBar(min = 0, max = nrow(ucsc_refgene), style = 3)
# for (row in 1:nrow(ucsc_refgene)){
#   ucsc_refgene$exonStarts[row] <- paste(paste(c(ucsc_refgene$cdsStart[row],unlist(strsplit(ucsc_refgene$exonStarts[row], "\\,"))[which(unlist(strsplit(ucsc_refgene$exonStarts[row], "\\,")) > ucsc_refgene$cdsStart[row] & unlist(strsplit(ucsc_refgene$exonStarts[row], "\\,")) < ucsc_refgene$cdsEnd[row])]), collapse=","), ",", sep = "")
#   ucsc_refgene$exonEnds[row] <- paste(paste(c(unlist(strsplit(ucsc_refgene$exonEnds[row], "\\,"))[which(unlist(strsplit(ucsc_refgene$exonEnds[row], "\\,")) < ucsc_refgene$cdsEnd[row] & unlist(strsplit(ucsc_refgene$exonEnds[row], "\\,")) > ucsc_refgene$cdsStart[row])], ucsc_refgene$cdsEnd[row]), collapse=","), ",", sep = "")
#   setTxtProgressBar(pb, row)
# }
# close(pb)
# ucsc_refgene <- ucsc_refgene %>% select(name2, name, chrom, exonStarts, exonEnds)
# ucsc_refgene <- do.call(rbind, lapply_pb(split(ucsc_refgene, rownames(ucsc_refgene)), function(x) cbind(x[,1],x[,2],x[,3],  unlist(strsplit(x[,4],"\\,")), unlist(strsplit(x[,5],"\\,"))))) %>% data.frame()
# colnames(ucsc_refgene) <- c("Symbol", "Transcript.id", "chrom", "exon_start", "exon_end")
# ucsc_refgene$exon_start <- as.numeric(ucsc_refgene$exon_start)
# ucsc_refgene$exon_end <- as.numeric(ucsc_refgene$exon_end)
# ucsc_refgene = ddply(ucsc_refgene, .(Symbol, Transcript.id), mutate, n=1:length(exon_start), .progress = 'text')
# hg19_refgene = ucsc_refgene
# save(hg19_refgene, file="Rdata/hg19_refgene.Rdata")

n = as.list(names(muts))
muts=mapply(function(x,y){x$sample=y; return(x)}, muts, n, SIMPLIFY = F )

df_mut = as.data.frame(do.call("rbind", muts), stringsAsFactors=F)
df_mut = unrowname(df_mut)
nsyn = subset(df_mut, ExonicFunc.refGene=="nonsynonymous")
syn = subset(df_mut, ExonicFunc.refGene=="synonymous")

tmp = do.call('rbind', strsplit(sapply(strsplit(nsyn$AAChange.refGene,"\\,"), function(x) x[1] ),"\\:"))
tmp = data.frame(
  symbol_19014 = tmp[,1],
  symbol =  tmp[,1],
  Transcript.id =  tmp[,2],
  exon =  tmp[,3],
  nChange =  tmp[,4],
  aaChange =  tmp[,5],
  key = with(nsyn, paster(chr,'.',start,'.',end,'.',ref,'.',alt,'.', sample)),
  aa.position = gsub("[^0-9]", "", tmp[,5])
)
write.table(tmp, file="nsyn.onco", row.names = F, quote = F, sep = '\t')

tmp = do.call('rbind', strsplit(sapply(strsplit(syn$AAChange.refGene,"\\,"), function(x) x[1] ),"\\:"))
tmp = data.frame(
  symbol_19014 = tmp[,1],
  symbol =  tmp[,1],
  Transcript.id =  tmp[,2],
  exon =  tmp[,3],
  nChange =  tmp[,4],
  aaChange =  tmp[,5],
  key = with(syn, paster(chr,'.',start,'.',end,'.',ref,'.',alt,'.', sample)),
  aa.position = gsub("[^0-9]", "", tmp[,5])
)
write.table(tmp, file="syn.onco", col.names=T, row.names = F, quote = F, sep = '\t')

# cd /home/FC/Software/OncodriveClust_0.4.1
# source env/bin/activate
# module load python
# oncodriveclust -m 3 -c --cgc /home/mourikisa/data/CGC_phen.tsv ~/OAC/nsyn.onco ~/OAC/syn.onco /home/mourikisa/data/transcript_length.tsv -o ~/OAC/oncodriveclust-results.tsv

onco_out <- read.table("oncodriveclust-results.tsv", header = T, sep = "\t")
onco_out <- onco_out %>% subset(QVALUE<=0.1)
onco_in <- read.table("nsyn.onco", header = T, sep = "\t")
onco_in$Sample <- apply(onco_in, 1, function(x) unlist(strsplit(x[7], "\\."))[6])

## Each cluster a separate row - it will be easier to gather patients and check if clusters are costant across cancer types
onco_out <- onco_out %>% mutate(CLUST_COORDS=strsplit(CLUST_COORDS, "\\,\\[")) %>% 
  unnest(CLUST_COORDS) %>% mutate(CLUST_COORDS=gsub("\\[|\\]", "", CLUST_COORDS)) %>%
  separate(CLUST_COORDS, into=c("CLUST_COORDS_START", "CLUST_COORDS_END"), sep="\\,") %>% 
  separate(CLUST_COORDS_END, into=c("CLUST_COORDS_END", "NUMBER_OF_MUTS_IN_CLUST"), sep="\\:")

## Find which samples have mutations in the clusters
find_samples <- function(gene, start, end, onco_in){
  keys <- onco_in %>% subset(symbol_19014==gene & as.numeric(aa.position) >= start & as.numeric(aa.position) <= end) %>% select(key)
  keys <- paste(keys$key, collapse=",")
  return(keys)
}
onco_out$KEY <- apply(onco_out, 1, function(x) find_samples(x[1], as.numeric(x["CLUST_COORDS_START"]), as.numeric(x["CLUST_COORDS_END"]), onco_in))

## Expand the key components in the table
onco_out <- onco_out %>% mutate(KEY=strsplit(KEY, "\\,")) %>% unnest(KEY)
onco_out$SAMPLE <- apply(onco_out, 1, function(x) unlist(strsplit(x[14], "\\."))[6])

df_mut$key = with(df_mut, paster(chr,'.',start,'.',end,'.',ref,'.',alt,'.', sample))
df_mut$oncodriveClust = df_mut$key%in%onco_out$KEY

muts =  dlply( df_mut, .(sample) )

save(df_mut, muts, file="Rdata/dataset_mutations_damaging.Rdata")

# adding 19014

load("Rdata/19014_GeneSymbols.Rdata")
genes_19014 = unique(geneSymbols[,c("chromosome","tstart","tstop","Entrez","Symbol")])
i19 = with(genes_19014,GRanges(seqnames = chromosome, IRanges(as.numeric(tstart), as.numeric(tstop)) ))
im  = with(df_mut,GRanges(seqnames = chr, IRanges(start, end) ))

ov = findOverlaps( im, i19 )

df_mut$symbol_19014=NA
df_mut$entrez_19014=NA

df_mut[queryHits(ov), c("symbol_19014", "entrez_19014") ] = genes_19014[ subjectHits(ov),c("Symbol","Entrez")]

muts =  dlply( df_mut, .(sample) )

save(df_mut, muts, file="Rdata/dataset_mutations_damaging.Rdata")




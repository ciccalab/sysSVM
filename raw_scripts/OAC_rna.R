# RNA ======
#system("~/Lavoro/Software/bedtools2/bin/intersectBed -a genes.19014.bed -b extra/hg19_UCSCrefSeqGRCh37hg19.gtf -wa -wb >extra/hg19_UCSCrefSeqGRCh37hg19.19014.bed", wait = TRUE)
# system("~/Lavoro/Software/bedtools2/bin/intersectBed -a genes.19014.bed -b extra/RefSeq.gff -wa -wb > extra/RefSeq.19014.bed", wait = TRUE)

refseq  = read.delim("extra/RefSeq.19014.bed", h=F) # ucsc = read.delim("extra/hg19_UCSCrefSeqGRCh37hg19.19014.bed",h=F)

ref = (refseq[,ncol(refseq)])
ref = gsub("gene_id ", "", ref)
ref = gsub("transcript_id ", "", ref)
ref = as.data.frame(do.call(rbind, strsplit(ref, "; ")))
colnames(ref) = c("gene_id", "transcript_id")
ref$transcript_id = gsub(";", "", ref$transcript_id)
ref = cbind(refseq[,4:5], ref)
colnames(ref)[1:2] = c("entrez", "symbol")
ref = unique(ref)
ref$type=sapply(strsplit(ref$gene_id, "\\_"), function(x) x[1])

tlen = subset(refseq, V8=='exon')
tlen$len = with(tlen, V10-V9)
x = (tlen[,14])
x = gsub("gene_id ", "", x)
x = gsub("transcript_id ", "", x)
x = as.data.frame(do.call(rbind, strsplit(x, "; ")))
colnames(x) = c("gene_id", "transcript_id")
x$transcript_id = gsub(";", "", x$transcript_id)
tlen = cbind(tlen[,c(4:5,15)], x)
colnames(tlen)[1:2] = c("entrez", "symbol")
tlen = ddply(tlen, .(entrez, symbol, gene_id), summarise, len=sum(len), .progress = 'text')

t_max_clen = ddply(tlen, .(symbol), summarise, len=max(len), .progress = 'text')

save(ref, tlen, t_max_clen, file="Rdata/refseq_mapping_id.Rdata")


# USE DESEQ DATA which is similar to what thanos did for TCGA =====

load("Rdata/refseq_mapping_id.Rdata")

counts = read.delim("cufflinks/bygene_count.tab")
s.counts = t(apply(counts[,1:19], 2, quantile, probs=seq(0,1,.05)))
write.xlsx(s.counts, file = "Results/RNA_analysis.xlsx", "DESEQ Dist." )

counts$symbol_19014 = ref[match(rownames(counts), ref$gene_id),"symbol"]
s.counts.19014 =t(apply(counts[which(!is.na(counts$symbol_19014)),1:19], 2, quantile, probs=seq(0,1,.05)))
write.xlsx(s.counts.19014, file = "Results/RNA_analysis.xlsx", "DESEQ Dist.19014" , append=T)

counts$symbol_19014[which(is.na(counts$symbol_19014))] = rownames(counts)[which(is.na(counts$symbol_19014))]

m = melt(subset(counts, !is.na(symbol_19014)), id.vars = "symbol_19014")
m = ddply(m, .(variable, symbol_19014), summarise, value=sum(value), .progress = 'text')
l = dlply(m, .(variable))

g = unique(l[[1]]$symbol_19014)

for(i in 1:length(l)){
  names(l[[i]]) = unique(l[[i]]$variable)
  l[[i]]  = l[[i]][,3]
}

counts_bygene = as.data.frame(do.call('cbind',l))
rownames(counts_bygene) = g

s.counts.bygene = t(apply(counts_bygene[match(row.names(counts_bygene),unique(ref$symbol)),1:19], 2, quantile, probs=seq(0,1,.05), na.rm=T))
write.xlsx(s.counts.bygene, file = "Results/RNA_analysis.xlsx", "DESEQ Dist.19014 collapsed", append=T )


len <- t_max_clen[match(row.names(counts_bygene),t_max_clen$symbol), "len"]

# convert by experiment
fpkms <- apply(counts_bygene, 2, countToFpkm, len)
tpms = apply(counts_bygene, 2, fpkmToTpm)
save(tpms, counts_bygene, counts, file="Rdata/RNA_dataset.Rdata")
# t(apply(tpms, 2, quantile, probs=seq(0,1,.05)))

head(tpms)
t = as.data.frame(tpms)
t$symbol_19014 = rownames(t)
t = melt(t, measure.vars = colnames(tpms), id.vars = "symbol_19014")
x = read.delim("cufflinks/rna_seq_mapping.txt", h = F, stringsAsFactors = F)
x[, 3] = gsub("\\-","\\.",x[,3])
t$sample = x[ match(as.character(t$variable), x[,3]), 2 ]

l = dlply(t, .(sample)) 
x = read.delim('rna_samples.txt', h=T)
df_rna = as.data.frame(do.call(rbind, l[ names(l)!="LP6007400-DNA_A01" ])) # LP6007400-DNA_A01 normal

mNorm = l[['LP6007400-DNA_A01']]$value
names(mNorm) =   l[['LP6007400-DNA_A01']]$symbol_19014

colnames(df_rna)[3] = "TPM"
df_rna$meanNormal = mNorm[df_rna$symbol_19014]
df_rna = ddply(df_rna, .(sample),  mutate, q2=quantile(TPM, seq(0,1,.25))["25%"], q10=quantile(TPM, seq(0,1,.25))["75%"], .progress="text")

df_rna$expr_class=NA                                                
df_rna$expr_class[ which( df_rna$TPM > df_rna$q10 ) ]                      = "HE"
df_rna$expr_class[ which( df_rna$TPM > df_rna$q2 & df_rna$TPM <= df_rna$q10 ) ]  = "ME"
df_rna$expr_class[ which( df_rna$TPM>0.1 & df_rna$TPM<=df_rna$q2 ) ]          = "LE"
df_rna$expr_class[ which( df_rna$TPM <= 0.1 & df_rna$meanNormal > 0.1 ) ]  = "not expr"
df_rna$expr_class[ which( df_rna$TPM <= 0.1 & df_rna$meanNormal <= 0.1 ) ] = "not expr NT"

df_rna$logFC=with(df_rna, log2( (TPM+.00000001)/(meanNormal+.00000001 )))

save(df_rna, tpms, counts_bygene, counts, file="Rdata/dataset_rna.Rdata")

write.xlsx(ddply(df_rna, .(sample), function(x) summary(x$TPM)), file = "Results/RNA_analysis.xlsx", "DESEQ TPM 19014 collapsed" , append=T)



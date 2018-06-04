source("config.R")

#***********************************
#           CNVS
#***********************************

a = list.files(path="ascat/")
a = as.list(a)
b = lapply(strsplit(list.files(path="ascat/"),"_vs_"), function(x) x[1])
  
l = vector("list",length(a))
for(i in 1:length(a)){
  l[[i]] =  read.csv( paster("ascat/",a[[i]],"/1.3/results/",b[[i]],".copynumber.caveman.csv"))
  l[[i]]$Chromosome = paster("chr",l[[i]]$Chromosome)
  l[[i]]$Major_CN = with(l[[i]], Total_CN-Minor_CN)
  l[[i]] = l[[i]][,c('Chromosome','Start','End', 'Total_CN', 'Major_CN', 'Minor_CN' )]
  l[[i]]$CNV_type = NA
  l[[i]]$CNV_type[ which(l[[i]]$Total_CN>2) ] =  "Gain"
  l[[i]]$CNV_type[ which(l[[i]]$Total_CN<2) ] =  "Loss"
  l[[i]]$CNV_type[ which(l[[i]]$Total_CN==2 & l[[i]]$Minor_CN==0) ] =  "Loss"
 
  l[[i]]$real_CNV_type = NA
  l[[i]]$real_CNV_type[ which(l[[i]]$Total_CN>2) ] =  "Gain"
  l[[i]]$real_CNV_type[ which(l[[i]]$Total_CN<2) ] =  "Loss"
  l[[i]]$real_CNV_type[ which(l[[i]]$Total_CN>=2 & l[[i]]$Minor_CN==0) ] =  "LOH"

}

names(l) = unlist(b)
l = mapply( function(x,y){ x$sample=y; return(x)}, l, b, SIMPLIFY = F)
df_cnv = as.data.frame(do.call("rbind",l))

chr = read.delim('/Volumes/FC/DB/hg19/hg19.chrom.sizes', h=F)
rownames(chr)=chr[,1]

df_cnv$len = with(df_cnv, End-Start)
df_cnv$perc = df_cnv$len/chr[df_cnv$Chromosome,2]
df_cnv$CNV_cut = df_cnv$Total_CN
df_cnv$CNV_cut[which(df_cnv$Total_CN>=5)] = 5

df_cnv = subset(df_cnv, Total_CN!=2 )

st = ddply(df_cnv, .(sample,CNV_cut), summarise, p = sum(len)/3095677412, .progress = 'text')
st$CNV_cut = factor(st$CNV_cut, levels=c("0","1","3","4","5"))
pdf(file="Results/CNV_percentage_of_genome.pdf", h=6,w=8)
ggplot(st, aes(x=sample,y=p, fill=CNV_cut))+geom_bar(stat='identity', color="black")+theme_boss()+scale_fill_manual(values=c("darkblue","lightblue","pink","orange","tomato"))+xlab("")+
  scale_y_continuous(limits=c(0,1))+ylab("% of Genome")+coord_flip()
dev.off()

st = ddply(st, .(sample), summarise, p = sum(p), .progress = 'text')
write.xlsx(st, "Results/Variant_analysis.xlsx","CNV_genome_coverage", append=T, row.names=F)

save(df_cnv,l, file="Rdata/dataset_cnvs.Rdata")

annotate_genes = function(x){
  write.bed(x[,c("Chromosome","Start","End","sample")], file="tmp.bed")
  system("~/Lavoro/Software/bedtools2/bin/intersectBed -a genes.19014.bed -b tmp.bed -f 0.25 -wa -wb > tmp.genes.bed", wait = TRUE)
  tmp = read.table("tmp.genes.bed", h=F)
  tmp$key=paster(tmp$V6,".",tmp$V7,'.',tmp$V8,'.',tmp$V9)
  key=paster(x$Chromosome,".",x$Start,'.',x$End,'.',x$sample)
  y = cbind(tmp[,1:5], x[match(tmp$key, key),])
  colnames(y)[1:5] = c("chrom","start","end","entrez_19014","symbol_19014")
  system("rm tmp.bed")
  system("rm tmp.genes.bed")
  y
}

df_cnv_19014 = annotate_genes(df_cnv)

save(df_cnv,l, df_cnv_19014, file="Rdata/dataset_cnvs.Rdata")

x = ddply(df_cnv_19014,.(symbol_19014, chrom,start, CNV_type), summarise, n=length(symbol_19014), .progress = 'text')
x$chrom = factor(x$chrom, levels=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'))

pdf(file="Results/area_gain.pdf",w=15,h=6)
p1=ggplot(subset(x,CNV_type=="Gain"), aes(x=start,y=n))+geom_line()+geom_area(fill="tomato",alpha=.5)+theme_boss()+facet_wrap(~chrom, ncol=24)+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+xlab('')
p2=ggplot(subset(x,CNV_type=="Loss"), aes(x=start,y=-n))+geom_line()+geom_area(fill="lightblue",alpha=.5)+theme_boss()+facet_wrap(~chrom, ncol=24)+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+xlab('')
grid.arrange(p1,p2, nrow=2)
dev.off()


# === subsetting to Thanos' cutoffs =======

df_cnv = subset(df_cnv, Total_CN>=4 | Total_CN<2)

st = ddply(df_cnv, .(sample,CNV_cut), summarise, p = sum(len)/3095677412, .progress = 'text')
st = ddply(st, .(sample), summarise, p = sum(p), .progress = 'text')
write.xlsx(st, "Results/Variant_analysis.xlsx","CNV_genome_coverage_CNV>=4_CNV<2", append=T, row.names=F)

# st = ddply(df_cnv, .(sample,Chromosome), summarise, gain=sum(CNV_type=='Gain'), p.gain = sum(perc[CNV_type=='Gain']),loss=sum(CNV_type=='Loss'), p.loss = sum(perc[CNV_type=='Loss']), .progress = 'text')
# m = melt(st, id.vars = c("sample", "Chromosome"), measure.vars = c("p.gain", 'p.loss'))
# m$Chromosome = factor(m$Chromosome, levels=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'))
# pdf(file="Results/CNV_regions_per_samples.pdf", h=20,w=15)
# ggplot(m, aes(x=Chromosome,y=value, fill=variable))+geom_bar(stat='identity', position=position_dodge())+theme_boss()+facet_wrap(~sample, nrow=18)
# dev.off()

# === subsetting to Recurrent GISTIC Regions =======

# a = c("chr13-40100001-45200000","chr20-58400001-63025520",
# "chr1-1-2300000","chr2-102700001-106000000","chr3-192300001-198022430","chr9-137400001-141213431","chr8-6200001-12700000","chr1-142600001-147000000","chr1-155000001-156500000")
# 
# d = c("chr4-187100001-191154276","chr8-2200001-6200000","chr21-14300001-16400000")

load("Rdata/GISTIC.Rdata")

a = as.data.frame(do.call("rbind", strsplit(a[,1], "[[:punct:]]")))
d = as.data.frame(do.call("rbind", strsplit(d[,1], "[[:punct:]]")))

load("Rdata/dataset_cnvs.Rdata")

df_cnv$gistic=FALSE
df_cnv_19014$gistic=FALSE

ia = GRanges(seqnames = a[,1], IRanges(as.numeric(a[,2]),as.numeric(a[,3])))

ix = GRanges(seqnames = df_cnv[,1], IRanges(as.numeric(df_cnv[,2]),as.numeric(df_cnv[,3])))
ov = findOverlaps(ix,ia)
df_cnv$gistic[ queryHits(ov) ] = TRUE
df_cnv$gistic[ df_cnv$CNV_type=="Loss" ] = FALSE

ix = GRanges(seqnames = df_cnv_19014[,1], IRanges(as.numeric(df_cnv_19014[,2]),as.numeric(df_cnv_19014[,3])))
ov = findOverlaps(ix,ia)
df_cnv_19014$gistic[ queryHits(ov) ] = TRUE
df_cnv_19014$gistic[ df_cnv_19014$CNV_type=="Loss" ] = FALSE

id = GRanges(seqnames = d[,1], IRanges(as.numeric(d[,2]),as.numeric(d[,3])))

ix = GRanges(seqnames = df_cnv[,1], IRanges(as.numeric(df_cnv[,2]),as.numeric(df_cnv[,3])))
ov = findOverlaps(ix,id)
tmp = which(df_cnv$CNV_type=="Loss")
df_cnv$gistic[ intersect(queryHits(ov), tmp) ] = TRUE

ix = GRanges(seqnames = df_cnv_19014[,1], IRanges(as.numeric(df_cnv_19014[,2]),as.numeric(df_cnv_19014[,3])))
ov = findOverlaps(ix,id)
tmp = which(df_cnv_19014$CNV_type=="Loss")
df_cnv_19014$gistic[ intersect(queryHits(ov), tmp) ] = TRUE

save(df_cnv,l, df_cnv_19014, file="Rdata/dataset_cnvs.Rdata")


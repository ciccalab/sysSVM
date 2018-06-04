SAMPLEID ="LP2000105-DNA_A01_vs_LP2000102-DNA_A01"
SNPS     ="OAC/input/results/LP2000102-DNA_A01.germlineHETpos.bed"
FNORM    ="OAC/input/results/LP2000102-DNA_A01.readcounts.txt"
FTUMOR   ="OAC/input/results/LP2000105-DNA_A01.readcounts.txt"
FOUT     ="OAC/LP2000105-DNA_A01_vs_LP2000102-DNA_A01.CHAT.txt"
DIROUT   ="OAC/output/"

library("CHAT")

get_seg_mat=function(SAMPLEID, SNPS, FNORM, FTUMOR
                     # , THRESHOLD_COVERAGE=10
                     ){
  s = read.delim(SNPS, h=F)
  n = read.delim(FNORM)
  t = read.delim(FTUMOR)

  id   = paste0(s[,1],'.', s[,3])
  id.n = paste0(n[,1],'.', n[,2])
  id.t = paste0(t[,1],'.', t[,2])

  x = which(id.n%in%id); if(length(x)>0) n = n[x, ] else stop(message("ERROR: no snps in normal"))
  x = which(id.t%in%id); if(length(x)>0) t = t[x, ] else stop(message("ERROR: no snps in tumor"))

  # n = subset(n, Reference.base%in%c(''))

  get_variants = function(x){
    ref = x[3]
    cols = c(   "A"="A.count"
                , "C"="C.count"
                , "G"="G.count"
                , "T"="T.count")

    cols = cols[names(cols)!=ref]

    if(sum(x[cols]==0)){
      return(c(x[c('Chromosome','Position','Depth')],'Depth.variant'=0))
    }else{

      return(c(x[c('Chromosome','Position','Depth')],'Depth.variant'=max( x[cols] )))
    }
  }

  p.n = t( apply(n,1,get_variants) )
  p.n = as.data.frame(p.n, stringsAsFactors = F)
  p.n$Depth=as.numeric(p.n$Depth)
  p.n$Depth.variant=as.numeric(p.n$Depth.variant)

  p.t = t( apply(t,1,get_variants) )
  p.t = as.data.frame(p.t, stringsAsFactors = F)
  p.t$Depth=as.numeric(p.t$Depth)
  p.t$Depth.variant=as.numeric(p.t$Depth.variant)

  rm(t,n)

  p.n$Position = as.numeric(p.n$Position)
  p.t$Position = as.numeric(p.t$Position)

  p.n = p.n[order(p.n$Chromosome, p.n$Position, decreasing = F),]
  p.t = p.t[order(p.t$Chromosome, p.t$Position, decreasing = F),]

  id.n = paste0(p.n[,1],'.', p.n[,2])
  id.t = paste0(p.t[,1],'.', p.t[,2])

  p.t = p.t[match(id.n, id.t),]

  # From https://sourceforge.net/p/clonalhetanalysistool/wiki/Extract%20LRR%20and%20BAF%20signals%20from%20next%20generation%20sequencing%20data/
  # vv.germ=which(a0>0&b0>0)
  # vv.cov=which(s>=thr.cov&s0>=thr.cov)
  # vv=intersect(vv.germ,vv.cov)
  #

  # alt = which( p.n$Depth.variant>0 & p.t$Depth.variant>0 )
  # tot = which( p.n$Depth>THRESHOLD_COVERAGE & p.n$Depth>THRESHOLD_COVERAGE )
  #
  # select = intersect(alt,tot)
  # p.n = p.n[select,]
  # p.t = p.t[select,]

  seg.mat=cbind( SAMPLEID
                 ,gsub("chr","",p.n[,1])
                 ,p.n[,2]
                 ,log2(p.t$Depth/p.n$Depth)
                 ,p.t$Depth.variant/p.t$Depth
                 ,p.n$Depth.variant/p.n$Depth)
  colnames(seg.mat)=c('sampleID','chr','pos','LRR','BAF','BAF-n')

  return(seg.mat)
}

getSegChr.Seq <- function(seg.mat,sampleid='Sample'
                          ,bin=1000 #number of markers contained in each bin
                          ,cbs=FALSE
                          ,thr.hets=0.15 #lower threshold of calling homozygous markers. BAF<=thr.hets or BAF>=1-thr.hets are considered homozygous.
                          ){
  dd.dat=c()
  id=seg.mat[1,1]
  for(cc in 1:22){
    vv=which(seg.mat[,2]==cc)
    bb.chr=seg.mat[vv,c(2,2,3,5,6)]
    ll.chr=cbind(seg.mat[vv,c(2,2,3,4)],rep(0,length(vv)))
    mode(bb.chr)=mode(ll.chr)='numeric'
    colnames(bb.chr)[4:5]=colnames(ll.chr)[4:5]=c(sampleid,paste(sampleid,'-normal',sep=''))
    if(cbs){
      dat.chr=getSegChr.CBS(bb.chr,ll.chr,sam.col=4, thr.hets=thr.hets,data.type='log', bin=bin)
      } else {
      dat.chr=getSegChr(bb.chr,ll.chr,sam.col=4, thr.hets=thr.hets,data.type='log', bin=bin)
      }
    dd.dat=rbind(dd.dat,dat.chr)
  }
  rownames(dd.dat)=dd.dat[,1]
  dd.dat = dd.dat[,2:8]
  mode(dd.dat)='numeric'

  return(dd.dat)
}

# RUN =============

seg.mat = get_seg_mat(SAMPLEID, SNPS, FNORM, FTUMOR, THRESHOLD_COVERAGE=30)

dd.dat  = getSegChr.Seq(na.omit(seg.mat), SAMPLEID, bin=5e4, thr.hets=.2 )

save(dd.dat,file=paste0(DIROUT,SAMPLEID,'.seg.Rdata'))


para=getPara()
para$datafile = paste0(DIROUT,SAMPLEID,'.seg.Rdata')
para$savefile = paste0(DIROUT,SAMPLEID,'.AGP.txt')
para$pngdir = DIROUT
# para$BAFfilter=100 # DNA segments with BAF markers below this value are removed to reduce noise.
# para$std.BAF=0.1
# para$std.LRR=0.1
para$is.normalize=FALSE
para$thr.penalty=1000

## AGP estimation
getAGP(para=para)


para.s=getPara.sAGP()
para.s$inputdata  = paste0(DIROUT,SAMPLEID,'.seg.Rdata')
para.s$purityfile = paste0(DIROUT,SAMPLEID,'.AGP.txt')
para.s$savedata   = paste0(DIROUT,SAMPLEID,'.sAGP.Rdata')


## sAGP estimation
getsAGP(para=para.s)


## SegPurity
oo       = getOrigin(dd.dat,para=para)
sAGP.dat = getSegPurity(dd.dat,oo, AGP=1, para=para.s)
write.csv(sAGP.dat, file=paste0(DIROUT,SAMPLEID,".csv"), row.names = F)

cnvs = read.delim("OAC/LP2000105-DNA_A01_CNVs.tsv")

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
genes <- getBM(attributes = c( 'hgnc_symbol','chromosome_name','start_position','end_position'),
                   filter =  c("entrezgene"),
                   values = cnvs$entrez,
                   mart = ensembl)
genes = subset(genes, chromosome_name%in%as.character(1:22))

# cnv_filename       = "Example/OAC/cnv_file_259_OAC.simple.tsv"
# cnvs = read.cnv.file(cnv_filename)
# cnvs=subset(cnvs, sample==SAMPLEID)

library(GenomicRanges)
x = GRanges(seqnames = sAGP.dat[,1], IRanges(start=sAGP.dat[,2], end=sAGP.dat[,3]) )
y = GRanges(seqnames = genes[,'chromosome_name'], IRanges(start=genes[,'start_position'], end=genes[,'end_position']) )


o =findOverlaps(y,x, minoverlap = 10)

clon=cbind(genes[queryHits(o),],new.dd[subjectHits(o),])
# clon = cbind(clon, cnvs[match(clon$hgnc_symbol,cnvs[,'symbol']),])

write.csv(clon, file=paste0(DIROUT,SAMPLEID,".purity.csv"), row.names = F)


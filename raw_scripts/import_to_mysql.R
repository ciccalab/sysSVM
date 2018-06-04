source("~/Lavoro/pancancer_model_training/config.R")
load("Rdata/EAC.15_02_16.euristics.Rdata")

cols=c('Cancer_type', 'Sample','Entrez','no_NSI_muts', 'no_TRUNC_muts', 'no_NTDam_muts', 'no_NTDamFunction_muts' ,
'no_NTDamCons_muts','no_NTDamSC_muts','no_GOF_muts'   ,'Copy_number' ,'TPM'    ,'vogel'  ,'Genic'  ,'inCDD'  ,      'alldomains'  ,  
'private_domains','commondomains'   ,'memberofcomplex' ,'High'      ,'Low'      ,'Medium' , 'NotExpressed' , 'Length_fullrefseq',    
'WGD'  ,      'degree'  ,   'betweenness','hub' ,     'central' , 'CNVGain',  'CNVLoss',  'ExpT_ME',  'ExpT_HE' ,  
'ExpT_LE' , 'ExpT_NE'  ,'ExpT_NET','old'   ,'young'  ,'luca'     ,'eukaryotes', 'metazoans'  ,'vertebrates' ,'opisthokonts',       
'mammals' ,'selective' ,'always_expressed','middle'  , 'one_tissue'   ,  'never_expressed' ,  'tot_tissues'  )      

eac.table$Cancer_type='OAC'
colnames(eac.table)[1] = 'Sample'
colnames(eac.table)[3] = 'Entrez'
colnames(eac.table) = gsub("\\.","_", colnames(eac.table))
cols[which(cols%nin%colnames(eac.table))]

df = eac.table[,cols]
df$TPM[is.na(df$TPM)] = -1

write.table(df, file="/Volumes/ceredam/novel_driver_prediction/OAC.tsv", col.names=F, row.names=F, quote=F, sep="\t")


# load data local infile '/home/ceredam/novel_driver_prediction/OAC.tsv' into table OAC;


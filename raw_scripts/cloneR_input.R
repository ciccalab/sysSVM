## Script to prepare the input for CloneR

##Excluded samples
excl = c("LP6005690-DNA_E02_vs_LP6005689-DNA_E02",
"LP6008280-DNA_F02_vs_LP6008264-DNA_F02",
"LP6008202-DNA_F01_vs_LP6008201-DNA_F01",
"LP6005935-DNA_C01_vs_LP6005934-DNA_C01",
"LP6008031-DNA_E03_vs_LP6008032-DNA_A04")

## Load histopathology data
sample = read.table("~/athena/data/OAC/Clinical_data/release_24/sample.tsv", header = T, sep="\t")
sex = read.table("~/athena/data/OAC/Clinical_data/release_24/donor.tsv", header = T, sep="\t")
## Get also missing genders (I dont have all the clinical info)
sex2 = read.table("~/athena/data/OAC/Clinical_data/release_24/missing_genders.tsv", header = F, sep="\t") %>% 
    rename(sample=V1, gender=V2) %>% 
    separate(sample, into=c("tumour", "normal"), sep="_vs_", remove = F) %>% 
    select(tumour, gender)

sex = sample %>% select(icgc_donor_id, submitted_sample_id) %>% left_join(sex%>%select(icgc_donor_id, donor_sex)) %>% 
    select(submitted_sample_id, donor_sex) %>% rename(tumour=submitted_sample_id) %>% 
    mutate(gender=ifelse(donor_sex=="female", "F", "M")) %>% 
    select(tumour, gender)
sex = rbind(sex, sex2)

## Get the mutations
load("~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/clonality/muts_129_66_71_OACs_annovar_dbnsfp_oncodriveClust.Rdata")


samples = muts %>% select(sample) %>% unique %>% subset(!sample%in%excl) %>% 
    separate(sample, into=c("tumour", "normal"), sep="_vs_", remove = F) %>% 
    left_join(sex) %>% mutate(gender2=ifelse(gender=="M", "T", "F"))

## Read in the tumour content as well
load("~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/histopathology_tumour_content/tumour_content_261_OACs.Rdata")
dataset2 <- read_excel("~/Desktop/Copy of Updates on OAC analysis.xlsx", skip = 1)
colnames(dataset2) = c("sample", "X1", "whole_slide", "marked_area")
dataset2 = dataset2 %>% mutate(Cellularity=paste0(whole_slide, "/", marked_area)) %>% subset(!(is.na(whole_slide) & is.na(marked_area)))
dataset <- read_excel("~/Desktop/Cellularities.xlsx")
tumour_content = dataset %>% select(sample, Cellularity) %>% full_join(dataset2%>%select(sample, Cellularity)) %>% rename(tumour=sample)
samples = samples %>% left_join(tumour_content)
## Here I just exported and assigned manual the correct value cause in some cases there is plain text
samples = samples%>%select(sample, gender2, Cellularity)
colnames(samples) = c("sample", "gender", "tumour_content")
write.table(samples, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/clonality/sample_file_259_OAC.tsv", sep="\t", row.names = F, quote = F)



## Mutation file
## Get only the exonic mutations
muts = muts %>% subset(Func.refGene!="" &
                               !grepl("downstream", muts$Func.refGene) &
                               !grepl("upstream", muts$Func.refGene) &
                               !grepl("intergenic", muts$Func.refGene) &
                               !grepl("ncRNA", muts$Func.refGene) &
                               !grepl("intronic", muts$Func.refGene) &    
                               !grepl("UTR", muts$Func.refGene))


# Select Columns: sample, chr, position , reference, variant , allele frequency, gene symbol, mutation type (i.e. SNV, InDel)
muts_df = muts %>% select(sample, chr, start, end, ref, alt, freq, entrez_19014, symbol_19014) %>% subset(!(is.na(entrez_19014) & is.na(symbol_19014)))
## Add mutation type
muts_df = muts_df %>% mutate(mutation_type=ifelse(nchar(ref)>1 | nchar(alt)>1, "InDel", "SNV"))
## Export only the samples for which we have tumour content
muts_df = muts_df %>% subset(sample%in%samples$sample)
write.table(muts_df, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/clonality/mutation_file_259_OAC.tsv", sep="\t", row.names = F, quote = F)

## CNV file
load("~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/clonality/cnvs_129_66_71_OACs.Rdata")
## Following 2 lines run on athena - they are here for reference (that's why they are commented out)
#cnvsPlusRaw = getRawCN(cnvs) ## function can be found in a separate file under /athena/data/OAC
#save(cnvsPlusRaw, file="/home/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/Rdata/cnvPlusRaw_129_66_71_OACs.Rdata")

cnvs_df = cnvsPlusRaw %>% select(sample, Chromosome, Start, End, CNV_type, Major_CN, Minor_CN, Raw_CN) %>% unique
cnvs_df = cnvs_df %>% subset(sample%in%samples$sample)
write.table(cnvs_df, file="~/athena/data/OAC/Combined_ICGC_cohort_129_66_71/clonality/cnv_file_259_OAC.tsv", sep="\t", row.names = F, quote = F)



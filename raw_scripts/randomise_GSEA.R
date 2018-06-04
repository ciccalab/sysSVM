rm(list=ls())

## Define command line arguments
args = commandArgs(trailingOnly = TRUE)

## Parameters
N = as.numeric(args[1])
TIMES = as.numeric(args[2])
RES_DIR = args[3]
CHUNK = args[4]

if (!dir.exists(RES_DIR)){
    dir.create(RES_DIR)
}
    
geneInfo_fn="~/data/geneInfoNCG5.Rdata"
gene_sets_dir="~/data/geneSets/reactome_12_12_16/"
gene_sets_gene_start=4
gene_sets_identifiers = "reactome"
sub_identifier="symbol"
scores_fn="~/data/OAC/129_OAC/ML/OAC/syscans_4_with_gains_in_sys_cans.Rdata"
cohort_fn="~/data/OAC/129_OAC/Rdata/mainCohort.Rdata"
remove.amplifications=T
remove.cgc=T
    
## Library cleaning
require(plyr)
require(tidyr)
require(dplyr)
    
if(is.null(N)){
    stop("No number of genes")
}
    
## pathway parameter in the following function can be kegg/reactome/all
## sub parameter can be entrez/symbol
getPathway = function(pathway="kegg", sub=sub_identifier, genes_start=gene_sets_gene_start,
                        path_dir=gene_sets_dir,
                        gi_fn=geneInfo_fn){
        
        make_pathway_list = function(pathMsigDbFile) {
            inputFile <- pathMsigDbFile
            con  <- file(inputFile, open = "r")
            c = 1
            pathway.list <- vector(mode="list",length=0)
            message("Loading gene sets....")
            while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
                myVector <- do.call("rbind",strsplit(oneLine, "\t"))
                t = vector(mode="list",length=1)
                t[[1]] = myVector[genes_start:length(myVector)]
                names(t) = myVector[1]
                pathway.list = c(pathway.list,t)
                c = c+1
            } 
            
            close(con)
            return(pathway.list)
        }
        
        fns = list.files(path_dir)
        pattern=paste0(pathway, ".*", sub)
        fn = fns[grepl(pattern, fns)]
        if(length(fn)>1 | length(fn)==0){
            stop("No or multiple files found. Please check the parameters provided")
        }
        
        message(paste0("Pathway file: ", fn))
        p = make_pathway_list(paste0(path_dir, fn))
        message(paste0("Number of pathways: ", length(p)))
        if(sub=="entrez"){
            genes_p = as.numeric(unique(unlist(p)))
        }else if(sub=="symbol"){
            genes_p = unique(unlist(p))
        }
        message(paste0("Number of genes: ", length(genes_p)))
        ## Get gene info to check overlap with 19014
        load(gi_fn)
        geneInfo = geneInfo %>% subset(duplicability!=-1) %>% select(symbol, entrez)
        if(sub=="entrez"){
            message(paste0("Number of genes (19,014): ", length(genes_p[genes_p%in%geneInfo$entrez])), "\n")
        }else if(sub=="symbol"){
            message(paste0("Number of genes (19,014): ", length(genes_p[genes_p%in%geneInfo$symbol])), "\n")
        }
        return(p)
        
    }
    
## Get the genes to sample from
load(scores_fn)
load(geneInfo_fn)
## Make sure you update the symbols to work with gene sets with symbols
syscan = syscan %>% select(-symbol) %>% left_join(geneInfo%>%select(entrez, symbol))
## Exclude "CGC discarded"
syscan = syscan %>% subset(gene_type!="CGC_discarded")
    
## ------ Add mutational signature information to the scores ---------
load(cohort_fn)
mainCohort = mainCohort %>% rename(sample=Primary_tumour_solid_tissue)
syscan = syscan %>% left_join(mainCohort%>%select(sample, Group)%>%unique)
    
## Filters
if(remove.amplifications){
    syscan = syscan %>% subset(CNVGain!=1)
}
if(remove.cgc){
    cgcs = syscan %>% subset(!is.na(Group)) %>% subset(gene_type=="cgc") %>% select(symbol) %>% unique %>% .$symbol
    syscan = syscan %>% subset(gene_type!="cgc") 
}
## Exclude samples that do not correspond to 
syscan = syscan %>% subset(!is.na(Group))
    
## get genes for sampling
genes_to_sample = syscan %>% select(symbol) %>% unique %>% .$symbol
## get gene sets
gene_sets = getPathway(pathway = gene_sets_identifiers)
    
    ## x is the gene set and y your gene list
run.hypergeometric = function(x, y, universe=19014){
    P = phyper(length(intersect(x,y)), length(x), universe-length(x), length(y), lower.tail = F)
    p = data.frame(p.value=P, gene.set.length=length(x), genes.in.gene.set=length(intersect(x,y)))
    return(p)
}

print("Creating random samples of genes...")
   
ex = replicate(TIMES, sample(genes_to_sample, N))
ex = split(ex, rep(1:ncol(ex), each = nrow(ex)))
## Add CGCs in the random samples
ex_sampled = ex
ex = lapply(ex, function(x) x=c(x, cgcs))


print("Running GSEA for random samples...")
rgsea = lapply(ex, function(x){
        
    l = lapply(gene_sets, function(y) run.hypergeometric(y, x))
    l = lapply(seq_along(l), function(z, n, i){
        ldply(z[[i]]) %>% mutate(pathway=n[[i]])
    }, z=l, n=names(l))
    l = do.call(rbind, l) %>% rename(value=V1, type=.id)
    l = l %>% unique %>% spread(type, value)
        
    l$fdr = p.adjust(l$p.value, method="BH")
    l$number.of.genes.tested = length(x)
    return(l)
})
    
## Make the list of gene set enrichments a table
rgsea = lapply(seq_along(rgsea), function(y,n,i){
    y[[i]] %>% mutate(iter=n[[i]])
}, y=rgsea, n=names(rgsea))    
    
rgsea = do.call(rbind, rgsea)
    
res = list(gsea=rgsea, gene_samples=ex_sampled, length_cgcs=length(cgcs), length_sys_cans=length(genes_to_sample))

## Save results
save(res, file=paste0(RES_DIR, "/random_sampling_GSEA_", CHUNK, ".Rdata"))
print("END")



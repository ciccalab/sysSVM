## You source this file as soon as you start working on this project to load all the
## relevant tables and functions

pkgsR  =  c('plyr',
            'tidyr',
            'caret',
            'ggplot2',
            'grid',
            'gridExtra',
            'tibble',
            'reshape2',
            'doMC',
            'xlsx',
            'readxl',
            'scales',
            'doParallel',
            'foreach',
            'RMySQL',
            'pROC',
            'e1071',
            'scales',
            'RColorBrewer',
            'sampling',
            'data.table',
            'dplyr')
for (pkgR in pkgsR){
    if (!pkgR %in% rownames(installed.packages())) {
        install.packages(pkgR, dependencies=TRUE, INSTALL_opts = c('--no-lock'), repos = "https://cran.ma.imperial.ac.uk/")
        library(pkgR, character.only = TRUE, quietly = TRUE)
    } else {
        library(pkgR, character.only = TRUE, quietly = TRUE)
    }
}

source("http://bioconductor.org/biocLite.R")

pkgs  =  c("RDRToolbox")
for (pkg in pkgs){
    if (!pkg %in% rownames(installed.packages())) {
        biocLite(pkg)
        library(pkg, character.only = TRUE, quietly = TRUE)
    } else {
        library(pkg, character.only = TRUE, quietly = TRUE)
    }
}

##update.packages(repos=biocinstallRepos(), ask=FALSE)

rm(pkg,pkgR,pkgs,pkgsR)
## Load dplyr always last
detach("package:dplyr", unload=TRUE)
library(dplyr)

theme_boss <- function(base_size = 12, base_family = "sans"){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border=element_blank()
    )
}

## Needed for some functions
registerDoMC(3)
cancer_types <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM",
                  "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD",
                  "LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM",
                  "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

drivers_secrier = c("SMYD3", "RUNX1", "CTNNA3", "RBFOX1", "AGBL4", "CDKN2A",
                    "CDKN2B", "SAMD5", "CDK14", "KIF26B", "THADA", "SASH1",
                    "MECOM", "JUP", "IKZF3", "FHIT", "WWOX", "MACROD2", "IMMP2L",
                    "CCSER1", "PDE4D", "NAALADL2", "PARK2", "PARD3B", "PRKG1",
                    "TP53", "SMAD4", "ARID1A", "KCNQ3", "CCDC102B", "CYP7B1")


## Function to get KEGG pathways using REST API
getKEGG = function(r=NULL, fn="Desktop/kegg_test.txt", app=T){
    require(KEGGREST)
    path_ids = unique(unlist(keggLink("pathway", "hsa")))
    cat(paste0("Pathways found: ", length(path_ids)), "\n")

    ## I placed this here in case there is a limitation in the number of times you can query KEGG
    if(is.null(r)){
        r = 1:length(path_ids)
    }

    for(i in r){
        kegg_path = keggGet(path_ids[i]) ## This returns an object with name, description etc
        kegg_path = keggGet(path_ids[i]) ## This returns an object with name, description etc
        nm = kegg_path[[1]]$NAME
        cat(paste0(i, ") ", nm), "\n")
        if("GENE"%in%names(kegg_path[[1]])){
            genes = unlist(lapply(strsplit(kegg_path[[1]]$GENE[seq(2, length(kegg_path[[1]]$GENE), 2)], ";"), function(x) x[1]))
            cat(paste(nm, paste(genes, collapse="\t"), sep = "\t"), file = fn, append = app, "\n")
            #paths[[nm]] = genes
        }else{
            cat("skipped", "\n")
            next
        }
    }
}


getGeneCoordsHg19 = function(){
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    ## Get gene coordinates
    g = genes(txdb)
    return(g)
}

## The wat I use this function is the following:
## I give the list I want and it returns a list with the number of common elements for all combinations of the elements of the list
## Prepare the gene sets for the intersection (I did the intersection only on the mutated genes initially)
# paths = gsep %>% select(pathway) %>% unique %>% .$pathway
# gene_sets = list()
# for(p in paths){
#     genes = gsep %>% subset(pathway==p) %>% .$symbol %>% unique
#     gene_sets[[p]] = genes
# }
getListIntersection = function(l, separator=".", exclude=NULL){

    ## excldue argument helps to remove elements from the lists during each comparison
    ## not implemented yet
    if(!is.null(exclude)){
        l = lapply(l, function(x) x=x[!x%in%exclude])
    }

    require(dplyr)
    require(tidyr)

    # Get the combinations of names of list elements
    nms <- combn( names(l) , 2 , FUN = paste0 , collapse = separator , simplify = FALSE )

    # Make the combinations of list elements
    print("Making combinations")
    ll <- combn( l , 2 , simplify = FALSE )
    print(length(ll))
    # Intersect the list elements
    print("starting intersections")
    out <- lapply( ll , function(x) length( intersect( x[[1]] , x[[2]] ) ) )
    out = setNames( out , nms )
    return(out)

}
## And then to take the minimum number of pathway from the above intersections for gene sets I do the following
## Prepare the table
# geneSets2length = lapply(gene_sets, function(x) length(x))
# geneSets2length = data.frame(pathway=names(geneSets2length), pathway.length=unlist(geneSets2length))
# rownames(geneSets2length) = NULL
# out = lapply(seq_along(out), function(y,n,i){
#     ldply(y[[i]]) %>% mutate(pathway=n[[i]])
# }, y=out, n=names(out))
# out =do.call(rbind, out)
# out = out %>% separate(pathway, into=c("pathway1", "pathway2"), sep="\\.")
# out = out %>% left_join(geneSets2length%>%rename(pathway1=pathway, pathway1.length=pathway.length)) %>% left_join(geneSets2length%>%rename(pathway2=pathway, pathway2.length=pathway.length))
# out = out %>% mutate(pathway1.perc=(V1/pathway1.length)*100, pathway2.perc=(V1/pathway2.length)*100)
## Add also samples and genes

## Old function to get containers
getContainers2 = function(gsea=gsea.noAmp.top10.plusCGC.pathways[["genes2paths"]], cutoff=80, sample.filter=0){ ## d here is a data frame with all the necessary columns (like length of pathways, percentage of overlap etc)
    require(dplyr)

    gsea = gsea %>% subset(!is.na(pathway) & Level1==F & Level2==F & samples.pathway>sample.filter)
    cat(paste0("Merging ", length(unique(gsea$pathway)), " pathways"), "\n")

    uniq_paths = unique(gsea$pathway)


    ## Now I need to construct again the new similarity table
    similarity_table = combn(uniq_paths, 2) %>% t() %>% data.frame()
    cat(paste0("Pairs: ", nrow(similarity_table)), "\n")
    colnames(similarity_table) = c("pathway1", "pathway2")
    similarity_table$overlap = apply(similarity_table, 1, function(x){
        length(intersect(gsea %>% subset(pathway==x[1]) %>% .$symbol %>% unique, gsea %>% subset(pathway==x[2]) %>% .$symbol %>% unique))
    })

    unique_paths = unique(c(similarity_table$pathway1, similarity_table$pathway2))
    additional_info = NULL
    for(p in unique_paths){
        s = gsea %>% subset(pathway==p) %>% .$sample %>% unique %>% length
        l = gsea %>% subset(pathway==p) %>% .$symbol %>% unique %>% length
        d = data.frame(pathway=p, samples=s, length=l)
        additional_info = rbind(additional_info, d)
    }

    ## bring them in and calculate overlap
    similarity_table = similarity_table %>% left_join(additional_info %>% rename(pathway1=pathway, length1=length, samples.pathway1=samples)) %>%
        left_join(additional_info %>% rename(pathway2=pathway, length2=length, samples.pathway2=samples)) %>%
        mutate(pathway1.perc=(overlap/length1)*100) %>%
        mutate(pathway2.perc=(overlap/length2)*100)

    ## Here, I need to start the iterative processes
    ## Reset collapsed pathways
    containers = list()
    for(p in unique_paths){ ## Meaning that there is still a combination of pathways/containers
        d = similarity_table %>% subset((pathway1==p | pathway2==p) & (pathway1.perc>=cutoff | pathway2.perc>=cutoff))
        if(nrow(d)>0){
            d = data.frame(pathway=unique(c(d$pathway1, d$pathway2))) %>% left_join(additional_info) %>% arrange(desc(samples))
            container = d$pathway[1]
            if(container %in% names(containers)){ ## Bring in new pathways to the container
                containers[[container]] = unique(c(containers[[container]], p))
            }else{
                containers[[container]] = d$pathway
            }

        }else{
            containers[[p]] = p ## Add those that do not have similarity with anything
        }

    }

    ## return the two lists
    return(list(containers=containers, similarity_table=similarity_table))

}

getContainers = function(gsea=gsea.withAmp.top10.plusCGC.pathways[["genes2paths"]], cutoff=c(60), weight=T, sample.filter=0){ ## d here is a data frame with all the necessary columns (like length of pathways, percentage of overlap etc)
    require(plyr)
    require(dplyr)
    require(tidyr)
    require(igraph)

    ## save file to export it
    export = gsea
    gsea = gsea %>% subset(!is.na(pathway) & Level1==F & Level2==F & samples.pathway>sample.filter)
    cat(paste0("Merging ", length(unique(gsea$pathway)), " pathways"), "\n")

    uniq_paths = unique(gsea$pathway)


    ## Now I need to construct again the new similarity table
    cat("Calculating similarity matrix...", "\n")
    similarity_table = combn(uniq_paths, 2) %>% t() %>% data.frame()
    cat(paste0("Pairs: ", nrow(similarity_table)), "\n")
    colnames(similarity_table) = c("pathway1", "pathway2")
    similarity_table$overlap = apply(similarity_table, 1, function(x){
        length(intersect(gsea %>% subset(pathway==x[1] & !is.na(pathway)) %>% .$symbol %>% unique, gsea %>% subset(pathway==x[2] & !is.na(pathway)) %>% .$symbol %>% unique))
    })

    unique_paths = unique(c(as.character(similarity_table$pathway1), as.character(similarity_table$pathway2)))
    additional_info = NULL
    for(p in unique_paths){
        s = gsea %>% subset(pathway==p) %>% .$sample %>% unique %>% length
        l = gsea %>% subset(pathway==p) %>% .$symbol %>% unique %>% length
        d = data.frame(pathway=p, samples=s, length=l)
        additional_info = rbind(additional_info, d)
    }

    ## bring them in and calculate overlap
    similarity_table = similarity_table %>% left_join(additional_info %>% rename(pathway1=pathway, length1=length, samples.pathway1=samples)) %>%
        left_join(additional_info %>% rename(pathway2=pathway, length2=length, samples.pathway2=samples)) %>%
        mutate(pathway1.perc=(overlap/length1)*100) %>%
        mutate(pathway2.perc=(overlap/length2)*100)

    res = list()
    for(co in cutoff){

        net = similarity_table %>% subset(pathway1.perc>=co | pathway2.perc>=co) %>% select(pathway1, pathway2)
        w = similarity_table %>% subset(pathway1.perc>=co | pathway2.perc>=co) %>% mutate(weight=ifelse(length1<length2, pathway1.perc, ifelse(length2<=length1, pathway2.perc, NA))) %>% .$weight
        cat(paste0(nrow(net), " pairs with similarity >= ", co), "\n")
        g = graph_from_data_frame(net, directed = F, vertices = uniq_paths)
        E(g)$weight = w

        ## Get clusters by resolving bridges
        cat("Getting pathway clusters...", "\n")
        if(w){
            cls = cluster_edge_betweenness(g, weights = E(g)$weight)
        }else{
            cls = cluster_edge_betweenness(g)
        }


        cat(paste0(length(cls), " containers found"), "\n")
        print(sizes(cls))

        ## give the name of the container to the pathway covering more samples
        containers = list()
        for(i in 1:length(cls)){
            cont = gsea %>% select(pathway, samples.pathway) %>% unique %>% subset(pathway%in%cls[[i]]) %>% arrange(desc(samples.pathway))
            containers[[cont$pathway[1]]] = cls[[i]]
        }

        l = lapply(seq_along(containers), function(y,n,i){
            ldply(y[[i]]) %>% mutate(container=n[[i]])
        }, y=containers, n=names(containers))
        l = do.call(rbind, l) %>% rename(pathway=V1)

        ## join it with export and return it
        print(nrow(export))
        exp = export %>% left_join(l)
        print(nrow(export))

        res[[as.character(co)]] = list(containers=containers, clusters=cls, network=g, genes2paths=exp)
    }

    res[["similarity_table"]] = similarity_table
    return(res)

}

getPatientCommonPathways = function(genes2paths=gsea.withAmp.top10.plusCGC.containers[["genes2paths"]],
                                     normal = T
){

    require(dplyr)
    require(tidyr)

    d = genes2paths

    ## Decide if you want normalisation or not
    normal = normal

    d = d %>% subset(!is.na(pathway))
    ## Get a graph where nodes are samples and connections are number of containers
    sample_combns = do.call(rbind, combn(unique(d$sample), 2, simplify = FALSE))

    compaths = function(s1, s2, normal=T){
        if(s1==s2){
            stop("Sample twice")
        }
        common.pathways = d %>% select(sample, pathway) %>% unique %>% subset(sample%in%c(s1, s2)) %>% count(pathway) %>% subset(n>1) %>% .$pathway
        union_pathways = d %>% select(sample, pathway) %>% unique %>% subset(sample%in%c(s1, s2)) %>% select(pathway) %>% unique %>% nrow
        if(normal){
            common.pathways=length(common.pathways)/union_pathways
        }else{
            common.pathways=length(common.pathways)
        }
        return(common.pathways)
    }

    sample_combns = sample_combns %>% data.frame() %>% rename(sample1=X1, sample2=X2)
    cat("Calculating common pathways among samples...")
    sample_combns$common.pathways = apply(sample_combns, 1, function(x) compaths(x[1], x[2]))

    ## And add self versus self
    for (s in unique(d$sample)){
        common.pathways = d %>% subset(sample==s) %>% select(pathway) %>% unique %>% .$pathway
        if(normal){
            sample_combns = rbind(sample_combns, data.frame(sample1=s, sample2=s, common.pathways=1))
        }else{
            sample_combns = rbind(sample_combns, data.frame(sample1=s, sample2=s, common.pathways=length(common.pathways)))
        }
    }

    return(sample_combns)
}


getPatientCommonProcesses = function(genes2paths=gsea.withAmp.top10.plusCGC.containers[["genes2paths"]],
                                     normal = T
                                     ){

    require(dplyr)
    require(tidyr)

    d = genes2paths

    ## Decide if you want normalisation or not
    normal = normal

    d = d %>% subset(!is.na(container))
    ## Get a graph where nodes are samples and connections are number of containers
    sample_combns = do.call(rbind, combn(unique(d$sample), 2, simplify = FALSE))

    comconts = function(s1, s2, normal=T){
        if(s1==s2){
            stop("Sample twice")
        }
        common.containers = d %>% select(sample, container) %>% unique %>% subset(sample%in%c(s1, s2)) %>% count(container) %>% subset(n>1) %>% .$container
        union_processes = d %>% select(sample, container) %>% unique %>% subset(sample%in%c(s1, s2)) %>% select(container) %>% unique %>% nrow
        if(normal){
            common.containers=length(common.containers)/union_processes
        }else{
            common.containers=length(common.containers)
        }
        return(common.containers)
    }

    sample_combns = sample_combns %>% data.frame() %>% rename(sample1=X1, sample2=X2)
    cat("Calculating common processes among samples...")
    sample_combns$common.containers = apply(sample_combns, 1, function(x) comconts(x[1], x[2]))

    ## And add self versus self
    for (s in unique(d$sample)){
        common.containers = d %>% subset(sample==s) %>% select(container) %>% unique %>% .$container
        if(normal){
            sample_combns = rbind(sample_combns, data.frame(sample1=s, sample2=s, common.containers=1))
        }else{
            sample_combns = rbind(sample_combns, data.frame(sample1=s, sample2=s, common.containers=length(common.containers)))
        }
    }

    return(sample_combns)
}

## check the overlap of the containers (genes)
checkContainers = function(container_list=containers[["containers"]],
                           genes2paths=gsea.noAmp.top10.plusCGC.pathways[["genes2paths"]]){

    require(dplyr)
    require(tidyr)

    combs = combn(names(container_list), 2) %>% t() %>% data.frame()
    colnames(combs) = c("container1", "container2")

    getSharedGenes = function(container1, container2){
        paths1 = container_list[[container1]]
        paths2 = container_list[[container2]]

        genes1 = genes2paths %>% subset(pathway%in%paths1) %>% .$symbol %>% unique
        genes2 = genes2paths %>% subset(pathway%in%paths2) %>% .$symbol %>% unique

        return(paste0( length(intersect(genes1, genes2)), ";", length(genes1), ";", length(genes2) ))
    }

    combs$inter = apply(combs, 1, function(x) getSharedGenes(x[1], x[2]))
    combs = combs %>% separate(inter, into=c("overlap", "len1", "len2"), sep=";") %>% mutate(perc1=(as.numeric(overlap)/as.numeric(len1))*100, perc2=(as.numeric(overlap)/as.numeric(len2))*100)
}


## This function gets 2 entrez and returns the distance of the genes in hg19
geneDistanceHg19 <- function(g=NULL, gene1=NULL, gene2=NULL){
    ## Convert genes in strings (metadata column is string)
    if(is.null(g) | is.null(gene1) | is.null(gene2)){
        stop("Gene missing")
    }
    gene1 = as.character(gene1)
    gene2 = as.character(gene2)
    d = distance(g[names(g)==gene1], g[names(g)==gene2], ignore.strand=TRUE)
    return(d)
}

pathwayEnrichmentRanks = function(geneInfo_fn="/Users/fc-8s-imac/athena/data/geneInfoNCG5.Rdata",
                             gene_sets_dir="/Users/fc-8s-imac/athena/data/geneSets/reactome_12_12_16/",
                             gene_sets_gene_start=4,
                             gene_sets_identifiers = "reactome",
                             sub_identifier="symbol",
                             path_length=F,
                             min_path_length=10,
                             max_path_length=500,
                             scores_fn="/Users/fc-8s-imac/athena/data/OAC/129_OAC/ML/OAC/syscans_4_with_gains_in_sys_cans.Rdata",
                             cohort_fn="/Users/fc-8s-imac/athena/data/OAC/129_OAC/Rdata/mainCohort.Rdata",
                             remove.amplifications=T,
                             remove.cgc=T,
                             null_distribution=F,
                             sample.by.sample=F){

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

    ## x is a vector of genes
    ## y is a list of pathways
    ## identifier can be entrez/symbol
    updatePathway = function(x, y, identifier=sub_identifier, min_len=min_path_length, max_len=max_path_length, len=path_length,
                             exclude_pattern="CANCER|LEUKEMIA|MELANOMA|GLIOMA",
                             gi_fn=geneInfo_fn){

        ## This function removes all genes that are not in the gene list of interest
        ## from the pathways

        if (identifier=="entrez"){
            y = lapply(y, as.numeric)
        }else if (identifier=="symbol"){
            y = lapply(y, as.character)
        }

        ## Update for the list of genes that are in our gene set
        y = lapply(y, function(g) g[g%in%x])
        ## Remove empty sets
        y = y[unlist(lapply(y, function(x) length(x)>0))]

        if(len==TRUE){
            y = y[unlist(lapply(y, function(x) length(x)>min_len))]
            y = y[unlist(lapply(y, function(x) length(x)<max_len))]
        }
        y = y[names(y)[!grepl(exclude_pattern, names(y))]]

        message(paste0("Number of pathways: ", length(y)))
        if(identifier=="entrez"){
            genes_p = as.numeric(unique(unlist(y)))
        }else if(identifier=="symbol"){
            genes_p = unique(unlist(y))
        }
        message(paste0("Number of genes: ", length(genes_p)))
        ## Get gene info to check overlap with 19014
        load(gi_fn)
        geneInfo = geneInfo %>% subset(duplicability!=-1) %>% select(symbol, entrez)
        if(identifier=="entrez"){
            message(paste0("Number of genes (19,014): ", length(genes_p[genes_p%in%geneInfo$entrez])))
        }else if(identifier=="symbol"){
            message(paste0("Number of genes (19,014): ", length(genes_p[genes_p%in%geneInfo$symbol])))
        }


        return(list(gene_sets=y, genes_in_p=length(genes_p)))
    }

    ## x is the pathway
    ## g is the data frame with ranks
    run.singed.ks = function(x, g, func_fn="/Users/fc-8s-imac/athena/data/OAC/signed_ks.R"){
        source(func_fn)
        if(sub_identifier=="entrez"){
            s = g %>% subset(entrez%in%x) %>% select(rank) %>% .$rank
        }else if(sub_identifier=="symbol"){
            s = g %>% subset(symbol%in%x) %>% select(rank) %>% .$rank
        }

        k <- ks.test.2(s, (1:max(g$rank))[-s], maxCombSize=10^10)
        res = list(es=k$ES, p.value=k$p, no.genes=length(x), ranks=paste(s, collapse=","))
        return(res)
    }

    getNullDist = function(gene_sets, g, times=1000, ks_func_fn="/Users/fc-8s-imac/athena/data/OAC/signed_ks.R"){
        ## Get the null distribution
        source(ks_func_fn)
        l = unlist(lapply(gene_sets, function(x) length(x)))
        l = unique(l)
        ranks = c(1:max(g$rank))
        ex = replicate(times, sapply(l, function(x) {
            z = sample(ranks, x)
            t = ks.test.2(z, (1:max(g$rank))[-z], maxCombSize=10^10)
            ## For testing I return also the vector
            paste(t$ES, t$p, x, paste(z, collapse="/"), sep=",")
        }, simplify = F), simplify = F) %>% unlist()
        ex = data.frame(ex=ex) %>% separate(ex , into=c("es", "p.value", "len", "sample.ranks"), sep=",")
        return(ex)
    }

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
    if(remove.cgc){
        syscan = syscan %>% subset(gene_type!="cgc")
    }
    if(remove.amplifications){
        syscan = syscan %>% subset(CNVGain!=1)
    }

    ## get pathways
    gsea = NULL
    null_dists = NULL
    paths = gene_sets_identifiers

    if(sample.by.sample==F){
        message("Groups-collapsed mode")
        ## Get the genes for each group
        caDominant = syscan %>% subset(Group=="C>A/T dominant") %>% subset(gene_type!="CGC_discarded") %>%
            subset(type=="P") %>% select(entrez, symbol, score) %>% group_by(entrez, symbol) %>% summarise(score=mean(score)) %>% ungroup()
        addition = syscan %>% subset(Group=="C>A/T dominant") %>% subset(gene_type!="CGC_discarded") %>%
            subset(type=="C") %>% select(entrez, symbol) %>% unique %>% mutate(score=max(caDominant$score)+1)
        caDominant = rbind(caDominant, addition)
        caDominant = caDominant %>% arrange(desc(score)) %>% mutate(rank=row_number())
        if(sub_identifier=="entrez"){
            caDominant = caDominant %>% select(entrez, score, rank) %>% subset(!is.na(entrez))
        }else if(sub_identifier=="entrez"){
            caDominant = caDominant %>% select(symbol, score, rank) %>% subset(!is.na(symbol))
        }

        ddr = syscan %>% subset(Group=="DDR impaired") %>% subset(gene_type!="CGC_discarded") %>%
            subset(type=="P") %>% select(entrez, symbol, score) %>% group_by(entrez, symbol) %>% summarise(score=mean(score)) %>% ungroup()
        addition = syscan %>% subset(Group=="DDR impaired") %>% subset(gene_type!="CGC_discarded") %>%
            subset(type=="C") %>% select(entrez, symbol) %>% unique %>% mutate(score=max(ddr$score)+1)
        ddr = rbind(ddr, addition)
        ddr = ddr %>% arrange(desc(score)) %>% mutate(rank=row_number())
        if(sub_identifier=="entrez"){
            ddr = ddr %>% select(entrez, score, rank) %>% subset(!is.na(entrez))
        }else if(sub_identifier=="entrez"){
            ddr = ddr %>% select(symbol, score, rank) %>% subset(!is.na(symbol))
        }

        mutagenic = syscan %>% subset(Group=="Mutagenic") %>% subset(gene_type!="CGC_discarded") %>%
            subset(type=="P") %>% select(entrez, symbol, score) %>% group_by(entrez, symbol) %>% summarise(score=mean(score)) %>% ungroup()
        addition = syscan %>% subset(Group=="Mutagenic") %>% subset(gene_type!="CGC_discarded") %>%
            subset(type=="C") %>% select(entrez, symbol) %>% unique %>% mutate(score=max(mutagenic$score)+1)
        mutagenic = rbind(mutagenic, addition)
        mutagenic = mutagenic %>% arrange(desc(score)) %>% mutate(rank=row_number())
        if(sub_identifier=="entrez"){
            mutagenic = mutagenic %>% select(entrez, score, rank) %>% subset(!is.na(entrez))
        }else if(sub_identifier=="entrez"){
            mutagenic = mutagenic %>% select(symbol, score, rank) %>% subset(!is.na(symbol))
        }

        for(p in paths){
            ## C>A/T dominant
            gene_sets = getPathway(pathway = p)
            ## Just in case there are duplicated pathways
            dups = names(gene_sets)[duplicated(names(gene_sets))]
            for (d in dups){
                idx = which(names(gene_sets)==d)
                gene_sets[idx[2:length(idx)]] <- NULL
            }
            l = lapply(gene_sets, function(x) length(x))
            l = lapply(seq_along(l), function(y, n, i){
                ldply(y[[i]]) %>% mutate(pathway=n[[i]])
            }, y=l, n=names(l))
            l = do.call(rbind, l) %>% rename(pathway=V1)
            ## Fix names cause then we split on \\. and the folowing characters create problems
            # names(gene_sets) = gsub("\\.", "_", names(gene_sets))
            # names(gene_sets) = gsub("-", "_", names(gene_sets))
            # names(gene_sets) = gsub(" ", "_", names(gene_sets))
            # names(gene_sets) = gsub("\\'", "", names(gene_sets))
            # names(gene_sets) = gsub("\\/", "_", names(gene_sets))
            # names(gene_sets) = gsub("\\(", "", names(gene_sets))
            # names(gene_sets) = gsub("\\)", "", names(gene_sets))
            # names(gene_sets) = gsub(":", "_", names(gene_sets))
            # names(gene_sets) = gsub("&", "", names(gene_sets))
            # names(gene_sets) = gsub(",", "", names(gene_sets))

            if(sub_identifier=="entrez"){
                gene_sets_up = updatePathway(caDominant$entrez, gene_sets)
            }else if(sub_identifier=="symbol"){
                gene_sets_up = updatePathway(caDominant$symbol, gene_sets)
            }

            gene_sets_up = gene_sets_up[["gene_sets"]]

            if(null_distribution){
                ## Get the null distributions
                message("Getting null distribution")
                nullDist = getNullDist(gene_sets = gene_sets_up, g = caDominant)
                nullDist$group = "C>A/T dominant"
                null_dists = rbind(null_dists, nullDist)
            }


            res = lapply(gene_sets_up, function(x) run.singed.ks(x, caDominant))
            ## Make res data frame to work with
            res = lapply(seq_along(res), function(y, n, i){
                    ldply(y[[i]]) %>% mutate(pathway=n[[i]])
                }, y=res, n=names(res))
            res = do.call(rbind, res) %>% spread(.id, V1) %>% mutate(group="C>A/T dominant", gene_set=p)
            res$fdr = p.adjust(res$p.value, method="BH")
            res$bonf = p.adjust(res$p.value, method="bonferroni")
            res = res %>% left_join(l)
            gsea = rbind(gsea, res)

            ## DDR impaired
            if(sub_identifier=="entrez"){
                gene_sets_up = updatePathway(ddr$entrez, gene_sets)
            }else if(sub_identifier=="symbol"){
                gene_sets_up = updatePathway(ddr$symbol, gene_sets)
            }

            gene_sets_up = gene_sets_up[["gene_sets"]]
            ## Get the null distributions
            if(null_distribution){
                message("Getting null distribution")
                nullDist = getNullDist(gene_sets = gene_sets_up, g = ddr)
                nullDist$group = "DDR impaired"
                null_dists = rbind(null_dists, nullDist)
            }

            res = lapply(gene_sets_up, function(x) run.singed.ks(x, ddr))
            ## Make res data frame to work with
            res = lapply(seq_along(res), function(y, n, i){
                ldply(y[[i]]) %>% mutate(pathway=n[[i]])
            }, y=res, n=names(res))
            res = do.call(rbind, res) %>% spread(.id, V1) %>% mutate(group="DDR impaired", gene_set=p)
            res$fdr = p.adjust(res$p.value, method="BH")
            res$bonf = p.adjust(res$p.value, method="bonferroni")
            res = res %>% left_join(l)
            gsea = rbind(gsea, res)

            ## Mutagenic
            if(sub_identifier=="entrez"){
                gene_sets_up = updatePathway(mutagenic$entrez, gene_sets)
            }else if(sub_identifier=="symbol"){
                gene_sets_up = updatePathway(mutagenic$symbol, gene_sets)
            }

            gene_sets_up = gene_sets_up[["gene_sets"]]
            ## Get the null distributions
            if(null_distribution){
                message("Getting null distribution")
                nullDist = getNullDist(gene_sets = gene_sets_up, g = mutagenic)
                nullDist$group = "Mutagenic"
                null_dists = rbind(null_dists, nullDist)
            }

            res = lapply(gene_sets_up, function(x) run.singed.ks(x, mutagenic))
            res = lapply(seq_along(res), function(y, n, i){
                ldply(y[[i]]) %>% mutate(pathway=n[[i]])
            }, y=res, n=names(res))
            res = do.call(rbind, res) %>% spread(.id, V1) %>% mutate(group="Mutagenic", gene_set=p)
            res$fdr = p.adjust(res$p.value, method="BH")
            res$bonf = p.adjust(res$p.value, method="bonferroni")
            res = res %>% left_join(l)
            gsea = rbind(gsea, res)
        }

        return(list(gsea=gsea, null_dists=null_dists))
    }else if(sample.by.sample){
        ## Sample-by-sample pathway enrichment
        message("Sample-by-sample mode")
        for(p in paths){
            genes2samples = NULL
            genes_mapped = NULL
            gene_sets = getPathway(pathway = p)
            dups = names(gene_sets)[duplicated(names(gene_sets))]
            track_samples = NULL
            for (d in dups){
                idx = which(names(gene_sets)==d)
                gene_sets[idx[2:length(idx)]] <- NULL
            }
            l = lapply(gene_sets, function(x) length(x))
            l = lapply(seq_along(l), function(y, n, i){
                ldply(y[[i]]) %>% mutate(pathway=n[[i]])
            }, y=l, n=names(l))
            l = do.call(rbind, l) %>% rename(pathway_length=V1)
            for (s in unique(syscan$sample)){

                mg = syscan %>% subset(sample==s) %>% .$Group %>% unique

                g = syscan %>% subset(sample==s) %>% subset(gene_type!="CGC_discarded") %>%
                    select(entrez, symbol, score) %>% arrange(desc(score)) %>% mutate(rank=row_number())

                ## Check point for NAs in score column
                if(g%>%subset(is.na(score))%>%nrow>0){
                    stop("NAs found in sample scores")
                }
                if(sub_identifier=="entrez"){
                    gene_sets_up = updatePathway(g$entrez, gene_sets)
                }else if(sub_identifier=="symbol"){
                    gene_sets_up = updatePathway(g$symbol, gene_sets)
                }
                genes_in_p = gene_sets_up[["genes_in_p"]]
                gene_sets_up = gene_sets_up[["gene_sets"]]
                genes_mapped = c(genes_mapped, unique(unlist(gene_sets_up)))
                ## Keep a vector with genes in each sample
                genes2samples = c(genes2samples, genes_in_p)
                message(paste0(s, ": ", genes_in_p, " in gene sets"))

                if(genes_in_p==0){
                    next
                }

                ## Get the null distributions
                if(null_distribution){
                    message("Getting null distribution")
                    nullDist = getNullDist(gene_sets = gene_sets_up, g = g)
                    nullDist = nullDist %>% mutate(sample=s, group=mg, gene_set=p)
                    null_dists = rbind(null_dists, nullDist)
                }

                res = lapply(gene_sets_up, function(x) run.singed.ks(x, g))
                res = lapply(seq_along(res), function(y, n, i){
                    ldply(y[[i]]) %>% mutate(pathway=n[[i]])
                }, y=res, n=names(res))
                res = do.call(rbind, res) %>% spread(.id, V1) %>% mutate(sample=s, group=mg, gene_set=p)
                res$fdr = p.adjust(res$p.value, method="BH")
                res$bonf = p.adjust(res$p.value, method="bonferroni")
                res = res %>% left_join(l)
                gsea = rbind(gsea, res)

                ## Populate track samples table to check the results
                print(s)
                print(genes_in_p)
                print(unique(g$Group))
                print(p)
                ts  = data.frame(sample=s, number.of.genes=nrow(g), genes.in.pathways=genes_in_p, group=mg, geneset=p)
                track_samples = rbind(track_samples, ts)
                cat("\n")
            }
            message(paste0("Summary of genes in samples: ", p))
            print(summary(genes2samples))
            message(paste0("Unique genes mapped across all samples: ", length(unique(genes_mapped))))
        }
        return(list(gsea=gsea, null_dists=null_dists, track_samples=track_samples))
    }
}

## use this function to map samples to enriched pathways
mapPathways2Samples = function(gsea, efdr=0.01,
                               geneInfo_fn="~/Mountpoints/rosalind_lustre/mourikisa/data/geneInfoNCG5.Rdata",
                               gene_sets_dir="~/Mountpoints/rosalind_lustre/mourikisa/data/geneSets/reactome_12_12_16/",
                               gene_sets_gene_start=4,
                               gene_sets_identifiers = "reactome",
                               sub_identifier="symbol"){


    ## Get genes in pathways
    gene_sets = gsea[["paths"]]

    ## Get genes and pathways from the gsea
    enriched_paths = gsea[["gsea"]] %>% subset(fdr<efdr) %>% .$pathway %>% unique
    cat(paste0("Enriched pathways: ", length(enriched_paths)))

    ## Get only enriched pathways without
    gene_sets = gene_sets[names(gene_sets)%in%enriched_paths]

    ## Get gene sets ready for left joining
    l = lapply(seq_along(gene_sets), function(y, n, i){
        ldply(y[[i]]) %>% mutate(pathway=n[[i]])
    }, y=gene_sets, n=names(gene_sets))
    l = do.call(rbind, l) %>% rename(symbol=V1)

    ## Add pathways in genes
    genes = gsea[["genes"]] %>% left_join(l)

    paths2s2g = genes %>% subset(!is.na(pathway)) %>% group_by(pathway) %>% summarise(samples.pathway=length(unique(sample)), genes.pathway=length(unique(symbol))) %>% ungroup()

    paths2s2g = paths2s2g %>% left_join(gsea[["gsea"]]%>%select(pathway, Level1, Level2))
    ## Add also Level in reactome
    genes = genes %>% left_join(paths2s2g)

    ## return the enriched gsea list of objects
    return(list(gsea=gsea[["gsea"]], genes=gsea[["genes"]], genes2paths=genes, paths2s2g=paths2s2g))
}

## Use this function to find minimal set of processes explain each group
getREACTOMEsubgraph = function(hierarchy.file="~/athena/data/geneSets/reactome_12_12_16/pathway_hierarchy.tsv",
                               identifier.file="~/athena/data/geneSets/reactome_12_12_16/complete_list_of_pathways.tsv",
                               pathways.enriched="~/athena/enriched_paths2samples_CGC_plus_top10syscans_noGroups_fdr_0.01_pathwayFilters_corrected.tsv"){

    require(igraph)
    require(xlsx)
    require(gplots)
    require(RColorBrewer)
    require(RedeR)
    require(ComplexHeatmap)


    hier = read.table(hierarchy.file, sep = "\t", header = F)
    colnames(hier) = c("path1", "path2")
    ## Careful here cause some pathways have apostrophes in their names
    path_ids = read.table(identifier.file, sep = "\t", header = F, quote = "")
    colnames(path_ids) = c("path_id", "path_name", "species")
    path_ids = path_ids %>% subset(species=="Homo sapiens") %>% rename(pathway=path_name) %>% select(-species)

    ## Exclude leading nad trailing spaces from pathway names
    path_ids$pathway = gsub("^\\s+|\\s+$", "", path_ids$pathway)
    path_ids$path_id = gsub("^\\s+|\\s+$", "", path_ids$path_id)

    ## Add info to herarchy table and subset for the homo sapiens interactions
    hier = hier %>% subset(path1%in%path_ids$path_id | path2%in%path_ids$path_id) %>% left_join(path_ids%>%select(path_id, pathway)%>%rename(path1=path_id, pathway1=pathway)) %>%
        left_join(path_ids%>%select(path_id, pathway)%>%rename(path2=path_id, pathway2=pathway))

    ## Trials to exclude level 1&2 nodes
    level_1_nodes = c("Developmental Biology",
                      "Reproduction",
                      "Circadian Clock",
                      "Hemostasis",
                      "Neuronal System",
                      "Immune System",
                      "Signal Transduction",
                      "Disease",
                      "DNA Repair",
                      "Metabolism",
                      "Gene Expression",
                      "Chromatin organization",
                      "Transmembrane transport of small molecules",
                      "DNA replication",
                      "Metabolism of proteins",
                      "Muscle contraction",
                      "Cell Cycle",
                      "Organelle biogenesis and maintenance",
                      "Mitophagy",
                      "Vesicle-mediated transport",
                      "Programmed Cell Death",
                      "Cellular responses to stress",
                      "Extracellular matrix organization",
                      "Cell-Cell communication")

    level_2_nodes = hier %>% subset(pathway1%in%level_1_nodes | pathway2%in%level_1_nodes)
    level_2_nodes = unique(c(level_2_nodes$pathway1, level_2_nodes$pathway2))
    level_2_nodes = level_2_nodes[!level_2_nodes%in%level_1_nodes]


    ## Read in the group-specific enriched pathways
    gsep = read.table(pathways.enriched, sep="\t", header = T, quote="")

    ## fix name quotes etc, rgsea.counts comes from randomisation
    gsep$pathway = gsub("\"", "", gsep$pathway)
    gsep$pathway = gsub("\'", "", gsep$pathway)
    rgsea.counts$pathway = gsub("\'", "", rgsea.counts$pathway)
    rgsea.counts = rgsea.counts %>% mutate(inCGCplusTop10Enriched=ifelse(pathway%in%unique(gsep$pathway), TRUE, FALSE))


    paths2samples2genes = gsep %>% group_by(pathway) %>% summarise(samples=length(unique(sample)), genes=length(unique(symbol)))
    gsep = gsep %>% left_join(paths2samples2genes) %>% mutate(label=paste0(pathway, "(", samples, ";", genes, ")"), inGS=1)

    ## Refinement of pathways
    gsep = gsep %>% subset(!pathway%in%level_1_nodes & !pathway%in%level_2_nodes)
    gsep = gsep %>% subset(samples>11)

    ## Refine based on randomisation
    truly_enriched = rgsea.counts %>% subset(ns_perc>95 & inCGCplusTop10Enriched==TRUE) %>% .$pathway
    gsep = gsep %>% subset(pathway%in%truly_enriched)

    ## Extract subnetwork
    subnet = hier %>% subset(pathway1%in%unique(gsep$pathway) & pathway2%in%unique(gsep$pathway))
    #subnet = hier %>% subset(path1%in%gsep_ids | path2%in%gsep_ids)

    ## Bring names for path1 & path2
    subnet = subnet %>% left_join(gsep%>%select(pathway, label)%>%unique%>%rename(pathway1=pathway, label1=label)) %>%
        left_join(gsep%>%select(pathway, label)%>%unique%>%rename(pathway2=pathway, label2=label))


    ## this is for ehn we select interactions with both nodes to be pathways enriched
    lonely_nodes = gsep %>% subset(!label%in%subnet$label1 & !label%in%subnet$label2) %>% select(label) %>% unique

    ## Number of nodes
    length(unique(c(subnet$label1, subnet$label2)))
    length(unique(c(subnet$path1, subnet$path2)))

    ## Matrix for plotting
    toPlot = gsep %>% select(label, sample) %>% unique %>% table()

    # toPlot = gsep %>% select(label, sample) %>% subset(grepl("ER to Golgi Antero", label) | grepl("Signaling by Robo receptor", label) |
    #                                               grepl("RAF/MAP kinase cascade", label) | grepl("Golgi-to-ER retrograde transport", label) |
    #                                               grepl("Interleukin-3, 5 and GM", label) | grepl("Ion channel transport", label) |
    #                                               grepl("N-glycan antennae elongation in the medial/trans-Golgi", label) |
    #                                               grepl("PIP3 activates AKT signaling", label) | grepl("Rap1 signalling", label ) |
    #                                               grepl("Translation initiation complex formation", label)) %>% unique %>% table()

    ## Plot the network and color the group-specif enriched pathways
    ## http://kateto.net/network-visualization
    nodes = data.frame(pathway=unique(c(subnet$label1, subnet$label2, lonely_nodes$label)))%>%mutate(inGS=ifelse(pathway%in%gsep$label, 1, 0))
    #nodes = data.frame(pathway=unique(c(subnet$label1, subnet$label2)))%>%mutate(inGS=ifelse(pathway%in%gsep$label, 1, 0))
    net = graph_from_data_frame(d=subnet%>%select(label1, label2), vertices=nodes,directed=F)
    cls = groups(components(net))
    cls = lapply(seq_along(cls), function(y, n, i){
        ldply(y[[i]]) %>% mutate(cluster=n[[i]])
    }, y=cls, n=names(cls))
    cls = do.call(rbind, cls)

    cls = cls %>% subset(V1%in%rownames(toPlot))
    rnames = cls$V1
    cls = cls %>% select(-V1)
    rownames(cls) = rnames

    ## For each cluster get the unique number of samples and genes
    ## Then exclude those samples and cluster and redo the same
    samples_covered = NULL
    sequence_of_clusters = NULL
    cls_to_check = cls %>% tibble::rownames_to_column() %>% mutate(cluster = as.numeric(cluster))## starting point, this will be updated in the while loop
    gsep_to_check = gsep
    #max_samples = gsep %>% subset(!pathway%in%exclusion_list) %>% select(sample) %>% unique %>% nrow
    max_samples = gsep %>% select(sample) %>% unique %>% nrow

    while(length(samples_covered)<max_samples){
        ## select the cluster with the biggest coverage of samples
        biggest_cluster = NULL
        samples_so_far = NULL
        for(i in unique(cls_to_check$cluster)){
            paths = cls_to_check %>% subset(cluster==i) %>% .$rowname
            samples = gsep_to_check %>% subset(label%in%paths) %>% select(sample) %>% unique %>% .$sample
            if(length(samples)>length(samples_so_far)){
                biggest_cluster = i
                samples_so_far = samples
            }else{
                next
            }
        }

        sequence_of_clusters = c(sequence_of_clusters, biggest_cluster)
        samples_covered = c(samples_covered, samples_so_far)
        print(length(samples_covered))
        ## Now exclude that samples from the gsep
        gsep_to_check = gsep_to_check %>% subset(!(sample%in%samples_covered))
        ## Update cls_to_check by removing the bigest cluster
        cls_to_check = cls_to_check %>% subset(cluster!=biggest_cluster)
    }


    # Generate colors based on media type:
    colrs <- c("gray50", "tomato")
    V(net)$color = ifelse(V(net)$inGS==1, "tomato", "gray50")

    #plot(net, vertex.size=7, vertex.label.cex=0.2)
    tkid <- tkplot(g) #tkid is the id of the tkplot that will open
    ll <- tkplot.getcoords(tkid) # grab the coordinates from tkplot
    plot(g, layout=ll,vertex.size=5, vertex.label.cex=1)

    ## Cluster patients based on the group-specific pathways
    toPlot = toPlot[match(rownames(cls), rownames(toPlot)), ]
    ## Make it table cause somthing wierd is going on with the names if the class of object is table
    attributes(toPlot)$class <- "matrix"


    ## Get all possible intersections of the binary matrix (toPlot)
    toPlot_list = list()
    for(c in colnames(toPlot)){
        toPlot_list[[c]] = names(toPlot[,c][toPlot[,c]==1])
    }


    ## Create annotation of samples
    # sample_ann = mainCohort %>% select(sample, specimen_donor_treatment_type,
    #                                    tumour_grade, donor_relapse_interval,
    #                                    donor_age_at_diagnosis, donor_survival_time, Group)
    sample_ann = mainCohort %>% select(sample,Group)
    sample_ann[sample_ann==""] = NA
    sample_ann = sample_ann %>% subset(sample%in%colnames(toPlot))
    sample_ann = sample_ann[match(sample_ann$sample, colnames(toPlot)), ]


    # ## Create cluser table based on top 10 plus CGC
    # path2clust = gsep %>% select(label, process) %>% unique
    # rnames = path2clust$label
    # path2clust = data.frame(process=path2clust$process)
    # rownames(path2clust) = rnames

    pdf(file = "~/Desktop/all_129_enriched_top_10.pdf", width = 15, height = 15)
    heatmap.2(toPlot, trace="none", col=c("grey50", "red"),margins=c(10,10), cexRow=0.1, key = FALSE, Rowv = FALSE, dendrogram = "column",
              main="All 129 top 15 gene enrichment", RowSideColors=cls$cluster, ColSideColors = sample_ann$col)
    dev.off()


    pdf(file = "~/Desktop/test.pdf", width = 25, height = 15)
    #cls.cols = c(brewer.pal(8, "Accent"), brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(12, "Set3"), brewer.pal(5, "Set1"))
    cls.cols = rainbow(length(unique(cls$cluster)))

    ha = HeatmapAnnotation(df = data.frame(group = sample_ann$Group),col = list(group = c("Mutagenic" =  "green", "DDR impaired" = "purple", "C>A/T dominant"="orange")))
    # ha = HeatmapAnnotation(df = sample_ann %>% select(-sample) %>% data.frame(), col = list(Group = c("Mutagenic" =  "green", "DDR impaired" = "purple", "C>A/T dominant"="orange"),
    #                                                                                         tumour_grade=c("3"="black", "2"="red", "1"="blue", "unknown"="grey50", "needs cpath clarification"="green")))
    Heatmap(toPlot,
            row_names_gp = gpar(fontsize=6),
            column_names_gp = gpar(fontsize=7),
            col=c("white", "red"),
            top_annotation = ha,
            #cluster_rows = FALSE,
            #cluster_columns = FALSE,
            #column_order = samples2containers$sample,
            #row_order = gsep_plus_containers %>% select(label, samples.container) %>% unique %>% arrange(desc(samples.container)) %>% .$label,
            column_title = 'GSEA on 120; No Groups considered; Top 10 genes plus CGC with amplifications; \n No level 1/2 reactome pathways; No pathways with less than 12 samples; Pathway containers') +
    Heatmap(cls$cluster, name="Pathway cluster",col = cls.cols, width = unit(5, "mm"), heatmap_legend_param = list(ncol=6))
    dev.off()



    ## Plot/cluster the network using RedeR package
    check = l[["container_overlap"]] %>% subset(V1>3)
    net = graph_from_data_frame(d=check%>%select(container1, container2),directed=F)
    E(net)$weight = check$V1
    plot(net, edge.width=check$V1, vertex.size=2)
    rdp=RedPort()
    calld(rdp)
    addGraph(rdp, net)
    relax(rdp)


    ## ----------------------------
    ## Plot the network as heatmap
    ## ----------------------------

    # # Re-generate dataframes for both nodes and edges, now containing
    # # calculated network attributes
    # node_list <- get.data.frame(net, what = "vertices")
    #
    # # Determine a community for each edge. If two nodes belong to the
    # # same community, label the edge with that community. If not,
    # # the edge community value is 'NA'
    # edge_list <- get.data.frame(net, what = "edges")
    # # Create a character vector containing every node name
    # all_nodes <- sort(node_list$name)
    #
    # # Adjust the 'to' and 'from' factor levels so they are equal
    # # to this complete list of node names
    # plot_data <- edge_list %>% mutate(
    #     to = factor(to, levels = all_nodes),
    #     from = factor(from, levels = all_nodes))
    #
    # # Create the adjacency matrix plot
    # ggplot(plot_data, aes(x = from, y = to, fill = group)) +
    #     geom_raster() +
    #     theme_bw() +
    #     # Because we need the x and y axis to display every node,
    #     # not just the nodes that have connections to each other,
    #     # make sure that ggplot does not drop unused factor levels
    #     scale_x_discrete(drop = FALSE) +
    #     scale_y_discrete(drop = FALSE) +
    #     theme(
    #         # Rotate the x-axis lables so they are legible
    #         axis.text.x = element_text(angle = 270, hjust = 0),
    #         # Force the plot into a square aspect ratio
    #         aspect.ratio = 1,
    #         # Hide the legend (optional)
    #         legend.position = "none")

    load("~/athena/239_OAC_enriched_pathways_plus_containers_80percent.Rdata")

    ## Decide if you want normalisation or not
    normal = T
    normal = F

    ## Exclude Tp53
    gsep_plus_containers = gsep_plus_containers %>% subset(symbol!="TP53")

    ## Get only CGCs
    gsep_plus_containers = gsep_plus_containers %>% subset(gene_type=="cgc")

    ## Exclude CGCs
    gsep_plus_containers = gsep_plus_containers %>% subset(gene_type!="cgc")

    ## Get a graph where nodes are samples and connections are number of containers
    sample_combns = do.call(rbind, combn(unique(gsep_plus_containers$sample), 2, simplify = FALSE))
    net = NULL
    for(i in 1:nrow(sample_combns)){
        cat(i, "\n")
        if(sample_combns[i,][1]==sample_combns[i,][2]){
            stop("Sample twice")
        }
        common.containers = gsep_plus_containers %>% select(sample, container) %>% unique %>% subset(sample%in%sample_combns[i,]) %>% count(container) %>% subset(n>1) %>% .$container
        union_processes = gsep_plus_containers %>% select(sample, container) %>% unique %>% subset(sample%in%sample_combns[i,]) %>% select(container) %>% unique %>% nrow
        if(normal){
            net = rbind(net, data.frame(sample1=sample_combns[i,][1], sample2=sample_combns[i,][2], common.containers=length(common.containers)/union_processes))
        }else{
            net = rbind(net, data.frame(sample1=sample_combns[i,][1], sample2=sample_combns[i,][2], common.containers=length(common.containers)))
        }
    }
    ## And add self versus self
    for (s in unique(gsep_plus_containers$sample)){
        common.containers = gsep_plus_containers %>% subset(sample==s) %>% select(container) %>% unique %>% .$container
        if(normal){
            net = rbind(net, data.frame(sample1=s, sample2=s, common.containers=1))
        }else{
            net = rbind(net, data.frame(sample1=s, sample2=s, common.containers=length(common.containers)))
        }
    }

    ## Get symmetric matrix
    g <- graph.data.frame(net%>% subset(common.containers>0), directed=FALSE)
    net2 = get.adjacency(g, attr="common.containers", sparse=FALSE)
    library(RColorBrewer)
    library(circlize)
    ## I select the breaks based on the distribution
    up = unname(quantile(c(net2), seq(0, 1, .05))[20])
    p =Heatmap(net2,colorRamp2(breaks = c(seq(0,up, up/10)),
                               colors = c("#053061","#2166ac","#4393c3","#92c5de","#d1e5f0","#f7f7f7","#fddbc7","#f4a582","#d6604d","#f03b20","#bd0026")),
               heatmap_legend_param = list(title="Shared processes/union", color_bar="continuous"),
               row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5))

    net2 = net2[unlist(row_dend(p)),unlist(column_dend(p))]

    p =Heatmap(net2,colorRamp2(breaks = c(seq(0,up, up/10)),
                               colors = c("#053061","#2166ac","#4393c3","#92c5de","#d1e5f0","#f7f7f7","#fddbc7","#f4a582","#d6604d","#f03b20","#bd0026")),
               heatmap_legend_param = list(title="Shared processes/union", color_bar="continuous"),
               row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
               cell_fun = function(j, i, x, y, w, h, col) {
                   if((i+j)>(nrow(net2))) {
                       grid.rect(x, y, w, h, gp = gpar(fill = "white", col="white"))
                   }
               }
               )


    ## Check groups
    ro = rownames(net2)[unlist(row_order(p))]
    high_sharing = ro[1:which(ro=="LP6005500-DNA_B03")]
    intermediate_sharing = ro[(which(ro=="LP6005500−DNA_B03")+1):which(ro=="LP6005334-DNA_G02")]
    no_sharing = ro[(which(ro=="LP6005500-DNA_B03")+1):length(ro)]

    t = gsep_plus_containers %>% select(sample, container) %>% unique %>% mutate(group=ifelse(sample%in%high_sharing, "high_sharing",
                                                                                              ifelse(sample%in%intermediate_sharing, "intermediate_sharing",
                                                                                                     ifelse(sample%in%no_sharing, "no_sharing",
                                                                                                            ifelse(sample%in%low_sharing, "low_sharing",NA))))) %>%
        group_by(group, container) %>% summarise(n=n()) %>% ungroup() %>% spread(group, n)


    ## Then do a proportion test on the above table
    cats = 2
    t = t_syscans %>% data.frame()
    total = data.frame(type=c("high_sharing", "low_sharing"), total=c(61, 59))
    stats=NULL
    for(i in 1:nrow(t)){
        cat(i, "\t")
        d = t[i,1:(cats+1)] %>% gather(type, value, -container)
        d[is.na(d)] = 0
        d = d %>% left_join(total)
        combs = combn(d$type, 2, simplify=F)
        for (y in 1:length(combs)){
            cat(y, "\t")
            dd = d %>% subset(type%in%combs[[y]])
            pt = prop.test(dd$value, dd$total)
            tb = data.frame(container=unique(dd$container), group1=dd$type[1], group2=dd$type[2],
                            group1_value=dd$value[1], group1_prop=pt$estimate[1], group2_value=dd$value[2],
                            group2_prop=pt$estimate[2], p.value=pt$p.value)
            rownames(tb) = NULL
            stats = rbind(stats, tb)
        }
        cat("\n")
    }

    ## create the table
    ## Get row and column order from Complex heatmap package
    ro = rownames(net2)[unlist(row_order(p))]
    co = rownames(net2)[rev(unlist(column_order(p)))]
    tb = net %>% subset(common.containers.number>0)
    ## reset the rownames
    rownames(tb) = 1:nrow(tb)
    orders = NULL
    for(i in 2:120){
        for(y in 1:(i-1)){
            orders <- c(orders,list(c(ro[i],co[y])))
        }
    }
    orders_final = lapply(orders, function(x){
                            r = rownames(tb %>% subset((sample1==x[1] & sample2==x[2]) | sample1==x[2] & sample2==x[1]))
                            })
    orders_final = unlist(orders_final)
    tb = tb[match(as.numeric(orders_final),rownames(tb)),]
    rownames(tb) = 1:nrow(tb)
    write.table(tb, file="Desktop/sample_shared_containers_network.tsv", sep="\t", quote=F)

    ## It seems that there are three groups (groups extracted by looking the heatmap)
    ## Extract the groups
    ## 1st group
    high_sharing = co[1:which(co=="LP6005500-DNA_H03")]
    limited_sharing = co[(which(co=="LP6005500-DNA_H03")+1):which(co=="LP6005935-DNA_C04")]
    intermediate_sharing = co[(which(co=="LP6005935-DNA_C04")+1):which(co=="SS6003119")]
    no_sharing = co[(which(co=="SS6003119")+1):length(co)]

    ## Number of processes per sample
    samples2info = rbind(gsep_plus_containers %>% subset(sample%in%high_sharing) %>% select(sample, container) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="high-sharing", info_type="containers"),
        gsep_plus_containers %>% subset(sample%in%intermediate_sharing) %>% select(sample, container) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="intermediate-sharing", info_type="containers"),
        gsep_plus_containers %>% subset(sample%in%limited_sharing) %>% select(sample, container) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="limited-sharing", info_type="containers"),
        gsep_plus_containers %>% subset(sample%in%no_sharing) %>% select(sample, container) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="no-sharing", info_type="containers"))

    samples2info = rbind(samples2info, gsep_plus_containers %>% subset(sample%in%high_sharing) %>% select(sample, pathway) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="high-sharing", info_type="pathways"),
                         gsep_plus_containers %>% subset(sample%in%intermediate_sharing) %>% select(sample, pathway) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="intermediate-sharing", info_type="pathways"),
                         gsep_plus_containers %>% subset(sample%in%limited_sharing) %>% select(sample, pathway) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="limited-sharing", info_type="pathways"),
                         gsep_plus_containers %>% subset(sample%in%no_sharing) %>% select(sample, pathway) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="no-sharing", info_type="pathways"))

    ## When you check the CGCs there are samples with 0s so we have to include them
    high_cgcs = gsea.top.10.plus.cgc[["genes"]] %>% subset(sample%in%high_sharing & gene_type=="cgc") %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="high-sharing", info_type="CGC genes")
    inter_cgcs = gsea.top.10.plus.cgc[["genes"]] %>% subset(sample%in%intermediate_sharing & gene_type=="cgc") %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="intermediate-sharing", info_type="CGC genes")
    limit_cgcs = gsea.top.10.plus.cgc[["genes"]] %>% subset(sample%in%limited_sharing & gene_type=="cgc") %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="limited-sharing", info_type="CGC genes")
    no_cgcs = gsea.top.10.plus.cgc[["genes"]] %>% subset(sample%in%no_sharing & gene_type=="cgc") %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="no-sharing", info_type="CGC genes")
    missing_samples = no_sharing[!no_sharing%in%no_cgcs$sample]
    no_cgcs = rbind(no_cgcs, data.frame(sample=missing_samples, value=rep(0, length(missing_samples)), sample_type=rep("no-sharing", length(missing_samples)), info_type=rep("CGC genes", length(missing_samples))))


    samples2info = rbind(samples2info,
                         high_cgcs,
                         inter_cgcs,
                         limit_cgcs,
                         no_cgcs
                         )

    ## Genes mapped in enriched pathways
    samples2info = rbind(samples2info, gsep_plus_containers %>% subset(sample%in%high_sharing) %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="high-sharing", info_type="Genes in enriched pathways"),
                         gsep_plus_containers %>% subset(sample%in%intermediate_sharing) %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="intermediate-sharing", info_type="Genes in enriched pathways"),
                         gsep_plus_containers %>% subset(sample%in%limited_sharing) %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="limited-sharing", info_type="Genes in enriched pathways"),
                         gsep_plus_containers %>% subset(sample%in%no_sharing) %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="no-sharing", info_type="Genes in enriched pathways"))

    ## CGCs mapped in enriched pathways
    samples2info = rbind(samples2info, gsep_plus_containers %>% subset(gene_type=="cgc") %>% subset(sample%in%high_sharing) %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="high-sharing", info_type="CGCs in enriched pathways"),
                         gsep_plus_containers %>% subset(gene_type=="cgc") %>% subset(sample%in%intermediate_sharing) %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="intermediate-sharing", info_type="CGCs in enriched pathways"),
                         gsep_plus_containers %>% subset(gene_type=="cgc") %>% subset(sample%in%limited_sharing) %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="limited-sharing", info_type="CGCs in enriched pathways"))
    ## Add the no-sharing (there are 3 samples with not CGCs mutated)
    no_cgcs = gsep_plus_containers %>% subset(gene_type=="cgc") %>% subset(sample%in%no_sharing) %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="no-sharing", info_type="CGCs in enriched pathways")
    missing_samples = no_sharing[!no_sharing%in%no_cgcs$sample]
    no_cgcs = rbind(no_cgcs, data.frame(sample=missing_samples, value=rep(0, length(missing_samples)), sample_type=rep("no-sharing", length(missing_samples)), info_type=rep("CGCs in enriched pathways", length(missing_samples))))
    samples2info = rbind(samples2info, no_cgcs)

    ## sys-candidates mapped in enriched pathways
    samples2info = rbind(samples2info, gsep_plus_containers %>% subset(gene_type!="cgc") %>% subset(sample%in%high_sharing) %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="high-sharing", info_type="sys-candidates in enriched pathways"),
                         gsep_plus_containers %>% subset(gene_type!="cgc") %>% subset(sample%in%intermediate_sharing) %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="intermediate-sharing", info_type="sys-candidates in enriched pathways"),
                         gsep_plus_containers %>% subset(gene_type!="cgc") %>% subset(sample%in%limited_sharing) %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="limited-sharing", info_type="sys-candidates in enriched pathways"),
                         gsep_plus_containers %>% subset(gene_type!="cgc") %>% subset(sample%in%no_sharing) %>% select(sample, symbol) %>% unique %>% group_by(sample) %>% summarise(value=n()) %>% mutate(sample_type="no-sharing", info_type="sys-candidates in enriched pathways"))


    ## Fix the order of the facets
    samples2info$info_type = factor(samples2info$info_type, levels=c("containers","pathways","CGC genes", "CGCs in enriched pathways", "Genes in enriched pathways", "sys-candidates in enriched pathways"))
    ggplot(samples2info, aes(x=sample_type, y=value, fill=sample_type)) + geom_boxplot() + facet_wrap(~info_type, ncol = 2, scales = "free")


    ## Check if on average CGCs, sys-candidates map to more/less pathways in each group
    ## CGCs mapped in enriched pathways
    genes2info = rbind(gsep_plus_containers %>% subset(gene_type=="cgc") %>% subset(sample%in%high_sharing) %>% select(symbol, container) %>% unique %>% group_by(symbol) %>% summarise(value=n()) %>% mutate(sample_type="high-sharing", info_type="Enriched pathways per CGCs"),
                         gsep_plus_containers %>% subset(gene_type=="cgc") %>% subset(sample%in%intermediate_sharing) %>% select(symbol, container) %>% unique %>% group_by(symbol) %>% summarise(value=n()) %>% mutate(sample_type="intermediate-sharing", info_type="Enriched pathways per CGCs"),
                         gsep_plus_containers %>% subset(gene_type=="cgc") %>% subset(sample%in%limited_sharing) %>% select(symbol, container) %>% unique %>% group_by(symbol) %>% summarise(value=n()) %>% mutate(sample_type="limited-sharing", info_type="Enriched pathways per CGCs"),
                         gsep_plus_containers %>% subset(gene_type=="cgc") %>% subset(sample%in%no_sharing) %>% select(symbol, container) %>% unique %>% group_by(symbol) %>% summarise(value=n()) %>% mutate(sample_type="no-sharing", info_type="Enriched pathways per CGCs"))


    ## sys-candidates mapped in enriched pathways
    genes2info = rbind(genes2info, gsep_plus_containers %>% subset(gene_type!="cgc") %>% subset(sample%in%high_sharing) %>% select(symbol, container) %>% unique %>% group_by(symbol) %>% summarise(value=n()) %>% mutate(sample_type="high-sharing", info_type="Enriched pathways per sys-candidates"),
                         gsep_plus_containers %>% subset(gene_type!="cgc") %>% subset(sample%in%intermediate_sharing) %>% select(symbol, container) %>% unique %>% group_by(symbol) %>% summarise(value=n()) %>% mutate(sample_type="intermediate-sharing", info_type="Enriched pathways per sys-candidates"),
                         gsep_plus_containers %>% subset(gene_type!="cgc") %>% subset(sample%in%limited_sharing) %>% select(symbol, container) %>% unique %>% group_by(symbol) %>% summarise(value=n()) %>% mutate(sample_type="limited-sharing", info_type="Enriched pathways per sys-candidates"),
                         gsep_plus_containers %>% subset(gene_type!="cgc") %>% subset(sample%in%no_sharing) %>% select(symbol, container) %>% unique %>% group_by(symbol) %>% summarise(value=n()) %>% mutate(sample_type="no-sharing", info_type="Enriched pathways per sys-candidates"))
    ggplot(genes2info, aes(x=sample_type, y=value, fill=sample_type)) + geom_boxplot() + facet_wrap(~info_type, ncol = 1, scales = "free")


    ## Percentage with TP53 in each group
    gsep_plus_containers %>% subset(sample%in%high_sharing & symbol=="TP53") %>% select(sample) %>% unique %>% nrow
    gsep_plus_containers %>% subset(sample%in%intermediate_sharing & symbol=="TP53") %>% select(sample) %>% unique %>% nrow
    gsep_plus_containers %>% subset(sample%in%limited_sharing & symbol=="TP53") %>% select(sample) %>% unique %>% nrow
    gsep_plus_containers %>% subset(sample%in%no_sharing & symbol=="TP53") %>% select(sample) %>% unique %>% nrow

    ## Number of processes that genes are mapped
    gsep_plus_containers %>% subset(sample%in%high_sharing) %>% select(symbol, container) %>% unique %>% group_by(symbol) %>% summarise(n=n()) %>% .$n %>% summary()
    gsep_plus_containers %>% subset(sample%in%intermediate_sharing) %>% select(symbol, container) %>% unique %>% group_by(symbol) %>% summarise(n=n()) %>% nrow
    gsep_plus_containers %>% subset(sample%in%limited_sharing) %>% select(symbol, container) %>% unique %>% group_by(symbol) %>% summarise(n=n()) %>% nrow
    gsep_plus_containers %>% subset(sample%in%no_sharing) %>% select(symbol, container) %>% unique %>% group_by(symbol) %>% summarise(n=n()) %>% nrow

    high_sharing_processes = gsep_plus_containers %>% subset(sample%in%high_sharing) %>% select(sample, container) %>% unique %>% count(container) %>% arrange(desc(n)) %>% rename(high_sharing=n)
    intermediate_sharing_processes = gsep_plus_containers %>% subset(sample%in%intermediate_sharing) %>% select(sample, container) %>% unique %>% count(container) %>% arrange(desc(n)) %>% rename(intermediate_sharing=n)
    groups2processes = high_sharing_processes %>% left_join(intermediate_sharing_processes)
    limited_sharing_processes = gsep_plus_containers %>% subset(sample%in%limited_sharing) %>% select(sample, container) %>% unique %>% count(container) %>% arrange(desc(n)) %>% rename(limited_sharing=n)
    groups2processes = groups2processes %>% left_join(limited_sharing_processes)
    no_sharing_processes = gsep_plus_containers %>% subset(sample%in%no_sharing) %>% select(sample, container) %>% unique %>% count(container) %>% arrange(desc(n)) %>% rename(no_sharing=n)
    groups2processes = groups2processes %>% left_join(no_sharing_processes)
    write.table(groups2processes, file="Desktop/groups2processes.tsv", row.names = F, sep = "\t", quote = F)

    ## Get percentages
    groups2processes_perc = groups2processes %>% mutate(high_sharing=high_sharing/43, intermediate_sharing=intermediate_sharing/19, limited_sharing=limited_sharing/33, no_sharing=no_sharing/25)
    rnames = groups2processes_perc$container
    groups2processes_perc = groups2processes_perc[,c(2:5)] %>% data.frame()
    rownames(groups2processes_perc) = rnames
    heatmap.2(as.matrix(groups2processes_perc), col="greenred",cexCol = 1, trace = "none",margins = c(10, 35), keysize = 1, main = "Fraction of samples with each \n container alter")

    ## Network visualisation
    netwrk =graph_from_data_frame(net%>%subset(common.containers>0)%>%select(sample1, sample2),directed=F)
    E(netwrk)$weigth = net%>%subset(common.containers>0) %>% .$common.containers
    E(netwrk)$width = net%>%subset(common.containers>0) %>% .$common.containers
    plot(netwrk, vertex.size=2, edge.width=E(netwrk)$weight)


    ggnet2(netwrk, edge.size = toPlot$common.containers)
    E(netwrk)$weight = net$common.containers
    plot(netwrk, edge.width=net$common.containers, vertex.size=2)
    rdp=RedPort()
    calld(rdp)
    addGraph(rdp, netwrk)
    relax(rdp)

}

getVennOfPathways = function(){
    ## Interection of group-specific pathways and the enriched pathways of all the samples together

    ## Top 10 only
    top5 = read.table("Desktop/enriched_paths2samples_CGC_plus_top5syscans_noGroups_fdr_0.01_pathwayFilters_corrected.tsv", quote = "", sep="\t", header = T)
    t1 = top5 %>% select(pathway) %>% unique %>% mutate(inCGCplusTop5=1)
    top10 = read.table("Desktop/enriched_paths2samples_CGC_plus_top10syscans_noGroups_fdr_0.01_pathwayFilters_corrected.tsv", quote = "", sep="\t", header = T)
    t2 = top10 %>% select(pathway) %>% unique %>% mutate(inCGCplusTop10=1)
    top15 = read.table("Desktop/enriched_paths2samples_CGC_plus_top15syscans_noGroups_fdr_0.01_pathwayFilters_corrected.tsv", quote = "", sep="\t", header = T)
    t3 = top15 %>% select(pathway) %>% unique %>% mutate(inCGCplusTop15=1)
    top25 = read.table("Desktop/enriched_paths2samples_CGC_plus_top25syscans_noGroups_fdr_0.01_pathwayFilters_corrected.tsv", quote = "", sep="\t", header = T)
    t4 = top25 %>% select(pathway) %>% unique %>% mutate(inCGCplusTop25=1)
    top50 = read.table("Desktop/enriched_paths2samples_CGC_plus_top50syscans_noGroups_fdr_0.01_pathwayFilters_corrected.tsv", quote = "", sep="\t", header = T)
    t5 = top50 %>% select(pathway) %>% unique %>% mutate(inCGCplusTop50=1)

    t = t1 %>% full_join(t2) %>% full_join(t3) %>% full_join(t4) %>% full_join(t5)
    t = t %>% left_join(procs)

    write.xlsx(t, file="Desktop/Pathway_intersection_CGCPlusTopX.xlsx", row.names = F)


    top = read.table("Desktop/enriched_paths2samples_CGC_plus_top10syscans_noGroups_fdr_0.01_pathwayFilters_corrected.tsv", quote = "", sep="\t", header = T)
    t1 = top10 %>% select(pathway) %>% unique %>% mutate(inCGCplusTop10=1)
    top10 = read.table("Desktop/enriched_paths2samples_CGC_plus_top10syscans_noGroups_fdr_0.01_pathwayFilters_corrected.tsv", quote = "", sep="\t", header = T)
    t1 = top10 %>% select(pathway) %>% unique %>% mutate(inCGCplusTop10=1)


    onlytop10 = onlytop10 %>% left_join(paths2s2g2)

    paths_rec2 = data.frame(pathway=unique(onlytop10$pathway), inOnlyTop10=1)

    ## CGC plus top10
    plusCGC = read.table("Desktop/enriched_paths2samples_CGC_plus_top10syscans_noGroups_fdr_0.01.txt", quote = "", sep="\t", header = T) %>% select(-kernels_predicted, -kernels_predicted_no)
    paths2s2g3 = plusCGC %>% group_by(pathway) %>% summarise(samplesCGCplusTop10=length(unique(sample)), genesCGCplusTop10=length(unique(entrez))) %>% ungroup()
    paths2s2g3 = paths2s2g3 %>% subset(samplesCGCplusTop10>11)
    plusCGC = plusCGC %>% subset(pathway%in%paths2s2g3$pathway)
    plusCGC = plusCGC %>% left_join(paths2s2g3)

    paths_rec3 = data.frame(pathway=unique(plusCGC$pathway), inCGCplusTop10=1)

    ## Top 10 including CGC
    withCGC = read.table("Desktop/enriched_paths2samples_top10syscans_withCGC_noGroups_fdr_0.01.txt", quote = "", sep="\t", header = T) %>% select(-kernels_predicted, -kernels_predicted_no)
    paths2s2g4 = withCGC %>% group_by(pathway) %>% summarise(samplesTop10withCGC=length(unique(sample)), genesTop10withCGC=length(unique(entrez))) %>% ungroup()
    paths2s2g4 = paths2s2g4 %>% subset(samplesTop10withCGC>11)
    withCGC = withCGC %>% subset(pathway%in%paths2s2g4$pathway)
    withCGC = withCGC %>% left_join(paths2s2g4)

    paths_rec4 = data.frame(pathway=unique(withCGC$pathway), inTop10withCGC=1)


    ## Secrier 31 drivers
    # secrier31 = read.table("Desktop/enriched_paths2samples_noGroups_Secrier31Drivers_fdr_0.01.tsv", quote = "", sep="\t", header = T) %>% select(-kernels_predicted, -kernels_predicted_no)
    # paths2s2g5 = secrier31 %>% group_by(pathway) %>% summarise(samplesSecrier31=length(unique(sample)), genesSecrier31=length(unique(entrez))) %>% ungroup()
    # paths2s2g5 = paths2s2g5 %>% subset(samplesSecrier31>11)
    # secrier31 = secrier31 %>% subset(pathway%in%paths2s2g5$pathway)
    # secrier31 = secrier31 %>% left_join(paths2s2g5)
    #
    # paths_rec5 = data.frame(pathway=unique(secrier31$pathway), inSecrier31=1)

    ## All genes
    all_genes = read.table("Desktop/enriched_paths2samples_noGroups_all_genes_plus_CGC_withoutAmplifications_fdr_0.01.tsv", quote = "", sep="\t", header = T) %>% select(-kernels_predicted, -kernels_predicted_no)
    paths2s2g5 = all_genes %>% group_by(pathway) %>% summarise(samplesAllgenes=length(unique(sample)), genesAllgenes=length(unique(entrez))) %>% ungroup()
    paths2s2g5 = paths2s2g5 %>% subset(samplesAllgenes>11)
    all_genes = all_genes %>% subset(pathway%in%paths2s2g5$pathway)
    all_genes = all_genes %>% left_join(paths2s2g5)

    paths_rec5 = data.frame(pathway=unique(all_genes$pathway), inAllgenes31=1)

    paths_rec = paths_rec %>% full_join(paths_rec2) %>% full_join(paths_rec3) %>% full_join(paths_rec4) %>% full_join(paths_rec5)
    ## Get frist level 1 and level 2 pathways from config
    paths_rec = paths_rec %>% mutate(Level_1=ifelse(pathway%in%level_1_nodes, TRUE, FALSE), Level_2=ifelse(pathway%in%level_2_nodes, TRUE, FALSE))

    paths_rec = paths_rec %>% left_join(paths2s2g1) %>% left_join(paths2s2g2) %>% left_join(paths2s2g3) %>% left_join(paths2s2g4) %>% left_join(paths2s2g5)

    write.table(paths_rec, file="Desktop/paths_rec.tsv", sep="\t", row.names = F, quote=F)
    write.table(onlyCGC, file="Desktop/onlyCGC.tsv", row.names = F, quote = F, sep="\t")
    write.table(onlytop10, file="Desktop/onlytop10.tsv", row.names = F, quote = F, sep="\t")
    write.table(plusCGC, file="Desktop/plusCGC.tsv", row.names = F, quote = F, sep="\t")
    write.table(withCGC, file="Desktop/withCGC.tsv", row.names = F, quote = F, sep="\t")
    write.table(secrier31, file="Desktop/secrier31.tsv", row.names = F, quote = F, sep="\t")
    ## Combine all 3 in one table
    t = onlyCGC %>% full_join(plusCGC) %>% full_join(onlytop10) %>% full_join(withCGC) %>% left_join(paths_rec)



    library(gplots)

    venn(list(onlyCGC=unique(onlyCGC%>% .$pathway),
              onlytop10=unique(onlytop10%>% .$pathway),
              top10plusCGC=unique(plusCGC%>% .$pathway),
              top10withCGC=unique(withCGC%>% .$pathway),
              Allgenes=unique(all_genes%>% .$pathway)))

    venn(list(onlyCGC=unique(onlyCGC%>%subset(!pathway%in%c(level_1_nodes, level_2_nodes))%>% .$pathway),
              onlytop10=unique(onlytop10%>%subset(!pathway%in%c(level_1_nodes, level_2_nodes))%>% .$pathway),
              top10plusCGC=unique(plusCGC%>%subset(!pathway%in%c(level_1_nodes, level_2_nodes))%>% .$pathway),
              top10withCGC=unique(withCGC%>%subset(!pathway%in%c(level_1_nodes, level_2_nodes))%>% .$pathway),
              Allgenes=unique(all_genes%>%subset(!pathway%in%c(level_1_nodes, level_2_nodes))%>% .$pathway)))



    ## Fixing the table (pathway filters/pathway groups/ranks of CGCs to NA etc)
    d = read.table("Desktop/enriched_paths2samples_CGC_plus_top5syscans_noGroups_fdr_0.01.txt", quote = "", sep="\t", header = T) %>% select(-kernels_predicted, -kernels_predicted_no)
    d1 = d %>% subset(gene_type=="cgc") %>% select(-patient.rank) %>% mutate(patient.rank=NA)
    d2 = d %>% subset(gene_type!="cgc") %>% select(-patient.rank)
    genes2rank = syscan %>% select(sample, symbol, score) %>% arrange(sample, score) %>% group_by(sample) %>% mutate(patient.rank=row_number(-score)) %>% ungroup %>% select(-score) %>% data.frame()
    d2 = d2 %>% left_join(genes2rank)

    paths2s2g = d %>% group_by(pathway) %>% summarise(samples=length(unique(sample)), genes=length(unique(entrez))) %>% ungroup()
    paths2s2g = paths2s2g %>% subset(samples>11)
    paths2s2g = paths2s2g %>% mutate(Level_1=ifelse(pathway%in%level_1_nodes, TRUE, FALSE), Level_2=ifelse(pathway%in%level_2_nodes, TRUE, FALSE))
    paths2s2g = paths2s2g %>% subset(Level_1==FALSE & Level_2==FALSE)
    d = rbind(d1, d2) %>% subset(pathway%in%paths2s2g$pathway) %>% left_join(paths2s2g)

    ## Add pathway groups
    paths2func = read.xlsx("Desktop/Pathways2Processes_CGC_AND_TOP10_2203.xlsx", 1) %>% rename(pathway=Pathways.enriched.in.CGC_Plus_TopTen)
    top10 = top10 %>% left_join(paths2func)
    top5 = top5 %>% left_join(paths2func)
    top15 = top15 %>% left_join(paths2func)


    write.table(d, file="Desktop/enriched_paths2samples_CGC_plus_top10syscans_noGroups_fdr_0.01_pathwayFilters_corrected.tsv", row.names = F, sep="\t", quote = F)

    t1 = top5 %>% select(pathway) %>% unique %>% mutate(in_Top5=1)
    t2 = top10 %>% select(pathway) %>% unique %>% mutate(in_Top10=1)
    t3 = top15 %>% select(pathway) %>% unique %>% mutate(in_Top15=1)
    t = t1 %>% full_join(t2) %>% full_join(t3)
    write.xlsx(t, file="Desktop/Pathway_intersection_CGCPlusTopX.xlsx", row.names = F)




}

## For the pathways that are the same among the three groups
## Map the genes altered onto the pathways
annotatePathways = function(paths=read.delim(pipe("pbpaste"), header=F), group="c2a", save.dir="~/Desktop/results", length.cutoff=500){
    ## pathway visualisation
    require(ReactomePA)

    ## Get the pathways
    if(group=="mutagenic"){
        oneList = read.delim("Desktop/oneList_mutagenic.txt", header = T, sep = "\t") %>% subset(!is.na(pathway))
        genes = oneList %>% dplyr::select(entrez, score) %>% group_by(entrez) %>% summarise(score=mean(score)) %>% ungroup()
        g = genes$score
        names(g) = genes$entrez
    }else if(group=="c2a"){
        oneList = read.delim("Desktop/oneList_c2adom.txt", header = T, sep = "\t") %>% subset(!is.na(pathway))
        genes = oneList %>% dplyr::select(entrez, score) %>% group_by(entrez) %>% summarise(score=mean(score)) %>% ungroup()
        g = genes$score
        names(g) = genes$entrez
    }else if(group=="ddr"){
        oneList = read.delim("Desktop/oneList_ddr.txt", header = T, sep = "\t") %>% subset(!is.na(pathway))
    }

    ## Create the output directory
    dir.create(save.dir, showWarnings = FALSE)

    for (p in paths$V1){
        pdf(paste0(save.dir, "/", gsub("\\/", "_", p), "_", group, ".pdf"))
        genes = oneList %>% subset(pathway==p) %>% dplyr::select(entrez, score) %>% group_by(entrez) %>% summarise(score=mean(score)) %>% ungroup()
        g = genes$score
        names(g) = genes$entrez
        viewPathway(p, readable=TRUE, main=paste0(p, "(", length(g), ")"), foldChange = g, vertex.label.font = 0.5, vertex.label.cex=0.4)
        dev.off()
    }
}

## I just dumped code here that I used to analyse pathways/mappings to samples etc
pathwayMetaAnalysis = function(){

    ## Check union of samples and genes in the enriched pathways
    setwd("~/Desktop/")

    oneList = read.delim("oneList_ddr.txt", header = T, sep = "\t")
    top10 = read.delim("top10_ddr.txt", header = T, sep = "\t") %>% subset(!is.na(pathway))

    ## Read in the pathways of the cluster - usually imported with copy from excel
    paths = read.delim(pipe("pbpaste"), header=T)
    paths
    ## Union of samples
    s1 = oneList %>% subset(pathway%in%paths$V1) %>% select(sample) %>% unique %>% .$sample
    s2 = top10 %>% subset(pathway%in%paths$V1) %>% select(sample) %>% unique %>% .$sample
    length(unique(c(s1, s2)))
    ## union of genes
    g1 = oneList %>% subset(pathway%in%paths$V1) %>% select(entrez) %>% unique %>% .$entrez
    oneList %>% subset(pathway%in%paths$V1) %>% select(symbol) %>% unique %>% nrow

    g2 = top10 %>% subset(pathway%in%paths$V1) %>% select(entrez) %>% unique %>% .$entrez
    top10 %>% subset(pathway%in%paths$V1) %>% select(symbol) %>% unique %>% nrow
    length(unique(c(g1, g2)))



    ## Get a venn-like graph for the pathway clusters
    d = read.delim("Desktop/pathway_clustering_c2a.txt", sep = "\t", header = T) %>% subset(pathway!="")
    oneList = read.delim("Desktop/oneList_c2adom.txt", header = T, sep = "\t")
    top10 = read.delim("Desktop/top10_c2adom.txt", header = T, sep = "\t") %>% subset(!is.na(pathway))

    ## Get a sample-based list with all the pathways
    s1 = oneList %>% subset(pathway%in%d$pathway) %>% dplyr::select(sample, pathway) %>% left_join(d%>%dplyr::select(pathway, Group.name))
    s2 = top10 %>% subset(pathway%in%d$pathway) %>% dplyr::select(sample, pathway)  %>% left_join(d%>%dplyr::select(pathway, Group.name))

    ## Before we collapse and count the combinations
    # clusters = c("MB", "GE", "ML", "CA")
    # samples = sort(unique(s$sample))
    # samples_covered = NULL
    # for (cl in clusters){
    #     cat(cl, "\n")
    #     add = s %>% subset(Group.name==cl) %>% select(sample) %>% unique %>% .$sample
    #     samples_covered = c(samples_covered, add)
    #     print(length(unique(samples_covered)))
    #     if(length(unique(samples_covered))==37){
    #         stop("reached maximum samples")
    #     }
    #
    # }

    s = rbind(s1, s2) %>% group_by(sample) %>% summarise(paths=paste(sort(unique(Group.name)), collapse=","))
    s$group.no = apply(s, 1, function(x) length(unlist(strsplit(x[2], split = ","))))
    s %>% count(paths) %>% data.frame()

    ## Create a list of all possible combinations of pathways
    i=1
    max_combs = 5
    l = list()
    while (i <= max_combs){
        cat(i, "\n")
        l[[i]] = split(combn(unique(d$Group.name), i), rep(1:ncol(combn(unique(d$Group.name), i)), each = nrow(combn(unique(d$Group.name), i))))
        i = i + 1
    }

    ## search for patterns and count how many time all patterns are found
    ## This is basically the intersection
    p = list()
    for (i in 1:length(l)){
        for(y in 1:length(l[[i]])){
            p[[paste(l[[i]][[y]], collapse = "&")]] = apply(s, 1, function(x) colSums(data.frame(sapply(l[[i]][[y]], grepl, x[2], simplify = T)))==length(l[[i]][[y]]))
        }
    }
    ## Keep only those with at least one TRUE
    p = p[unlist(lapply(p, function(x) length(x[x==TRUE])>=1))]
    p = lapply(p, function(x) length(x[x==TRUE]))
    p = lapply(seq_along(p), function(y,n,i){
        ldply(y[[i]]) %>% mutate(groups=n[[i]])
    }, y=p, n=names(p))
    p = do.call(rbind, p)

    ## lets do the union
    pu = list()
    for (i in 1:length(l)){
        for(y in 1:length(l[[i]])){
            pu[[paste(l[[i]][[y]], collapse = "&")]] = apply(s, 1, function(x) colSums(data.frame(sapply(l[[i]][[y]][l[[i]][[y]]!=""], grepl, x[2], simplify = T)))>=1)
        }
    }
    ## Keep only those with at least one TRUE
    pu = pu[unlist(lapply(pu, function(x) length(x[x==TRUE])>=1))]
    pu = lapply(pu, function(x) length(x[x==TRUE]))
    pu = lapply(seq_along(pu), function(y,n,i){
        ldply(y[[i]]) %>% mutate(groups=n[[i]])
    }, y=pu, n=names(pu))
    pu = do.call(rbind, pu)

    colnames(p) = c("intersection", "groups")
    colnames(pu) = c("union", "groups")
    pa = p %>% full_join(pu)
    pa$groups.no = apply(pa, 1, function(x) length(unlist(strsplit(x[2], split = "&"))))
    write.table(pa%>%dplyr::select(intersection, union, groups, groups.no), file="Desktop/pathway_cluster_combinations_c2adom.tsv", quote = F, row.names = F, sep = "\t")
}


## Check in which pathways and samples KAT2A, KAT2B and related genes are
analyzeKATs = function(geneInfo_fn="/Users/fc-8s-imac/athena/data/geneInfoNCG5.Rdata",
                       gene_sets_dir="/Users/fc-8s-imac/athena/data/geneSets/reactome_12_12_16/",
                       gene_sets_gene_start=4,
                       gene_sets_identifiers = "reactome",
                       sub_identifier="symbol",
                       scores_fn="/Users/fc-8s-imac/athena/data/OAC/129_OAC/ML/OAC/syscans_4_with_gains_in_sys_cans.Rdata",
                       cohort_fn="/Users/fc-8s-imac/athena/data/OAC/129_OAC/Rdata/mainCohort.Rdata",
                       remove.amplifications=T,
                       remove.cgc=T){

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
    if(remove.cgc){
        syscan = syscan %>% subset(gene_type!="cgc")
    }
    if(remove.amplifications){
        syscan = syscan %>% subset(CNVGain!=1)
    }

    ## Get the genes for each group
    caDominant = syscan %>% subset(Group=="C>A/T dominant") %>% subset(gene_type!="CGC_discarded") %>%
        subset(type=="P") %>% select(entrez, symbol, score) %>% group_by(entrez, symbol) %>% summarise(score=mean(score)) %>% ungroup()
    addition = syscan %>% subset(Group=="C>A/T dominant") %>% subset(gene_type!="CGC_discarded") %>%
        subset(type=="C") %>% select(entrez, symbol) %>% unique %>% mutate(score=max(caDominant$score)+1)
    caDominant = rbind(caDominant, addition)
    caDominant = caDominant %>% arrange(desc(score)) %>% mutate(rank=row_number())
    if(sub_identifier=="entrez"){
        caDominant = caDominant %>% select(entrez, score, rank) %>% subset(!is.na(entrez))
    }else if(sub_identifier=="entrez"){
        caDominant = caDominant %>% select(symbol, score, rank) %>% subset(!is.na(symbol))
    }

    ddr = syscan %>% subset(Group=="DDR impaired") %>% subset(gene_type!="CGC_discarded") %>%
        subset(type=="P") %>% select(entrez, symbol, score) %>% group_by(entrez, symbol) %>% summarise(score=mean(score)) %>% ungroup()
    addition = syscan %>% subset(Group=="DDR impaired") %>% subset(gene_type!="CGC_discarded") %>%
        subset(type=="C") %>% select(entrez, symbol) %>% unique %>% mutate(score=max(ddr$score)+1)
    ddr = rbind(ddr, addition)
    ddr = ddr %>% arrange(desc(score)) %>% mutate(rank=row_number())
    if(sub_identifier=="entrez"){
        ddr = ddr %>% select(entrez, score, rank) %>% subset(!is.na(entrez))
    }else if(sub_identifier=="entrez"){
        ddr = ddr %>% select(symbol, score, rank) %>% subset(!is.na(symbol))
    }

    mutagenic = syscan %>% subset(Group=="Mutagenic") %>% subset(gene_type!="CGC_discarded") %>%
        subset(type=="P") %>% select(entrez, symbol, score) %>% group_by(entrez, symbol) %>% summarise(score=mean(score)) %>% ungroup()
    addition = syscan %>% subset(Group=="Mutagenic") %>% subset(gene_type!="CGC_discarded") %>%
        subset(type=="C") %>% select(entrez, symbol) %>% unique %>% mutate(score=max(mutagenic$score)+1)
    mutagenic = rbind(mutagenic, addition)
    mutagenic = mutagenic %>% arrange(desc(score)) %>% mutate(rank=row_number())
    if(sub_identifier=="entrez"){
        mutagenic = mutagenic %>% select(entrez, score, rank) %>% subset(!is.na(entrez))
    }else if(sub_identifier=="entrez"){
        mutagenic = mutagenic %>% select(symbol, score, rank) %>% subset(!is.na(symbol))
    }

    gene_sets = getPathway(pathway = gene_sets_identifiers)
    genes = c("KAT2A", "KAT2B")

    mut_enriched_pathways = read.xlsx("Desktop/mutagenic_enriched_paths.xlsx", 1)
    ca_enriched_pathways = read.xlsx("Desktop/c2a2t_dominant_enriched_paths.xlsx", 1)
    ddr_enriched_pathways = read.xlsx("Desktop/ddr_impaired_enriched_paths.xlsx", 1)

    res = NULL
    for (g in genes){
        cat(g, "\n")
        gs = gene_sets[unlist(lapply(gene_sets, function(x) g%in%x))]
        gs_unlisted = unique(unlist(gs))
        ## Make the list dataframe
        gs = lapply(seq_along(gs), function(y, n, i){
            ldply(y[[i]]) %>% mutate(pathway=n[[i]])
        }, y=gs, n=names(gs))


        gs = do.call(rbind, gs) %>% rename(symbol=V1) %>% mutate(gene.tested=g)

        pats = syscan %>% select(sample) %>% unique %>% .$sample

        for (p in pats){
            cat(p, "\n")
            ## get patient's group
            pgroup = syscan %>% subset(sample==p) %>% select(Group) %>% unique %>% .$Group
            ## get patient's genes
            ## Check which genes are in the combined list because depending on
            ## whether CGCs/amplifications are included the genes will change
            if (pgroup=="Mutagenic" & !is.na(pgroup)){
                gp = syscan %>% subset(sample==p & entrez%in%mutagenic$entrez)
                gp = gp %>% left_join(mutagenic %>% select(-score) %>% rename(rank.in.group=rank))
                enriched_pathways = mut_enriched_pathways
            }else if(pgroup=="C>A/T dominant" & !is.na(pgroup)){
                gp = syscan %>% subset(sample==p & entrez%in%caDominant$entrez)
                gp = gp %>% left_join(caDominant %>% rename(rank.in.group=rank))
                enriched_pathways = ca_enriched_pathways
            }else if(pgroup=="DDR impaired" & !is.na(pgroup)){
                gp = syscan %>% subset(sample==p & entrez%in%ddr$entrez)
                gp = gp %>% left_join(ddr %>% rename(rank.in.group=rank))
                enriched_pathways = ddr_enriched_pathways
            }else if(is.na(pgroup)){
                next
            }

            ## Get ranks in the patient
            gp = gp %>% arrange(desc(score)) %>% mutate(rank=row_number())
            if(sub_identifier=="entrez"){
                gp = gp %>% subset(entrez %in% gs_unlisted) %>% select(sample, entrez, symbol, rank, rank.in.group) %>% rename(rank.in.patient=rank)
            }else if(sub_identifier=="symbol"){
                gp = gp %>% subset(symbol %in% gs_unlisted) %>% select(sample, entrez, symbol, rank, rank.in.group) %>% rename(rank.in.patient=rank)
            }

            ## Add group information
            gp = gp %>% mutate(group=pgroup)
            ## Now bring in information for the pathways
            gp = gp %>% left_join(gs)

            ## Check if the pathway is enriched
            gp = gp %>% mutate(enriched=ifelse(pathway%in%enriched_pathways$pathway, TRUE, FALSE))
            res = rbind(res, gp)
        }
    }
    return(res)

}
## I run this in athena
## Function to create the submission file for randomSamplingGSEA function
createSub = function(genes=804, iters=1000000, chunks=1000, range=NULL, ncore=1, save_dir="~/athena/data/OAC/randomise_GSEA_cmds.txt",
                     res_dir = "/mnt/lustre/users/k1469280/mourikisa/data/OAC/129_OAC/randomise_GSEA"){
    iters_per_fn = iters/chunks
    if(!is.null(range)){ ## Use range in case you want to rename the files cause you want to add more iterations
        if(length(range)==length(1:chunks)){
            chunks=range
        }else{
            stop("Vector length mismatch")
        }
    }else{
        chunks = 1:chunks
    }
    cmds = NULL
    for(i in chunks){
        sub <- paste0("qsub -N job", i," -pe threaded ",ncore," -q FCgroup.q@node042,FCgroup.q@node043 -l h_vmem=12G ./submitRandomiseGSEA.sh ", genes, " ", iters_per_fn, " ", res_dir, " ", i)
        cmds = rbind(cmds, sub)
    }

    write.table(cmds, sep="\t", quote=F, row.names=F, col.names=F, file=save_dir)
}


parseRandomisation = function(res_dir="/mnt/lustre/users/k1469280/mourikisa/data/OAC/129_OAC/randomise_GSEA/",
                              paths = "/home/mourikisa/enriched_paths2samples_CGC_plus_top10syscans_noGroups_fdr_0.01_pathwayFilters_corrected.tsv",
                              save_dir="/home/mourikisa"){
  fns = list.files(res_dir, full.names = TRUE)

  ## Get 239
  d = read.table(paths, header = T, sep="\t")
  d = d %>% select(pathway) %>% unique %>% .$pathway

  rgsea = NULL
  for(fn in fns){
    cat(fn, "\n")
    load(fn)
    res = res[["gsea"]]
    #res = res %>% subset(pathway%in%d)
    rgsea = rbind(rgsea, res)
    gc()
  }

  # pdf(paste0(save_dir, "/frd_dist.pdf"))
  # hist(rgsea$fdr, breaks=100, main = paste0(res%>%.$pathway%>%unique, " \nrandom sampling GSEA FDR"), xlab = "FDR")
  # legend(x="topright",legend=cbind(paste0(names(summary(rgsea$fdr)), ":", summary(rgsea$fdr))))
  # dev.off()

  save(rgsea, file=paste0(save_dir,"/rgsea.Rdata"))

}


## This function was written to a seprate file so it can be run in parallel in athena
## I kept it here for backup
randomSamplingGSEA = function(n=NULL, ## This is the number of genes I random sample
                              times=NULL,
                              geneInfo_fn="~/data/geneInfoNCG5.Rdata",
                              gene_sets_dir="~/data/geneSets/reactome_12_12_16/",
                              gene_sets_gene_start=4,
                              gene_sets_identifiers = "reactome",
                              sub_identifier="symbol",
                              scores_fn="~/data/OAC/129_OAC/ML/OAC/syscans_4_with_gains_in_sys_cans.Rdata",
                              cohort_fn="~/data/OAC/129_OAC/Rdata/mainCohort.Rdata",
                              remove.amplifications=T,
                              remove.cgc=T){

    ## Library cleaning
    detach("package:dplyr", unload = T)
    require(plyr)
    require(tidyr)
    require(dplyr)

    if(is.null(n)){
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

    ex = replicate(times, sample(genes_to_sample, n))
    ex = split(ex, rep(1:ncol(ex), each = nrow(ex)))
    ## Add CGCs in the random samples
    ex = lapply(ex, function(x) x=c(x, cgcs))

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

    return(list(gsea=rgsea, gene_samples=ex, length_cgcs=length(cgcs), length_sys_cans=length(genes_to_sample)))
}

## Use this function to determine enriched pathways given a gene list for OAC subgroups or the total cohort
getTopGenes = function(geneInfo_fn="~/Mountpoints/rosalind_lustre/mourikisa/data/geneInfoNCG5.Rdata",
                       scores_fn="~/Mountpoints/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/syscan_noNotExpressed.Rdata",
                       given_tb=NULL, ## If this is not NULL you do GSEA on all the genes in this table
                       remove.amplifications=F,
                       remove.cgc=T, ## redundant with CGC.plus
                       CGC.plus = F, ## Set it to False to do GSEA only on syscans
                       genes.per.patient=10 ## set it to NULL to do GSEA only on the CGCs
                    ){

    load(geneInfo_fn)
    if(is.null(given_tb)){
        load(scores_fn)
    }else if(!is.null(given_tb)){ ## I do that to give directly prediction or training set
        syscan = given_tb
        ## These tables are the raw input for SVM and they have no symbol
        syscan = syscan %>% tibble::rownames_to_column() %>% separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")
        syscan = syscan %>% mutate(entrez=as.numeric(entrez)) %>% left_join(geneInfo%>%select(entrez, symbol))
    }

    ## Filters
    if(remove.cgc){
        syscan = syscan %>% subset(gene_type!="cgc")
    }
    if(remove.amplifications){
        syscan = syscan %>% subset(CNVGain!=1)
    }


    genes=NULL
    for(s in unique(syscan$sample)){
        sd = syscan %>% subset(sample==s)
        if(is.null(given_tb)){
            if(CGC.plus & !is.null(genes.per.patient)){
                sd1 = sd %>% subset(gene_type=="cgc")
                sd2 = sd %>% subset(gene_type!="cgc") %>% arrange(desc(score)) %>% slice(1:genes.per.patient)
                sd = rbind(sd1, sd2)
            }else if(CGC.plus==F & !is.null(genes.per.patient)){
                sd = sd %>% subset(gene_type!="cgc") %>% arrange(desc(score)) %>% slice(1:genes.per.patient)
            }else if(is.null(genes.per.patient)){
                sd = sd %>% subset(gene_type=="cgc")
            }
        }else if(!is.null(given_tb)){
            sd = sd
        }
        genes = rbind(genes, sd)
    }

    ## Store this to return also the genes
    genes_to_return = genes

    if(is.null(given_tb)){
        ## Get distribution
        print("CGCs")
        print(genes %>% subset(gene_type=="cgc") %>% group_by(sample) %>% summarise(genes=length(unique(symbol))) %>% ungroup() %>%
                  full_join(syscan%>%select(sample)%>%unique) %>% mutate(genes=ifelse(is.na(genes), 0, genes)) %>% .$genes %>% summary())
        print(genes %>% subset(gene_type=="cgc") %>% select(symbol) %>% unique %>% nrow)
        print("Non CGCs")
        print(genes %>% subset(gene_type!="cgc") %>% group_by(sample) %>% summarise(genes=length(unique(symbol))) %>% ungroup() %>%
                  full_join(syscan%>%subset(gene_type!="cgc")%>%select(sample)%>%unique) %>% mutate(genes=ifelse(is.na(genes), 0, genes)) %>% .$genes %>% summary())
        print(genes %>% subset(gene_type!="cgc") %>% select(symbol) %>% unique %>% nrow)
        print("Samples")
        print(genes %>% select(sample) %>% unique %>% nrow)
    }

    genes = genes$symbol
    ## Make the vector unique because we may have the same gene in two samples
    cat(paste0("Genes extracted: ", length(genes)), "\n")
    genes = unique(genes)
    cat(paste0("Unique genes: ", length(genes)), "\n")

    return(genes_to_return)
}

## Use this function to determine enriched pathways given a gene list for OAC subgroups or the total cohort
pathwayEnrichmentTopGenes = function(geneInfo_fn="~/Mountpoints/rosalind_lustre/mourikisa/data/geneInfoNCG5.Rdata",
                                     gene_sets_dir="~/Mountpoints/rosalind_lustre/mourikisa/data/geneSets/reactome_12_12_16/",
                                     gene_sets_gene_start=4,
                                     gene_sets_identifiers = "reactome",
                                     sub_identifier="symbol",
                                     scores_fn="~/Mountpoints/rosalind_lustre/mourikisa/data/OAC/Combined_ICGC_cohort_129_66_71/sysSVM/OAC/syscan_noNotExpressed.Rdata",
                                     given_tb=NULL, ## If this is not NULL you do GSEA on all the genes in this table
                                     remove.amplifications=T,
                                     remove.cgc=F, ## redundant with CGC.plus
                                     CGC.plus = T, ## Set it to False to do GSEA only on syscans
                                     genes.per.patient=10, ## set it to NULL to do GSEA only on the CGCs
                                     by.groups=F,
                                     ranks=F,
                                     path_length=T,
                                     min_path_length=10,
                                     max_path_length=500,
                                     pathlvls_excl=T,
                                     identifier.file="~/Mountpoints/rosalind_lustre/mourikisa/data/geneSets/reactome_12_12_16/complete_list_of_pathways.tsv",
                                     hierarchy.file="~/Mountpoints/rosalind_lustre/mourikisa/data/geneSets/reactome_12_12_16/pathway_hierarchy.tsv"){

    ## pathway parameter in the following function can be kegg/reactome/all
    ## sub parameter can be entrez/symbol
    getPathway = function(pathway="kegg", sub=sub_identifier, genes_start=gene_sets_gene_start,
                          path_dir=gene_sets_dir,
                          len=path_length,
                          min_len=min_path_length,
                          max_len=max_path_length,
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


        if(len==TRUE & !is.null(min_len) & !is.null(max_len)){
            p = p[unlist(lapply(p, function(x) length(x)>=min_len))]
            p = p[unlist(lapply(p, function(x) length(x)<=max_len))]
        }

        if(pathlvls_excl){
            hier = read.table(hierarchy.file, sep = "\t", header = F)
            hier = read.table(hierarchy.file, sep = "\t", header = F)
            colnames(hier) = c("path1", "path2")
            ## Careful here cause some pathways have apostrophes in their names
            path_ids = read.table(identifier.file, sep = "\t", header = F, quote = "")
            colnames(path_ids) = c("path_id", "path_name", "species")
            path_ids = path_ids %>% subset(species=="Homo sapiens") %>% rename(pathway=path_name) %>% select(-species)

            ## Exclude leading nad trailing spaces from pathway names
            path_ids$pathway = gsub("^\\s+|\\s+$", "", path_ids$pathway)
            path_ids$path_id = gsub("^\\s+|\\s+$", "", path_ids$path_id)

            ## Add info to herarchy table and subset for the homo sapiens interactions
            hier = hier %>% subset(path1%in%path_ids$path_id | path2%in%path_ids$path_id) %>% left_join(path_ids%>%select(path_id, pathway)%>%rename(path1=path_id, pathway1=pathway)) %>%
                left_join(path_ids%>%select(path_id, pathway)%>%rename(path2=path_id, pathway2=pathway))

            ## Trials to exclude level 1&2 nodes
            level_1_nodes = c("Developmental Biology",
                              "Reproduction",
                              "Circadian Clock",
                              "Hemostasis",
                              "Neuronal System",
                              "Immune System",
                              "Signal Transduction",
                              "Disease",
                              "DNA Repair",
                              "Metabolism",
                              "Gene Expression",
                              "Chromatin organization",
                              "Transmembrane transport of small molecules",
                              "DNA replication",
                              "Metabolism of proteins",
                              "Muscle contraction",
                              "Cell Cycle",
                              "Organelle biogenesis and maintenance",
                              "Mitophagy",
                              "Vesicle-mediated transport",
                              "Programmed Cell Death",
                              "Cellular responses to stress",
                              "Extracellular matrix organization",
                              "Cell-Cell communication")

            level_2_nodes = hier %>% subset(pathway1%in%level_1_nodes | pathway2%in%level_1_nodes)
            level_2_nodes = unique(c(level_2_nodes$pathway1, level_2_nodes$pathway2))
            level_2_nodes = level_2_nodes[!level_2_nodes%in%level_1_nodes]
            p = p[!names(p)%in%level_1_nodes]
            p = p[!names(p)%in%level_2_nodes]

        }

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

    ## x is a vector of genes
    ## y is a list of pathways
    ## identifier can be entrez/symbol
    updatePathway = function(x, y, identifier=sub_identifier, min_len=min_path_length, max_len=max_path_length, len=path_length,
                             gi_fn=geneInfo_fn){

        ## This function removes all genes that are not in the gene list of interest
        ## from the pathways

        if (identifier=="entrez"){
            y = lapply(y, as.numeric)
        }else if (identifier=="symbol"){
            y = lapply(y, as.character)
        }

        ## Update for the list of genes that are in our gene set
        y = lapply(y, function(g) g[g%in%x])
        ## Remove empty sets
        y = y[unlist(lapply(y, function(x) length(x)>0))]

        if(len==TRUE){
            y = y[unlist(lapply(y, function(x) length(x)>min_len))]
            y = y[unlist(lapply(y, function(x) length(x)<max_len))]
        }

        message(paste0("Number of pathways: ", length(y)))
        if(identifier=="entrez"){
            genes_p = as.numeric(unique(unlist(y)))
        }else if(identifier=="symbol"){
            genes_p = unique(unlist(y))
        }
        message(paste0("Number of genes: ", length(genes_p)))
        ## Get gene info to check overlap with 19014
        load(gi_fn)
        geneInfo = geneInfo %>% subset(duplicability!=-1) %>% select(symbol, entrez)
        if(identifier=="entrez"){
            message(paste0("Number of genes (19,014): ", length(genes_p[genes_p%in%geneInfo$entrez])))
        }else if(identifier=="symbol"){
            message(paste0("Number of genes (19,014): ", length(genes_p[genes_p%in%geneInfo$symbol])))
        }


        return(list(gene_sets=y, genes_in_p=length(genes_p)))
    }

    load(geneInfo_fn)
    if(is.null(given_tb)){
        load(scores_fn)
    }else if(!is.null(given_tb)){ ## I do that to give directly prediction or training set
        syscan = given_tb
        ## These tables are the raw input for SVM and they have no symbol
        syscan = syscan %>% tibble::rownames_to_column() %>% separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")
        syscan = syscan %>% mutate(entrez=as.numeric(entrez)) %>% left_join(geneInfo%>%select(entrez, symbol))
    }

    ## Filters
    if(remove.cgc){
        syscan = syscan %>% subset(gene_type!="cgc")
    }
    if(remove.amplifications){
        syscan = syscan %>% subset(CNVGain!=1)
    }

    ## Exclude samples that do not correspond to
    #syscan = syscan %>% subset(!is.na(Group))

    gene_sets = getPathway(pathway = gene_sets_identifiers)
    genesets2length = lapply(gene_sets, function(x) length(x))
    genesets2length = data.frame(pathway=names(genesets2length), pathway.length=unlist(genesets2length))

    ## Checkpoint here whether the analysis should consider OAC subgroups or not
    if(ranks==F){
        if(by.groups){
            group_enrichment = NULL
            for(g in unique(syscan$Group[!is.na(syscan$Group)])){
                cat(g, "\n")
                gd = syscan %>% subset(Group==g)

                group_genes = NULL
                ## Get the top genes in each sample
                for(s in unique(gd$sample)){
                    sd = gd %>% subset(sample==s)
                    if(CGC.plus){
                        sd1 = sd %>% subset(gene_type=="cgc")
                        sd2 = sd %>% subset(gene_type!="cgc") %>% arrange(desc(score)) %>% slice(1:genes.per.patient)
                        sd = rbind(sd1, sd2)
                    }else{
                        sd = sd %>% arrange(desc(score)) %>% slice(1:genes.per.patient)
                    }
                    group_genes = c(group_genes, sd$symbol)
                }

                ## Make the vector unique because we may have the same gene in two samples
                group_genes = unique(group_genes)

                ## x is the gene set and y your gene list
                run.hypergeometric = function(x, y, universe=19014){
                    P = phyper(length(intersect(x,y)), length(x), universe-length(x), length(y), lower.tail = F)
                    p = data.frame(p.value=P, gene.set.length=length(x), genes.in.gene.set=length(intersect(x,y)))
                    return(p)
                }

                l = lapply(gene_sets, function(x) run.hypergeometric(x, group_genes))
                l = lapply(seq_along(l), function(y, n, i){
                    ldply(y[[i]]) %>% mutate(pathway=n[[i]])
                }, y=l, n=names(l))
                l = do.call(rbind, l) %>% rename(value=V1, type=.id)
                l = l %>% unique %>% spread(type, value)

                l$fdr = p.adjust(l$p.value, method="BH")
                l$number.of.genes.tested = length(group_genes)
                l$group = g
                group_enrichment = rbind(group_enrichment, l)
                return(group_enrichment)
            }
        }else{ ## Do GSEA without considering subgroups
            genes=NULL
            for(s in unique(syscan$sample)){
                sd = syscan %>% subset(sample==s)
                if(is.null(given_tb)){
                    if(CGC.plus & !is.null(genes.per.patient)){
                        sd1 = sd %>% subset(gene_type=="cgc")
                        sd2 = sd %>% subset(gene_type!="cgc") %>% arrange(desc(score)) %>% slice(1:genes.per.patient)
                        sd = rbind(sd1, sd2)
                    }else if(CGC.plus==F & !is.null(genes.per.patient)){
                        sd = sd %>% subset(gene_type!="cgc") %>% arrange(desc(score)) %>% slice(1:genes.per.patient)
                    }else if(is.null(genes.per.patient)){
                        sd = sd %>% subset(gene_type=="cgc")
                    }
                }else if(!is.null(given_tb)){
                    sd = sd
                }
                genes = rbind(genes, sd)
            }

            ## Store this to return also the genes
            genes_to_return = genes

            if(is.null(given_tb)){
                ## Get distribution
                print("CGCs")
                print(genes %>% subset(gene_type=="cgc") %>% group_by(sample) %>% summarise(genes=length(unique(symbol))) %>% ungroup() %>%
                          full_join(syscan%>%select(sample)%>%unique) %>% mutate(genes=ifelse(is.na(genes), 0, genes)) %>% .$genes %>% summary())
                print(genes %>% subset(gene_type=="cgc") %>% select(symbol) %>% unique %>% nrow)
                print(genes %>% subset(gene_type=="cgc" & symbol%in%unique(unlist(gene_sets))) %>% select(symbol) %>% unique %>% nrow)
                print("Non CGCs")
                print(genes %>% subset(gene_type!="cgc") %>% group_by(sample) %>% summarise(genes=length(unique(symbol))) %>% ungroup() %>%
                          full_join(syscan%>%subset(gene_type!="cgc")%>%select(sample)%>%unique) %>% mutate(genes=ifelse(is.na(genes), 0, genes)) %>% .$genes %>% summary())
                print(genes %>% subset(gene_type!="cgc") %>% select(symbol) %>% unique %>% nrow)
                print(genes %>% subset(gene_type!="cgc" & symbol%in%unique(unlist(gene_sets))) %>% select(symbol) %>% unique %>% nrow)
                print("Samples")
                print(genes %>% select(sample) %>% unique %>% nrow)
            }

            genes = genes$symbol
            ## Make the vector unique because we may have the same gene in two samples
            cat(paste0("Genes extracted: ", length(genes)), "\n")
            genes = unique(genes)
            cat(paste0("Unique genes: ", length(genes)), "\n")

            ## Get how many of those are mapped to pathways
            genesInGeneSets = sum(genes%in%unique(unlist(gene_sets)))
            cat(paste0("Unique genes in gene sets: ", genesInGeneSets), "\n")

            ## x is the gene set and y your gene list
            run.hypergeometric = function(x, y, universe=19014){
                P = phyper(length(intersect(x,y)), length(x), universe-length(x), length(y), lower.tail = F)
                p = data.frame(p.value=P, gene.set.length=length(x), genes.in.gene.set=length(intersect(x,y)))
                return(p)
            }

            l = lapply(gene_sets, function(x) run.hypergeometric(x, genes))
            l = lapply(seq_along(l), function(y, n, i){
                ldply(y[[i]]) %>% mutate(pathway=n[[i]])
            }, y=l, n=names(l))
            l = do.call(rbind, l) %>% rename(value=V1, type=.id)
            l = l %>% unique %>% spread(type, value)

            l$fdr = p.adjust(l$p.value, method="BH")
            l$number.of.genes.tested = length(genes)

            ## Add here annotation for the pathway level in reactome
            hier = read.table(hierarchy.file, sep = "\t", header = F)
            hier = read.table(hierarchy.file, sep = "\t", header = F)
            colnames(hier) = c("path1", "path2")
            ## Careful here cause some pathways have apostrophes in their names
            path_ids = read.table(identifier.file, sep = "\t", header = F, quote = "")
            colnames(path_ids) = c("path_id", "path_name", "species")
            path_ids = path_ids %>% subset(species=="Homo sapiens") %>% rename(pathway=path_name) %>% select(-species)

            ## Exclude leading nad trailing spaces from pathway names
            path_ids$pathway = gsub("^\\s+|\\s+$", "", path_ids$pathway)
            path_ids$path_id = gsub("^\\s+|\\s+$", "", path_ids$path_id)

            ## Add info to herarchy table and subset for the homo sapiens interactions
            hier = hier %>% subset(path1%in%path_ids$path_id | path2%in%path_ids$path_id) %>% left_join(path_ids%>%select(path_id, pathway)%>%rename(path1=path_id, pathway1=pathway)) %>%
            left_join(path_ids%>%select(path_id, pathway)%>%rename(path2=path_id, pathway2=pathway))

            ## Trials to exclude level 1&2 nodes
            level_1_nodes = c("Developmental Biology",
                                  "Reproduction",
                                  "Circadian Clock",
                                  "Hemostasis",
                                  "Neuronal System",
                                  "Immune System",
                                  "Signal Transduction",
                                  "Disease",
                                  "DNA Repair",
                                  "Metabolism",
                                  "Gene Expression",
                                  "Chromatin organization",
                                  "Transmembrane transport of small molecules",
                                  "DNA replication",
                                  "Metabolism of proteins",
                                  "Muscle contraction",
                                  "Cell Cycle",
                                  "Organelle biogenesis and maintenance",
                                  "Mitophagy",
                                  "Vesicle-mediated transport",
                                  "Programmed Cell Death",
                                  "Cellular responses to stress",
                                  "Extracellular matrix organization",
                                  "Cell-Cell communication")

            level_2_nodes = hier %>% subset(pathway1%in%level_1_nodes | pathway2%in%level_1_nodes)
            level_2_nodes = unique(c(level_2_nodes$pathway1, level_2_nodes$pathway2))
            level_2_nodes = level_2_nodes[!level_2_nodes%in%level_1_nodes]
            l = l %>% mutate(Level1=ifelse(pathway%in%level_1_nodes, TRUE, FALSE),
                            Level2=ifelse(pathway%in%level_2_nodes, TRUE, FALSE))

            return(list(gsea=l, genes=genes_to_return, paths=gene_sets))
        }
    }else if(ranks==T){

        ## x is the pathway
        ## g is the data frame with ranks
        run.singed.ks = function(x, g, func_fn="/Users/fc-8s-imac/athena/data/OAC/signed_ks.R"){
            source(func_fn)
            if(sub_identifier=="entrez"){
                s = g %>% subset(entrez%in%x) %>% select(rank) %>% .$rank
            }else if(sub_identifier=="symbol"){
                s = g %>% subset(symbol%in%x) %>% select(rank) %>% .$rank
            }

            k <- ks.test.2(s, (1:max(g$rank))[-s], maxCombSize=10^10)
            res = list(es=k$ES, p.value=k$p, no.genes=length(x), ranks=paste(s, collapse=","))
            return(res)
        }

        cat("Rank mode", "\n")
        genes=NULL
        for(s in unique(syscan$sample)){
            sd = syscan %>% subset(sample==s)
            sd = sd %>% arrange(desc(score)) %>% slice(1:genes.per.patient)
            genes = rbind(genes, sd)
        }

        ## Create one list for all the genes
        genes = genes %>% select(symbol, score) %>% group_by(symbol) %>% summarise(score=mean(score)) %>% ungroup()
        genes = genes %>% arrange(desc(score)) %>% mutate(rank=row_number())

        ## Just in case there are duplicated pathways
        dups = names(gene_sets)[duplicated(names(gene_sets))]
        for (d in dups){
            idx = which(names(gene_sets)==d)
            gene_sets[idx[2:length(idx)]] <- NULL
        }

        l = lapply(gene_sets, function(x) length(x))
        l = lapply(seq_along(l), function(y, n, i){
            ldply(y[[i]]) %>% mutate(pathway=n[[i]])
        }, y=l, n=names(l))
        l = do.call(rbind, l) %>% rename(pathway_length=V1)

        if(sub_identifier=="entrez"){
            gene_sets_up = updatePathway(genes$entrez, gene_sets)
        }else if(sub_identifier=="symbol"){
            gene_sets_up = updatePathway(genes$symbol, gene_sets)
        }

        gene_sets_up = gene_sets_up[["gene_sets"]]

        res = lapply(gene_sets_up, function(x) run.singed.ks(x, genes))
        ## Make res data frame to work with
        res = lapply(seq_along(res), function(y, n, i){
            ldply(y[[i]]) %>% mutate(pathway=n[[i]])
        }, y=res, n=names(res))
        res = do.call(rbind, res) %>% spread(.id, V1)
        res$fdr = p.adjust(res$p.value, method="BH")
        res$bonf = p.adjust(res$p.value, method="bonferroni")
        res = res %>% left_join(l)
        return(res)
    }
}

## function to take the distributions of genes mapped etc after the sample-by-sample GSEA
checkSampleBySample = function(gsea, fdr.cutoff=0.1){

    gsea = gsea[["gsea"]]
    ## Number of genes mapped to the pathways
    all_mappings = NULL
    for (g in unique(gsea$group)){
        if(is.na(g)){
            next
        }
        cat(g, "\n")
        group_mappings = NULL
        d = gsea %>% subset(group==g & fdr<fdr.cutoff & es>0)
        for (s in unique(d$sample)){
            cat(s, "\n")
            dd = d %>% subset(sample==s)
            u_genes = dd %>% mutate(ranks=strsplit(ranks, ",")) %>% unnest()
            u_genes = unique(dd$ranks)
            x = data.frame(sample=s, enriched_pathways=nrow(dd), genes_in_pathways=length(u_genes))
            group_mappings = rbind(group_mappings, x)
        }
        group_mappings = group_mappings %>% mutate(group=g)
        all_mappings = rbind(all_mappings, group_mappings)
        rm(group_mappings)
    }


}

comparePathways = function(gsea, r=2, c=2){

    ## Venn diagrams
    par(mfrow=c(r,c))
    gs = c("kegg", "reactome", "all")
    for(g in gs){
        d = gsea %>% subset(gene_set==g)
        d1 = d%>%subset(fdr<0.01 & es>0 & group=="C>A/T dominant")%>%.$pathway
        d2 = d%>%subset(fdr<0.01 & es>0 & group=="DDR impaired")%>%.$pathway
        d3 = d%>%subset(fdr<0.01 & es>0 & group=="Mutagenic")%>%.$pathway
        if(g=="kegg"){
            venn(list("C>A/T dominant_kegg"=d1, "DDR impaired_kegg"=d2, "Mutagenic_kegg"=d3))
        }else if(g=="reactome"){
            venn(list("C>A/T dominant_reactome"=d1, "DDR impaired_reactome"=d2, "Mutagenic_reactome"=d3))
        }else if(g=="all"){
            venn(list("C>A/T dominant_all"=d1, "DDR impaired_all"=d2, "Mutagenic_all"=d3))
        }

    }
    par(mfrow=c(1,1))
}

pathways2hallmarks = function(gsea,
                              gene_sets_dir="/Volumes/mourikisa/data/geneSets/reactome_12_12_16/",
                              gene_sets_gene_start=4,
                              geneInfo_fn="/Volumes/mourikisa/data/geneInfoNCG5.Rdata",
                              hallmarks_fn="/Volumes/mourikisa/data/geneSets/Supplementary Table S6 - hallmark curation keywords.xlsx"){
    require(xlsx)
    getPathway = function(pathway="reactome", sub="symbol", genes_start=gene_sets_gene_start,
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

    ## Add pathway hallmarks annotation
    halls = read.xlsx(hallmarks_fn, 1) %>%
        rename(keyword=Keyword.searched.in.pathway.name) %>%
        rename(keyword.gene=Keywork.searched.in.the.names.of.the.composing.genes) %>%
        mutate(keyword=gsub("^\\s+|\\s+$", "", keyword)) %>%
        mutate(keyword.gene=gsub("^\\s+|\\s+$", "", keyword.gene))

    searchPathName = halls %>% select(Hallmark, keyword) %>%
        mutate(keyword=strsplit(keyword, ",")) %>%
        unnest() %>% mutate(keyword=gsub("^\\s+|\\s+$", "", keyword))

    searchGenes = halls %>% select(Hallmark, keyword.gene) %>%
        mutate(keyword.gene=strsplit(keyword.gene, ",")) %>%
        unnest() %>% mutate(keyword.gene=gsub("^\\s+|\\s+$", "", keyword.gene))
    searchGenes = searchGenes %>% subset(!is.na(keyword.gene))

    pl = getPathway()

    hm_path_names = NULL
    message("Searching path names")
    for(hm in unique(searchPathName$Hallmark)){
      toMatch = searchPathName %>% subset(Hallmark==hm) %>% .$keyword
      matches = unique (grep(paste(toMatch,collapse="|"),
                              unique(names(pl)), value=TRUE))
      hms = data.frame(pathway=matches) %>% mutate(Hallmark=hm, evidence="name")
      hm_path_names = rbind(hm_path_names, hms)
    }
    message("Searching genes in each path")
    pl_df = ldply(pl, data.frame)
    colnames(pl_df) = c("pathway", "gene")
    hm_path_genes = pl_df %>% left_join(searchGenes, by=c("gene"="keyword.gene")) %>%
      subset(!is.na(Hallmark)) %>% select(pathway, Hallmark) %>% unique %>%
      mutate(evidence="gene")

    p2h = rbind(hm_path_names, hm_path_genes)
    p2h = p2h %>% group_by(pathway) %>% summarise(Hallmark=paste(unique(Hallmark), collapse=","),
                                                  evidence=paste(unique(evidence), collapse=","))

    gsea = gsea %>% left_join(p2h)
    annotation = gsea %>% select(sample, group) %>% unique
    annotation$group[is.na(annotation$group)] = "No group"

    kat_pathways = unique(names(pl[unlist(lapply(pl, function(x) "KAT2A"%in%x))]), names(pl[unlist(lapply(pl, function(x) "KAT2B"%in%x))]))
    kat_pathways = data.frame(pathway=kat_pathways, KAT_pathway=TRUE)
    gsea = gsea %>% left_join(kat_pathways)

    hm = "KAT2A/B related pathways"


    message(hm)
    toPlot = gsea %>% subset(p.value<0.05) %>% subset(KAT_pathway==TRUE) %>% mutate(p.value=ifelse(p.value==0, 1e-10, p.value)) %>% mutate(path_ratio=as.numeric(no.genes)/as.numeric(pathway_length))
    #toPlot = gsea %>% subset(p.value<0.05) %>% subset(grepl(hm, Hallmark)) %>% mutate(p.value=ifelse(p.value==0, 1e-10, p.value))
    toPlot = toPlot %>% select(pathway, sample, group, p.value, fdr) %>% mutate(value= -log10(as.numeric(p.value)))
    toPlot = toPlot %>% select(pathway, sample, path_ratio) %>% spread(pathway, path_ratio)
    toPlot = toPlot %>% left_join(annotation)
    rownames(toPlot) = toPlot$sample
    toPlot = toPlot[2:length(toPlot)]
    toPlot[is.na(toPlot)] = 0



        h1 = Heatmap(toPlot%>%select(-group),
                     column_title = hm,
                     #cluster_rows = FALSE,
                     #cluster_columns = FALSE,
                     row_names_gp = gpar(fontsize = 6),
                     column_names_gp = gpar(fontsize = 8),
                     heatmap_legend_param = list(title = "-log10(p.value)",
                                                 title_gp = gpar(fontsize = 12, fontface = "bold"),
                                                 color_bar = "continuous",
                                                 grid_height = unit(6, "mm"),
                                                 grid_width = unit(6, "mm"),
                                                 labels_gp = gpar(fontsize = 12)))

        h2 = Heatmap(toPlot[,"group"], name = "group",
                     heatmap_legend_param = list(title = "group",
                                                 title_gp = gpar(fontsize = 12, fontface = "bold"),
                                                 color_bar = "discrete",
                                                 grid_height = unit(6, "mm"),
                                                 grid_width = unit(6, "mm"),
                                                 labels_gp = gpar(fontsize = 12)),
                     width = unit(6, "mm"))


        draw(h1 + h2)

}
## Check sample coverage for genes in the enriched pathways (combined list for each group)
enrichedPathways2Samples = function(gsea, fdr.cutoff=0.01,
                                    path="reactome",
                                    sub_identifier="symbol",
                                    gene_sets_gene_start=4,
                                    gene_sets_dir="/Users/fc-8s-imac/athena/data/geneSets/reactome_12_12_16/",
                                    geneInfo_fn="/Users/fc-8s-imac/athena/data/geneInfoNCG5.Rdata",
                                    scores_fn="/Users/fc-8s-imac/athena/data/OAC/129_OAC/ML/OAC/syscans_4_with_gains_in_sys_cans.Rdata",
                                    cohort_fn="/Users/fc-8s-imac/athena/data/OAC/129_OAC/Rdata/mainCohort.Rdata",
                                    remove.amplifications=T,
                                    remove.cgc=T,
                                    path_length=F,
                                    es.positive=T,
                                    min_path_length=10,
                                    max_path_length=500){

    ## get patient data
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
    if(remove.cgc){
        syscan = syscan %>% subset(gene_type!="cgc")
    }
    if(remove.amplifications){
        syscan = syscan %>% subset(CNVGain!=1)
    }

    gsea = gsea[["gsea"]]

    getPathway = function(pathway=path, sub=sub_identifier, genes_start=gene_sets_gene_start,
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

    ## x is a vector of genes
    ## y is a list of pathways
    ## identifier can be entrez/symbol
    updatePathway = function(x, y, identifier=sub_identifier, min_len=min_path_length, max_len=max_path_length, len=path_length,
                             exclude_pattern="CANCER|LEUKEMIA|MELANOMA|GLIOMA",
                             gi_fn=geneInfo_fn){

        ## This function removes all genes that are not in the gene list of interest
        ## from the pathways

        if (identifier=="entrez"){
            y = lapply(y, as.numeric)
        }else if (identifier=="symbol"){
            y = lapply(y, as.character)
        }

        ## Update for the list of genes that are in our gene set
        y = lapply(y, function(g) g[g%in%x])
        ## Remove empty sets
        y = y[unlist(lapply(y, function(x) length(x)>0))]

        if(len==TRUE){
            y = y[unlist(lapply(y, function(x) length(x)>min_len))]
            y = y[unlist(lapply(y, function(x) length(x)<max_len))]
        }
        y = y[names(y)[!grepl(exclude_pattern, names(y))]]

        message(paste0("Number of pathways: ", length(y)))
        if(identifier=="entrez"){
            genes_p = as.numeric(unique(unlist(y)))
        }else if(identifier=="symbol"){
            genes_p = unique(unlist(y))
        }
        message(paste0("Number of genes: ", length(genes_p)))
        ## Get gene info to check overlap with 19014
        load(gi_fn)
        geneInfo = geneInfo %>% subset(duplicability!=-1) %>% select(symbol, entrez)
        if(identifier=="entrez"){
            message(paste0("Number of genes (19,014): ", length(genes_p[genes_p%in%geneInfo$entrez])))
        }else if(identifier=="symbol"){
            message(paste0("Number of genes (19,014): ", length(genes_p[genes_p%in%geneInfo$symbol])))
        }


        return(list(gene_sets=y, genes_in_p=length(genes_p)))
    }

    reactome = getPathway()

    ## Select enriched pathways based on the fdr
    genes_in_enriched_pathways = NULL
    if(es.positive==F){
        for(g in unique(gsea$group)){
            genes = syscan %>% subset(Group==g) %>% subset(gene_type!="CGC_discarded") %>%
                subset(type=="P") %>% select(entrez, symbol, score) %>% group_by(entrez, symbol) %>% summarise(score=mean(score)) %>% ungroup()
            addition = syscan %>% subset(Group==g) %>% subset(gene_type!="CGC_discarded") %>%
                subset(type=="C") %>% select(entrez, symbol) %>% unique %>% mutate(score=max(genes$score)+1)
            genes = rbind(genes, addition)
            d = gsea %>% subset(group==g & fdr<fdr.cutoff)
            gs = updatePathway(genes$symbol, reactome)
            gs = gs[["gene_sets"]]
            gs = gs[names(gs)[names(gs)%in%d$pathway]]

            ## Unlist all pathways and take unique genes
            gs = unique(unlist(gs))

            dd = syscan %>% subset(symbol%in% gs & Group==g) %>% count(sample) %>% mutate(group=g, genes_in_enriched_pathways=length(gs))
            genes_in_enriched_pathways = rbind(genes_in_enriched_pathways, dd)
        }
    }else{
        for(g in unique(gsea$group)){
            genes = syscan %>% subset(Group==g) %>% subset(gene_type!="CGC_discarded") %>%
                subset(type=="P") %>% select(entrez, symbol, score) %>% group_by(entrez, symbol) %>% summarise(score=mean(score)) %>% ungroup()
            addition = syscan %>% subset(Group==g) %>% subset(gene_type!="CGC_discarded") %>%
                subset(type=="C") %>% select(entrez, symbol) %>% unique %>% mutate(score=max(genes$score)+1)
            genes = rbind(genes, addition)
            d = gsea %>% subset(group==g & fdr<fdr.cutoff & es>0)
            gs = updatePathway(genes$symbol, reactome)
            gs = gs[["gene_sets"]]
            gs = gs[names(gs)[names(gs)%in%d$pathway]]

            ## Unlist all pathways and take unique genes
            gs = unique(unlist(gs))
            cat(g, "\n")
            cat(paste0("Enriched pathwas: ", nrow(d)), "\n")
            cat(paste0("Number of genes in enriched pathways: ", length(gs)), "\n")
            dd = syscan %>% subset(symbol%in% gs & Group==g) %>% count(sample) %>% mutate(group=g, genes_in_enriched_pathways=length(gs))
            genes_in_enriched_pathways = rbind(genes_in_enriched_pathways, dd)
        }
    }


    p = ggplot(genes_in_enriched_pathways, aes(x=group, y=n)) + geom_boxplot() +
        scale_y_continuous(limits = c(0, max(genes_in_enriched_pathways$n)), breaks = seq(0, max(genes_in_enriched_pathways$n), 20)) +
        ggtitle("Genes in enriched pathways per sample \n (one list per group approach)") +
        ylab("Genes (#)") + xlab("") +
        theme_boss() +
        theme(
            axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour = "black")
        )
    pdf("/Users/fc-8s-imac/athena/data/OAC/129_OAC/Tables_and_Plots/genes_in_enriched_pathways.pdf")
    print(p)
    dev.off()

}


checkRunning = function(dir="/Volumes/mourikisa/novel_driver_prediction/pancancer_model_training/automation/median_mode_imputation/"){
    cancer_types <- c("OAC", "BRCA", "COAD", "ESCA", "GBM",
                      "LAML", "LGG", "PAAD", "PCPG")
    for (ct in cancer_types){
        message(ct)
        res_dir = paste0(dir, "/", ct, "/Results/")
        kernels = c("linear", "polynomial", "radial", "sigmoid")
        for (k in kernels){
            fns = list.files(res_dir)
            fns = fns[grepl(k, fns)]
            message(paste0(k, ": ", length(fns)))
        }
    }
    cat("\n\n")
}

check_imputation = function(imputation_dir="/Volumes/mourikisa/data/"){
    require(mice)
    ## When I parallelize the imputation I store multiple objects objects
    ## I need to load one by one and gather the imputed datasets
    pars = 10 ## this is the number of objects run in parallel
    imp_datasets = 10 ## how many imputed datasets were produced by each process
    ## The imputed datasets have no entrez, you need to pull entrez from the geneProperties file
    ## Load geneProperties
    load("/Volumes/mourikisa/data/geneProperties_final.Rdata")
    entrez = geneProperties %>% select(Entrez)
    imputations = NULL ## gather all the imputed datasets and the initial one in a data frame
    count = 0
    for (i in 1:pars){
        cat(i, "\n")
        load(paste0(imputation_dir, "/geneProperties_final_miceImputed_object_", i, ".Rdata")) ## geneProperties_miceImputed
        if(i==1){ ## First time
            for (y in 0:imp_datasets){
                imputations = rbind(imputations, cbind(entrez, complete(geneProperties_miceImputed, y)) %>% data.frame() %>% mutate(iteration=count))
                count = count + 1
            }

        }else{
            for (y in 1:imp_datasets){
                imputations = rbind(imputations, cbind(entrez, complete(geneProperties_miceImputed, y)) %>% data.frame() %>% mutate(iteration=count))
                count = count + 1
            }
        }
    }
    rm(geneProperties_miceImputed)



    ## ---------------------------
    ## Check imputed features here
    ## ---------------------------
    ## Plot the distributions of each feature to make sure that the usage of the mean of the imputed values makes sense
    ## test properties to see what
    test = geneProperties%>%select(Entrez, age) %>% gather(feature, value, -Entrez) %>% subset(is.na(value)) %>% mutate(missing=T)
    for (i in 1:100){
        cat(i, "\n")
        d = imputations %>% subset(iteration==i) %>% select(Entrez, age) %>% gather(feature, value, -Entrez) %>% subset(Entrez%in%test$Entrez[test$missing==T]) %>% mutate(missing=F)
        test = rbind(test, d)
    }
    test = test %>% subset(!is.na(value)) %>% left_join(geneProperties%>%select(Entrez, cancer_type)) %>% count(cancer_type, value) %>% mutate(total=test%>%subset(!is.na(value))%>%nrow) %>% mutate(frac=n/total)
    ggplot(test, aes(x=value, y=frac, fill=cancer_type)) + geom_bar(stat="identity", position="dodge")


    ## I need to take the mean of the imputed values
    ## To do this I can utilize gather and spread functions
    ## I need to do that feature-by-feature because the type of variables are different
    ## The idea is i gather on each feature get the mean of the imputated values, then I spread and store all the data frames in a list
    ## Ten I merge all the data frames in the list to make a single data frame
    df_list = list()
    for (f in colnames(imputations)[2:(length(imputations)-1)]){
        cat(f, "\t")
        cat(class(imputations[,f]), "\n")
        ## I gather all the imputations without number 0 cause this is the missing data!
        if(class(imputations[,f])=="factor"){
            d = imputations %>% subset(iteration!=0) %>% select_("Entrez", f) %>% gather(feature, value, -Entrez)
            ## here we will take the majority
            d = d %>% group_by(Entrez, feature) %>% group_by(Entrez, feature) %>% summarise(value=names(which.max(table(value)))) %>% ungroup
        }else if(class(imputations[,f])%in%c("numeric", "integer")){
            d = imputations %>% subset(iteration!=0) %>% select_("Entrez", f) %>% gather(feature, value, -Entrez)
            ## here we will take the mean
            if (f=="betweenness"){ ## I dont round in betweenness
                d = d %>% group_by(Entrez, feature) %>% summarise(value=round(mean(value), digits=3)) %>% ungroup
            }else{
                d = d %>% group_by(Entrez, feature) %>% summarise(value=round(mean(value))) %>% ungroup
            }
        }
        d = d %>% spread(feature, value)
        df_list[[f]] = d
    }
    ## and now merge the tables into one
    merged = NULL
    for (i in 1:length(df_list)){
        cat(i, "\n")
        if(i==1){
            merged = df_list[[i]]
        }else{
            merged = merged %>% left_join(df_list[[i]])
        }
    }
    merged = merged %>% data.frame()
    ## sort it similarly to geneProperties
    merged = merged[order(match(merged[,"Entrez"],geneProperties[,"Entrez"])),]
    rownames(merged) <- seq(length=nrow(merged))
    save(geneProperties_miceImputed, file="/Volumes/mourikisa/data/geneProperties_final_miceImputed.Rdata")

    ## Load also gene info to get  cgc and rst annotation
    load("/Volumes/mourikisa/data/geneInfoNCG5.Rdata")
    geneInfo = geneInfo %>% subset(duplicability!=-1) %>% select(entrez, symbol, cancer_type, cancer_dom, cancer_rec)

    ## And the actual geneProperties cause we need the Entrez
    load("/Volumes/mourikisa/data/geneProperties_final.Rdata")
    geneProperties = geneProperties %>% left_join(geneInfo, by=c("Entrez"="entrez"))
    merged = merged %>% left_join(geneInfo, by=c("Entrez"="entrez"))

    ## We check the stats once for the geneProperties with the missing values and once for the imputed data
    overall_stats = NULL

    ## -----------------------------------
    ## For the raw data with mising values
    ## -----------------------------------
    ## Checking Length
    test = geneProperties %>% subset(!is.na(Length.fullrefseq) & cancer_type!="can") %>% select(Length.fullrefseq, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$Length.fullrefseq %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$Length.fullrefseq %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$Length.fullrefseq, test%>%subset(cancer_type=="rst")%>%.$Length.fullrefseq)$p.value
    tb = data.frame(data="Raw", description="Gene length (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of RefSeq length (aa) for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Check Duplicability
    test = geneProperties %>% subset(!is.na(duplicated) & cancer_type!="can" & (cancer_rec==1 | cancer_type=="rst")) %>% select(duplicated, cancer_type) %>% mutate(cancer_type=ifelse(cancer_type=="cgc", "rec", cancer_type))
    perc.cg = 100*(test %>% subset(duplicated==0 & cancer_type=="rec") %>% nrow)/(test %>% subset(cancer_type=="rec") %>% nrow)
    perc.rst = 100*(test %>% subset(duplicated==0 & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"duplicated"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Raw", description="Gene duplicability (CGC-recessive vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of singleton genes in CGC-recessive and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Checking WGD
    test = geneProperties %>% subset(!is.na(WGD) & cancer_type!="can") %>% select(WGD, cancer_type)
    perc.cg = 100*(test %>% subset(WGD==1 & cancer_type=="cgc") %>% nrow)/(test %>% subset(!is.na(WGD) & cancer_type=="cgc") %>% nrow)
    perc.rst = 100*(test %>% subset(WGD==1 & cancer_type=="rst") %>% nrow)/(test %>% subset(!is.na(WGD) & cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"WGD"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Raw", description="Ohnolog (CGC vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of genes underwent WGD in CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Number of domains between CGC and rest
    test = geneProperties %>% subset(!is.na(alldomains) & cancer_type!="can") %>% select(alldomains, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$alldomains %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$alldomains %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$alldomains, test%>%subset(cancer_type=="rst")%>%.$alldomains)$p.value
    tb = data.frame(data="Raw", description="Number of protein domains (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of alldomains for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Check hubs
    test = geneProperties %>% subset(!is.na(hub) & cancer_type!="can") %>% select(hub, cancer_type)
    perc.cg = 100*(test %>% subset(hub==1 & cancer_type=="cgc") %>% nrow)/(test %>% subset(cancer_type=="cgc") %>% nrow)
    perc.rst = 100*(test %>% subset(hub==1 & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"hub"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Raw", description="Hub (CGC vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of hubs in CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Check Central
    test = geneProperties %>% subset(!is.na(central) & cancer_type!="can") %>% select(central, cancer_type)
    perc.cg = 100*(test %>% subset(central==1 & cancer_type=="cgc") %>% nrow)/(test %>% subset(cancer_type=="cgc") %>% nrow)
    perc.rst = 100*(test %>% subset(central==1 & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"central"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Raw", description="Central (CGC vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of central proteins in CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Recessive cancer genes are of ancient evolutionary origin
    test = geneProperties %>% subset(!is.na(age) & cancer_type!="can" & (cancer_rec==1 | cancer_type=="rst")) %>% select(age, cancer_type) %>% mutate(cancer_type=ifelse(cancer_type=="cgc", "rec", cancer_type))
    perc.cg = 100*(test %>% subset(age=="old" & cancer_type=="rec") %>% nrow)/(test %>% subset(cancer_type=="rec") %>% nrow)
    perc.rst = 100*(test %>% subset(age=="old" & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"age"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Raw", description="Age (CGC-recessive vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of old genes in CGC-recessive and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Dominant cancer genes versus age
    test = geneProperties %>% subset(!is.na(age) & cancer_type!="can" & (cancer_dom==1 | cancer_type=="rst")) %>% select(age, cancer_type) %>% mutate(cancer_type=ifelse(cancer_type=="cgc", "dom", cancer_type))
    perc.cg = 100*(test %>% subset(age=="old" & cancer_type=="dom") %>% nrow)/(test %>% subset(cancer_type=="dom") %>% nrow)
    perc.rst = 100*(test %>% subset(age=="old" & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"age"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Raw", description="Age (CGC-dominant vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of old genes in CGC-dominant and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## All cancer genes versus age
    test = geneProperties %>% subset(!is.na(age)) %>% select(age, cancer_type) %>% mutate(cancer_type=ifelse(cancer_type=="cgc" | cancer_type=="can", "cancer_gene", cancer_type))
    perc.cg = 100*(test %>% subset(age=="old" & cancer_type=="cancer_gene") %>% nrow)/(test %>% subset(cancer_type=="cancer_gene") %>% nrow)
    perc.rst = 100*(test %>% subset(age=="old" & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"age"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Raw", description="Age (Cancer genes vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of old genes in Cancer genes and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Expression beadth
    test = geneProperties %>% subset(!is.na(tot.tissues) & cancer_type!="can") %>% select(tot.tissues, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$tot.tissues %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$tot.tissues %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$tot.tissues, test%>%subset(cancer_type=="rst")%>%.$tot.tissues)$p.value
    tb = data.frame(data="Raw", description="Number of tissues expressed (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of tissues in which the genes are expressed for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Checking HiC
    test = geneProperties %>% subset(!is.na(hic) & cancer_type!="can") %>% select(hic, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$hic %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$hic %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$hic, test%>%subset(cancer_type=="rst")%>%.$hic)$p.value
    tb = data.frame(data="Raw", description="Chromatin state (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of HiC for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Checking miRNA
    test = geneProperties %>% subset(!is.na(mirna) & cancer_type!="can") %>% select(mirna, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$mirna %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$mirna %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$mirna, test%>%subset(cancer_type=="rst")%>%.$mirna)$p.value
    tb = data.frame(data="Raw", description="miRNA target (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of miRNA for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)


    ## ------------------------------------------------------
    ## For the median/mode imputed data without mising values
    ## ------------------------------------------------------
    ## Checking Length
    test = geneProperties_mmImputed %>% subset(!is.na(Length.fullrefseq) & cancer_type!="can") %>% select(Length.fullrefseq, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$Length.fullrefseq %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$Length.fullrefseq %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$Length.fullrefseq, test%>%subset(cancer_type=="rst")%>%.$Length.fullrefseq)$p.value
    tb = data.frame(data="Imputed_MM", description="Gene length (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of RefSeq length (aa) for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Check Duplicability
    test = geneProperties_mmImputed %>% subset(!is.na(duplicated) & cancer_type!="can" & (cancer_rec==1 | cancer_type=="rst")) %>% select(duplicated, cancer_type) %>% mutate(cancer_type=ifelse(cancer_type=="cgc", "rec", cancer_type))
    perc.cg = 100*(test %>% subset(duplicated==0 & cancer_type=="rec") %>% nrow)/(test %>% subset(cancer_type=="rec") %>% nrow)
    perc.rst = 100*(test %>% subset(duplicated==0 & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"duplicated"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MM", description="Gene duplicability (CGC-recessive vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of singleton genes in CGC-recessive and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Checking WGD
    test = geneProperties_mmImputed %>% subset(!is.na(WGD) & cancer_type!="can") %>% select(WGD, cancer_type)
    perc.cg = 100*(test %>% subset(WGD==1 & cancer_type=="cgc") %>% nrow)/(test %>% subset(!is.na(WGD) & cancer_type=="cgc") %>% nrow)
    perc.rst = 100*(test %>% subset(WGD==1 & cancer_type=="rst") %>% nrow)/(test %>% subset(!is.na(WGD) & cancer_type=="rst") %>% nrow)
     cp = chisq.test(table(test[,"WGD"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MM", description="Ohnolog (CGC vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of genes underwent WGD in CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Number of domains between CGC and rest
    test = geneProperties_mmImputed %>% subset(!is.na(alldomains) & cancer_type!="can") %>% select(alldomains, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$alldomains %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$alldomains %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$alldomains, test%>%subset(cancer_type=="rst")%>%.$alldomains)$p.value
    tb = data.frame(data="Imputed_MM", description="Number of protein domains (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of alldomains for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Check hubs
    test = geneProperties_mmImputed %>% subset(!is.na(hub) & cancer_type!="can") %>% select(hub, cancer_type)
    perc.cg = 100*(test %>% subset(hub==1 & cancer_type=="cgc") %>% nrow)/(test %>% subset(cancer_type=="cgc") %>% nrow)
    perc.rst = 100*(test %>% subset(hub==1 & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"hub"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MM", description="Hub (CGC vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of hubs in CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Check Central
    test = geneProperties_mmImputed %>% subset(!is.na(central) & cancer_type!="can") %>% select(central, cancer_type)
    perc.cg = 100*(test %>% subset(central==1 & cancer_type=="cgc") %>% nrow)/(test %>% subset(cancer_type=="cgc") %>% nrow)
    perc.rst = 100*(test %>% subset(central==1 & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"central"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MM", description="Central (CGC vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of central proteins in CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Recessive cancer genes are of ancient evolutionary origin
    test = geneProperties_mmImputed %>% subset(!is.na(age) & cancer_type!="can" & (cancer_rec==1 | cancer_type=="rst")) %>% select(age, cancer_type) %>% mutate(cancer_type=ifelse(cancer_type=="cgc", "rec", cancer_type))
    perc.cg = 100*(test %>% subset(age=="old" & cancer_type=="rec") %>% nrow)/(test %>% subset(cancer_type=="rec") %>% nrow)
    perc.rst = 100*(test %>% subset(age=="old" & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"age"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MM", description="Age (CGC-recessive vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of old genes in CGC-recessive and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Dominant cancer genes versus age
    test = geneProperties_mmImputed %>% subset(!is.na(age) & cancer_type!="can" & (cancer_dom==1 | cancer_type=="rst")) %>% select(age, cancer_type) %>% mutate(cancer_type=ifelse(cancer_type=="cgc", "dom", cancer_type))
    perc.cg = 100*(test %>% subset(age=="old" & cancer_type=="dom") %>% nrow)/(test %>% subset(cancer_type=="dom") %>% nrow)
    perc.rst = 100*(test %>% subset(age=="old" & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"age"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MM", description="Age (CGC-dominant vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of old genes in CGC-dominant and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## All cancer genes versus age
    test = geneProperties_mmImputed %>% subset(!is.na(age)) %>% select(age, cancer_type) %>% mutate(cancer_type=ifelse(cancer_type=="cgc" | cancer_type=="can", "cancer_gene", cancer_type))
    perc.cg = 100*(test %>% subset(age=="old" & cancer_type=="cancer_gene") %>% nrow)/(test %>% subset(cancer_type=="cancer_gene") %>% nrow)
    perc.rst = 100*(test %>% subset(age=="old" & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"age"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MM", description="Age (Cancer genes vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of old genes in Cancer genes and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Expression beadth
    test = geneProperties_mmImputed %>% subset(!is.na(tot.tissues) & cancer_type!="can") %>% select(tot.tissues, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$tot.tissues %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$tot.tissues %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$tot.tissues, test%>%subset(cancer_type=="rst")%>%.$tot.tissues)$p.value
    tb = data.frame(data="Imputed_MM", description="Number of tissues expressed (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of tissues in which the genes are expressed for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Checking HiC
    test = geneProperties_mmImputed %>% subset(!is.na(hic) & cancer_type!="can") %>% select(hic, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$hic %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$hic %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$hic, test%>%subset(cancer_type=="rst")%>%.$hic)$p.value
    tb = data.frame(data="Imputed_MM", description="Chromatin state (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of HiC for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Checking miRNA
    test = geneProperties_mmImputed %>% subset(!is.na(mirna) & cancer_type!="can") %>% select(mirna, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$mirna %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$mirna %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$mirna, test%>%subset(cancer_type=="rst")%>%.$mirna)$p.value
    tb = data.frame(data="Imputed_MM", description="miRNA target (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of miRNA for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)


    ## ------------------------------------------------
    ## For the mice-imputed data without mising values
    ## ------------------------------------------------
    ## Checking Length
    test = geneProperties_miceImputed %>% subset(!is.na(Length.fullrefseq) & cancer_type!="can") %>% select(Length.fullrefseq, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$Length.fullrefseq %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$Length.fullrefseq %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$Length.fullrefseq, test%>%subset(cancer_type=="rst")%>%.$Length.fullrefseq)$p.value
    tb = data.frame(data="Imputed_MICE", description="Gene length (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of RefSeq length (aa) for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Check Duplicability
    test = geneProperties_miceImputed %>% subset(!is.na(duplicated) & cancer_type!="can" & (cancer_rec==1 | cancer_type=="rst")) %>% select(duplicated, cancer_type) %>% mutate(cancer_type=ifelse(cancer_type=="cgc", "rec", cancer_type))
    perc.cg = 100*(test %>% subset(duplicated==0 & cancer_type=="rec") %>% nrow)/(test %>% subset(cancer_type=="rec") %>% nrow)
    perc.rst = 100*(test %>% subset(duplicated==0 & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"duplicated"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MICE", description="Gene duplicability (CGC-recessive vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of singleton genes in CGC-recessive and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Checking WGD
    test = geneProperties_miceImputed %>% subset(!is.na(WGD) & cancer_type!="can") %>% select(WGD, cancer_type)
    perc.cg = 100*(test %>% subset(WGD==1 & cancer_type=="cgc") %>% nrow)/(test %>% subset(!is.na(WGD) & cancer_type=="cgc") %>% nrow)
    perc.rst = 100*(test %>% subset(WGD==1 & cancer_type=="rst") %>% nrow)/(test %>% subset(!is.na(WGD) & cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"WGD"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MICE", description="Ohnolog (CGC vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of genes underwent WGD in CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Number of domains between CGC and rest
    test = geneProperties_miceImputed %>% subset(!is.na(alldomains) & cancer_type!="can") %>% select(alldomains, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$alldomains %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$alldomains %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$alldomains, test%>%subset(cancer_type=="rst")%>%.$alldomains)$p.value
    tb = data.frame(data="Imputed_MICE", description="Number of protein domains (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of alldomains for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Check hubs
    test = geneProperties_miceImputed %>% subset(!is.na(hub) & cancer_type!="can") %>% select(hub, cancer_type)
    perc.cg = 100*(test %>% subset(hub==1 & cancer_type=="cgc") %>% nrow)/(test %>% subset(cancer_type=="cgc") %>% nrow)
    perc.rst = 100*(test %>% subset(hub==1 & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"hub"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MICE", description="Hub (CGC vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of hubs in CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Check Central
    test = geneProperties_miceImputed %>% subset(!is.na(central) & cancer_type!="can") %>% select(central, cancer_type)
    perc.cg = 100*(test %>% subset(central==1 & cancer_type=="cgc") %>% nrow)/(test %>% subset(cancer_type=="cgc") %>% nrow)
    perc.rst = 100*(test %>% subset(central==1 & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"central"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MICE", description="Central (CGC vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of central proteins in CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Recessive cancer genes are of ancient evolutionary origin
    test = geneProperties_miceImputed %>% subset(!is.na(age) & cancer_type!="can" & (cancer_rec==1 | cancer_type=="rst")) %>% select(age, cancer_type) %>% mutate(cancer_type=ifelse(cancer_type=="cgc", "rec", cancer_type))
    perc.cg = 100*(test %>% subset(age=="old" & cancer_type=="rec") %>% nrow)/(test %>% subset(cancer_type=="rec") %>% nrow)
    perc.rst = 100*(test %>% subset(age=="old" & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"age"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MICE", description="Age (CGC-recessive vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of old genes in CGC-recessive and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Dominant cancer genes versus age
    test = geneProperties_miceImputed %>% subset(!is.na(age) & cancer_type!="can" & (cancer_dom==1 | cancer_type=="rst")) %>% select(age, cancer_type) %>% mutate(cancer_type=ifelse(cancer_type=="cgc", "dom", cancer_type))
    perc.cg = 100*(test %>% subset(age=="old" & cancer_type=="dom") %>% nrow)/(test %>% subset(cancer_type=="dom") %>% nrow)
    perc.rst = 100*(test %>% subset(age=="old" & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"age"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MICE", description="Age (CGC-dominant vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of old genes in CGC-dominant and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## All cancer genes versus age
    test = geneProperties_miceImputed %>% subset(!is.na(age)) %>% select(age, cancer_type) %>% mutate(cancer_type=ifelse(cancer_type=="cgc" | cancer_type=="can", "cancer_gene", cancer_type))
    perc.cg = 100*(test %>% subset(age=="old" & cancer_type=="cancer_gene") %>% nrow)/(test %>% subset(cancer_type=="cancer_gene") %>% nrow)
    perc.rst = 100*(test %>% subset(age=="old" & cancer_type=="rst") %>% nrow)/(test %>% subset(cancer_type=="rst") %>% nrow)
    cp = chisq.test(table(test[,"age"], test[,"cancer_type"]))$p.value
    tb = data.frame(data="Imputed_MICE", description="Age (Cancer genes vs Rest)", cg=perc.cg, rst=perc.rst, test="chi-square test", p_value=cp, Note="cg and rst are percentages of old genes in Cancer genes and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Expression beadth
    test = geneProperties_miceImputed %>% subset(!is.na(tot.tissues) & cancer_type!="can") %>% select(tot.tissues, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$tot.tissues %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$tot.tissues %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$tot.tissues, test%>%subset(cancer_type=="rst")%>%.$tot.tissues)$p.value
    tb = data.frame(data="Imputed_MICE", description="Number of tissues expressed (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of tissues in which the genes are expressed for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Checking HiC
    test = geneProperties_miceImputed %>% subset(!is.na(hic) & cancer_type!="can") %>% select(hic, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$hic %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$hic %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$hic, test%>%subset(cancer_type=="rst")%>%.$hic)$p.value
    tb = data.frame(data="Imputed_MICE", description="Chromatin state (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of HiC for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    ## Checking miRNA
    test = geneProperties_miceImputed %>% subset(!is.na(mirna) & cancer_type!="can") %>% select(mirna, cancer_type)
    cg = test %>% subset(cancer_type=="cgc") %>% .$mirna %>% median
    rst = test %>% subset(cancer_type=="rst") %>% .$mirna %>% median
    wp = wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$mirna, test%>%subset(cancer_type=="rst")%>%.$mirna)$p.value
    tb = data.frame(data="Imputed_MICE", description="miRNA target (CGC vs Rest)", cg=cg, rst=rst, test="Wilcoxon", p_value=wp, Note="cg and rst are the median numbers of miRNA for CGC and Rest respectively")
    overall_stats = rbind(overall_stats, tb)

    overall_stats = overall_stats %>% mutate(significant = ifelse(p_value<0.05, 'Significant', "Not Significant"))


}

## Impute in the following function was the geneproperties data for the 19014
## At the moment I run it manually
imputeData = function(){

  load("/Volumes/mourikisa/data/geneProperties_final.Rdata")

  ## Make factors the categorical variables
  geneProperties$duplicated = factor(geneProperties$duplicated)
  geneProperties$WGD = factor(geneProperties$WGD)
  geneProperties$hub = factor(geneProperties$hub)
  geneProperties$central = factor(geneProperties$central)
  geneProperties$age = factor(geneProperties$age)
  geneProperties$origin = factor(geneProperties$origin)
  geneProperties$exp.breadth = factor(geneProperties$exp.breadth)

  ## Explore missing values and patterns
  md.pattern(geneProperties)
  library(VIM)
  aggr(geneProperties[,c(3:5,7,10,14:15)], col=c('navyblue','yellow'),
       numbers=TRUE, sortVars=TRUE,
       labels=names(geneProperties[,c(3:5,7,10,14:15)]), cex.axis=.7,
       gap=3, ylab=c("Missing data","Pattern"))

  ## Define a mode function to fill in missing values with mode
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  ## Median & Mode imputation
  geneProperties_mmImputed = geneProperties

  ## Bring geneInfo in cause we need CGC-specific median and Rest-specific median
  ## Load also gene info to get  cgc and rst annotation
  load("/Volumes/mourikisa/data/geneInfoNCG5.Rdata")
  geneInfo = geneInfo %>% subset(duplicability!=-1) %>% select(entrez, symbol, cancer_type, cancer_dom, cancer_rec)
  geneProperties_mmImputed = geneProperties_mmImputed %>% left_join(geneInfo, by=c("Entrez"="entrez"))

  ## Median/Mode imputation
  ## Fill WGD
  geneProperties_mmImputed$WGD[is.na(geneProperties_mmImputed$WGD) & geneProperties_mmImputed$cancer_type=="cgc"] = Mode(geneProperties_mmImputed$WGD[geneProperties_mmImputed$cancer_type=="cgc"])
  geneProperties_mmImputed$WGD[is.na(geneProperties_mmImputed$WGD) & geneProperties_mmImputed$cancer_type=="can"] = Mode(geneProperties_mmImputed$WGD[geneProperties_mmImputed$cancer_type=="rst"])
  geneProperties_mmImputed$WGD[is.na(geneProperties_mmImputed$WGD) & geneProperties_mmImputed$cancer_type=="rst"] = Mode(geneProperties_mmImputed$WGD[geneProperties_mmImputed$cancer_type=="rst"])

  ## Fill alldomains
  geneProperties_mmImputed$alldomains[is.na(geneProperties_mmImputed$alldomains) & geneProperties_mmImputed$cancer_type=="cgc"] = median(geneProperties_mmImputed$alldomains[geneProperties_mmImputed$cancer_type=="cgc"], na.rm = T)
  geneProperties_mmImputed$alldomains[is.na(geneProperties_mmImputed$alldomains) & geneProperties_mmImputed$cancer_type=="can"] = median(geneProperties_mmImputed$alldomains[geneProperties_mmImputed$cancer_type=="rst"], na.rm = T)
  geneProperties_mmImputed$alldomains[is.na(geneProperties_mmImputed$alldomains) & geneProperties_mmImputed$cancer_type=="rst"] = median(geneProperties_mmImputed$alldomains[geneProperties_mmImputed$cancer_type=="rst"], na.rm = T)

  ## Fill degree
  geneProperties_mmImputed$degree[is.na(geneProperties_mmImputed$degree) & geneProperties_mmImputed$cancer_type=="cgc"] = median(geneProperties_mmImputed$degree[geneProperties_mmImputed$cancer_type=="cgc"], na.rm = T)
  geneProperties_mmImputed$degree[is.na(geneProperties_mmImputed$degree) & geneProperties_mmImputed$cancer_type=="can"] = median(geneProperties_mmImputed$degree[geneProperties_mmImputed$cancer_type=="rst"], na.rm = T)
  geneProperties_mmImputed$degree[is.na(geneProperties_mmImputed$degree) & geneProperties_mmImputed$cancer_type=="rst"] = median(geneProperties_mmImputed$degree[geneProperties_mmImputed$cancer_type=="rst"], na.rm = T)

  ## Fill betweenness
  geneProperties_mmImputed$betweenness[is.na(geneProperties_mmImputed$betweenness) & geneProperties_mmImputed$cancer_type=="cgc"] = median(geneProperties_mmImputed$betweenness[geneProperties_mmImputed$cancer_type=="cgc"], na.rm = T)
  geneProperties_mmImputed$betweenness[is.na(geneProperties_mmImputed$betweenness) & geneProperties_mmImputed$cancer_type=="can"] = median(geneProperties_mmImputed$betweenness[geneProperties_mmImputed$cancer_type=="rst"], na.rm = T)
  geneProperties_mmImputed$betweenness[is.na(geneProperties_mmImputed$betweenness) & geneProperties_mmImputed$cancer_type=="rst"] = median(geneProperties_mmImputed$betweenness[geneProperties_mmImputed$cancer_type=="rst"], na.rm = T)

  ## Fill hub
  geneProperties_mmImputed$hub[is.na(geneProperties_mmImputed$hub) & geneProperties_mmImputed$cancer_type=="cgc"] = Mode(geneProperties_mmImputed$hub[geneProperties_mmImputed$cancer_type=="cgc"])
  geneProperties_mmImputed$hub[is.na(geneProperties_mmImputed$hub) & geneProperties_mmImputed$cancer_type=="can"] = Mode(geneProperties_mmImputed$hub[geneProperties_mmImputed$cancer_type=="rst"])
  geneProperties_mmImputed$hub[is.na(geneProperties_mmImputed$hub) & geneProperties_mmImputed$cancer_type=="rst"] = Mode(geneProperties_mmImputed$hub[geneProperties_mmImputed$cancer_type=="rst"])

  ## Fill central
  geneProperties_mmImputed$central[is.na(geneProperties_mmImputed$central) & geneProperties_mmImputed$cancer_type=="cgc"] = Mode(geneProperties_mmImputed$central[geneProperties_mmImputed$cancer_type=="cgc"])
  geneProperties_mmImputed$central[is.na(geneProperties_mmImputed$central) & geneProperties_mmImputed$cancer_type=="can"] = Mode(geneProperties_mmImputed$central[geneProperties_mmImputed$cancer_type=="rst"])
  geneProperties_mmImputed$central[is.na(geneProperties_mmImputed$central) & geneProperties_mmImputed$cancer_type=="rst"] = Mode(geneProperties_mmImputed$central[geneProperties_mmImputed$cancer_type=="rst"])

  ## Fill age
  geneProperties_mmImputed$age[is.na(geneProperties_mmImputed$age) & geneProperties_mmImputed$cancer_type=="cgc"] = Mode(geneProperties_mmImputed$age[geneProperties_mmImputed$cancer_type=="cgc"])
  geneProperties_mmImputed$age[is.na(geneProperties_mmImputed$age) & geneProperties_mmImputed$cancer_type=="can"] = Mode(geneProperties_mmImputed$age[geneProperties_mmImputed$cancer_type=="rst"])
  geneProperties_mmImputed$age[is.na(geneProperties_mmImputed$age) & geneProperties_mmImputed$cancer_type=="rst"] = Mode(geneProperties_mmImputed$age[geneProperties_mmImputed$cancer_type=="rst"])

  ## Fill origin
  geneProperties_mmImputed$origin[is.na(geneProperties_mmImputed$origin) & geneProperties_mmImputed$cancer_type=="cgc"] = Mode(geneProperties_mmImputed$origin[geneProperties_mmImputed$cancer_type=="cgc"])
  geneProperties_mmImputed$origin[is.na(geneProperties_mmImputed$origin) & geneProperties_mmImputed$cancer_type=="can"] = Mode(geneProperties_mmImputed$origin[geneProperties_mmImputed$cancer_type=="rst"])
  geneProperties_mmImputed$origin[is.na(geneProperties_mmImputed$origin) & geneProperties_mmImputed$cancer_type=="rst"] = Mode(geneProperties_mmImputed$origin[geneProperties_mmImputed$cancer_type=="rst"])

  ## Fill exp.breadth
  geneProperties_mmImputed$exp.breadth[is.na(geneProperties_mmImputed$exp.breadth) & geneProperties_mmImputed$cancer_type=="cgc"] = Mode(geneProperties_mmImputed$exp.breadth[geneProperties_mmImputed$cancer_type=="cgc"])
  geneProperties_mmImputed$exp.breadth[is.na(geneProperties_mmImputed$exp.breadth) & geneProperties_mmImputed$cancer_type=="can"] = Mode(geneProperties_mmImputed$exp.breadth[geneProperties_mmImputed$cancer_type=="rst"])
  geneProperties_mmImputed$exp.breadth[is.na(geneProperties_mmImputed$exp.breadth) & geneProperties_mmImputed$cancer_type=="rst"] = Mode(geneProperties_mmImputed$exp.breadth[geneProperties_mmImputed$cancer_type=="rst"])

  ## Fill tot.tissue
  geneProperties_mmImputed$tot.tissues[is.na(geneProperties_mmImputed$tot.tissues) & geneProperties_mmImputed$cancer_type=="cgc"] = median(geneProperties_mmImputed$tot.tissues[geneProperties_mmImputed$cancer_type=="cgc"], na.rm = T)
  geneProperties_mmImputed$tot.tissues[is.na(geneProperties_mmImputed$tot.tissues) & geneProperties_mmImputed$cancer_type=="can"] = median(geneProperties_mmImputed$tot.tissues[geneProperties_mmImputed$cancer_type=="rst"], na.rm = T)
  geneProperties_mmImputed$tot.tissues[is.na(geneProperties_mmImputed$tot.tissues) & geneProperties_mmImputed$cancer_type=="rst"] = median(geneProperties_mmImputed$tot.tissues[geneProperties_mmImputed$cancer_type=="rst"], na.rm = T)

  ## Fill hic
  geneProperties_mmImputed$hic[is.na(geneProperties_mmImputed$hic) & geneProperties_mmImputed$cancer_type=="cgc"] = median(geneProperties_mmImputed$hic[geneProperties_mmImputed$cancer_type=="cgc"], na.rm = T)
  geneProperties_mmImputed$hic[is.na(geneProperties_mmImputed$hic) & geneProperties_mmImputed$cancer_type=="can"] = median(geneProperties_mmImputed$hic[geneProperties_mmImputed$cancer_type=="rst"], na.rm = T)
  geneProperties_mmImputed$hic[is.na(geneProperties_mmImputed$hic) & geneProperties_mmImputed$cancer_type=="rst"] = median(geneProperties_mmImputed$hic[geneProperties_mmImputed$cancer_type=="rst"], na.rm = T)

  ## Fill mirna
  geneProperties_mmImputed$mirna[is.na(geneProperties_mmImputed$mirna) & geneProperties_mmImputed$cancer_type=="cgc"] = median(geneProperties_mmImputed$mirna[geneProperties_mmImputed$cancer_type=="cgc"], na.rm = T)
  geneProperties_mmImputed$mirna[is.na(geneProperties_mmImputed$mirna) & geneProperties_mmImputed$cancer_type=="can"] = median(geneProperties_mmImputed$mirna[geneProperties_mmImputed$cancer_type=="rst"], na.rm = T)
  geneProperties_mmImputed$mirna[is.na(geneProperties_mmImputed$mirna) & geneProperties_mmImputed$cancer_type=="rst"] = median(geneProperties_mmImputed$mirna[geneProperties_mmImputed$cancer_type=="rst"], na.rm = T)


  ## I have problems with member of complex cause there is only 1 observation
  save(geneProperties_mmImputed, file="/Volumes/mourikisa/data/geneProperties_final_mmImputed.Rdata")


  ## PMM imputation using mice - here we predict everything
  ## You need to conver characters to factors in order to be imputed - they are already
  ## You can check the manuscript her - https://stat.ethz.ch/education/semesters/ss2012/ams/paper/mice.pdf
  ## I chose 10 because they say that 10-20 is usually sufficient
  geneProperties$duplicated = factor(geneProperties$duplicated)
  geneProperties$WGD = factor(geneProperties$WGD)
  geneProperties$hub = factor(geneProperties$hub)
  geneProperties$central = factor(geneProperties$central)
  geneProperties$age = factor(geneProperties$age)
  geneProperties$origin = factor(geneProperties$origin)
  geneProperties$exp.breadth = factor(geneProperties$exp.breadth)
  ## I parallelized it here. Check the impute.R in automation file, I run the following commnd
  geneProperties_miceImputed <- mice(geneProperties[2:15], m=10, maxit = 20)
  summary(geneProperties_miceImputed)
  save(geneProperties_miceImputed, file="/Users/thmourikis/Desktop/geneProperties_final_miceImputed_object.Rdata")

  ## Then we replace the NAs with the mean of the imputed values
  ## For this I need to utilize spread and gather functions



  ## Continuous variables
  densityplot(geneProperties_miceImputed)
  ## Categorical variables
  categs = geneProperties %>% Filter(f = is.factor) %>% names
  categs = categs[!(categs%in%c("Genic", "memberofcomplex"))]
  mv = complete(geneProperties_miceImputed, 0)
  imputed_data = NULL
  for (i in 1:10){
    cat(i, "\n")
    d = complete(geneProperties_miceImputed, i)
    for (cat in categs){
      initial = mv %>% select_(cat) %>% count_(cat)
      dd = d %>% select_(cat) %>% count_(cat) %>% left_join(initial %>% rename(initial=n)) %>%
        mutate(change=abs(n-initial))
      dd = dd %>% mutate(imputation=i, feature=cat)
      colnames(dd)[1] = "levels"
      imputed_data = rbind(dd, imputed_data)
    }
  }

  ggplot(imputed_data, aes(x=levels, y=change, group=feature, fill=feature)) +
    geom_bar(stat="identity", position = "dodge") + facet_wrap(~imputation, ncol = 2) +
    theme(
      axis.text.x=element_text(angle = 90)
    )

}


getAllDriversTCGA = function(data_dir = "~/athena/novel_driver_prediction/pancancer_model_training/automation/median_mode_imputation/",
                             ct=NULL){

    ## Get gene type
    ## Get the gene names from NCG
    load("~/athena/data/geneInfoNCG5.Rdata")
    load("~/athena/data/cancerGenesNCG5.Rdata")
    ## Fix gene info table from NCG
    geneInfo = geneInfo %>% select(entrez, symbol, cancer_type, gene_categories) %>% unique %>% rename(gene_type=cancer_type)
    ## Get a cancer gene with all the associated primary sites and cancer sites
    cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
        group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
                                       cancer_site=paste(unique(cancer_site), collapse=",")) %>%
        ungroup
    geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))

    if(is.null(ct)){
        cancer_types <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM",
                          "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD",
                          "LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM",
                          "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
    }else{
        cancer_types=ct
    }



    ## get the total driver cohort

    total_drivers = NULL
    for (ct in cancer_types){
        message(ct)
        load(paste0(data_dir,ct, "/training_set_noScale.Rdata"))
        training_ns = training_ns %>% tibble::rownames_to_column() %>% separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")
        load(paste0(data_dir,ct, "/validation_set_noScale.Rdata"))
        validation_ns = validation_ns %>% tibble::rownames_to_column() %>% separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")
        cohort = rbind(training_ns%>%select(-type)%>%mutate(type="training"), validation_ns%>%mutate(type="prediction"))
        total_drivers = rbind(total_drivers, cohort %>% data.frame())
    }
    total_drivers$entrez = as.numeric(total_drivers$entrez)
    return(total_drivers)

}
## This function takes as input the data frame from getAllDriversTCGA and gives samples with drivers per gene
inferDriversPerGene = function(d = NULL, cgc_fn="/Volumes/mourikisa/novel_driver_prediction/518_CGC_annotation_24_10_16F.xlsx"){

    if (is.null(d)){
        stop("Please provide input data")
    }
    ## Store the CGC drivers after refinement of their driver alterations
    cgc_drivers = NULL
    if(!is.null(cgc_fn)){
        ## Load the CGC that we will consider drivers
        message("Getting CGC annotation...")
        wb = loadWorkbook(cgc_fn)
        cgcs = readWorksheet(wb, sheet=1)
        ## refine to those that will be retained
        retained_CGC = cgcs %>% subset(RET.DIS=="RET")

        ## Refine the training set to include only the ones we decided to coonsider drivers
        message("Refining CGC drivers...")
        for (i in 1:nrow(retained_CGC)){

            en = retained_CGC$entrez[i]
            gene_type= retained_CGC$GENETICS[i]
            symbol = retained_CGC$symbol[i]
            message(paste(symbol, gene_type, sep = "\t"))
            if (gene_type=="O"){
                dd = d %>% subset(entrez==en) %>% subset(CNVGain==1 | no_NTDam_muts!=0 | no_GOF_muts!=0)
                cgc_drivers = rbind(cgc_drivers, dd)
            }else if (gene_type=="TS"){
                dd = d %>% subset(entrez==en) %>% subset(Copy_number==0 | no_TRUNC_muts!=0 | no_NTDam_muts!=0 | (Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)))
                cgc_drivers = rbind(cgc_drivers, dd)
            }else if (gene_type=="O/TS"){
                dd = d %>% subset(entrez==en)
                cgc_drivers = rbind(cgc_drivers, dd)
            }else{
                stop("Unknown gene type")
            }
        }
    }

    ## And now after the refinement of cgc
    ## rbind it with type=="prediction" (prediction set without the CGCs) and count samples per gene
    if(!is.null(cgc_drivers)){
        message("CGC refined alterations mode...")
        allDrivers = rbind(cgc_drivers, d%>%subset(type=="prediction"))
        ## Count per gene
        allDrivers = allDrivers %>% select(sample, entrez) %>% unique %>% group_by(entrez) %>% summarise(samples_with_drivers=n()) %>% ungroup()
    }else if (is.null(cgc_drivers)){
        message("CGC non-refined alterations mode...")
        allDrivers = d %>% select(sample, entrez) %>% unique %>% group_by(entrez) %>% summarise(samples_with_drivers=n()) %>% ungroup()
    }


    return(allDrivers)

}


## Recalculate score
reScore = function(path){

    load(paste0(path, "/scores.Rdata"))
    ## Get the distance from the decision boundary in each kernel
    all_preds = scores[["all_scores"]]
    all_preds = all_preds %>% select(cancer_type, sample, entrez, dv, label, prob_platt, gene_type, kernel, training)

    ## Get the best models
    bm = read.table(paste0(path, "/best_models.tsv"), header = T, sep = "\t", stringsAsFactors = F)
    if (!("kernel"%in%colnames(bm))){
        bm = bm %>% separate(analysis, into=c("kernel", "params"), sep="\\.", extra="merge", remove = F) %>%
            select(-params) %>% select(kernel, mean, var) %>% mutate(ratio=mean/var)
    }else{
        bm = bm %>% select(kernel, mean, var) %>% mutate(ratio=mean/var)
    }

    ## Get sample-specific rank per kernel
    all_preds = all_preds %>% arrange(sample, dv) %>% group_by(sample, kernel) %>% mutate(rank=row_number(-dv)) %>% ungroup
    all_preds = all_preds %>% group_by(sample) %>% mutate(N=max(rank), Nlog10=log10(N)) %>% ungroup

    ## bring the kernel information into the scores
    all_preds = all_preds %>% left_join(bm) %>% rename(cvs=mean)

    ## calculate new score
    #all_preds = all_preds %>% mutate(score=-log10(rank/N)*cvs)
    all_preds = all_preds %>% mutate(score=-log10(rank/N)*ratio)

    ## Sum up the score of each gene across kernels
    scores = all_preds %>% group_by(cancer_type, sample, entrez, gene_type, N, Nlog10) %>%
        summarise(sum_of_scores=sum(score), kernels_used=length(unique(kernel))) %>%
        #mutate(score=(1/(kernels_used*Nlog10))*sum_of_scores) %>% ungroup %>% data.frame()
        ungroup %>% group_by(sample) %>% mutate(max_score=max(sum_of_scores)) %>% mutate(score=sum_of_scores/max_score) %>%
        ungroup %>% data.frame()

    ## Add from how many kernels they have predicted as TRUE
    preds2kernels = all_preds %>% subset(label) %>% group_by(cancer_type, sample, entrez) %>% summarise(kernels_predicted=paste(unique(kernel), collapse=","), kernels_predicted_no=length(unique(kernel))) %>% ungroup
    scores = scores %>% left_join(preds2kernels)
    scores$kernels_predicted_no[is.na(scores$kernels_predicted_no)] = 0

    return(list(all_scores=all_preds, scores=scores))

}

saveReport = function (path, cohort=NULL, title=NULL, out_fn="sypcad_report.pdf", tissue=NULL) {

    if(is.null(tissue)){
        message("Tissue argument is null!")
    }

    pdf(file=paste0( path, "/", out_fn,sep=""),width = unit(12, "inches"), height = unit(14, "inches"))

    pushViewport(viewport(layout=grid.layout(10, 5, widths=c(0.05, 1, 0.05, 0.05, 1), heights=c(.05,.01,.15,.01, 0.1, 0.28,0.01,0.1, 0.28,.05))))

    grid.text( title , vp = viewport(layout.pos.row = 1, layout.pos.col = 3), gp=gpar(fontface='bold'), just="center")

    grid.text("A. General information", vp = viewport(layout.pos.row = 2, layout.pos.col = 1), gp=gpar(fontface='bold'), just="left")

    ## Get the table witht hte general info for the cancer type
    if(!is.null(cohort)){
        ## Make a table for the genes and print the table in the PDF
        ## This stays the same always
        training_g = cohort %>% subset(type=="C") %>% select(entrez) %>% unique %>% .$entrez
        training_s = cohort %>% subset(type=="C") %>% select(sample) %>% unique %>% .$sample
        training_t = cohort %>% subset(type=="C")
        validation_g = cohort %>% subset(type=="P") %>% select(entrez) %>% unique %>% .$entrez
        validation_s = cohort %>% subset(type=="P") %>% select(sample) %>% unique %>% .$sample
        validation_t = cohort %>% subset(type=="P")
        samples_drivers = unique(cohort$sample)
        samples_unexplained = setdiff(validation_s, training_s)
        d1 = data.frame(info=c("Training(g/s/t)", "Validation(g/s/t)","Samples (drivers)", "Samples (unexplained)"),
                       values=c(paste(length(training_g), length(training_s), nrow(training_t), sep="/"), paste(length(validation_g), length(validation_s), nrow(validation_t), sep="/"),
                                length(samples_drivers), length(samples_unexplained)))
        g1 = tableGrob(d1, rows = NULL, cols=NULL, vp=viewport(layout.pos.row=3,layout.pos.col=1:2))
        grid.draw(g1)

        cgc = cohort %>% subset(sys_can & gene_type=="cgc") %>% select(entrez) %>% unique %>% .$entrez
        ncg_cans = cohort %>% subset(sys_can & gene_type=="can") %>% select(entrez) %>% unique %>% .$entrez
        fps = cohort %>% subset(sys_can & NCG_FP) %>% select(entrez) %>% unique %>% .$entrez
        rst = cohort %>% subset(sys_can & gene_type=="rst") %>% select(entrez) %>% unique %>% .$entrez
        ## Load the Cancer genes file to distinguish between cancer specific and not
        load("/Volumes/mourikisa/data/geneInfoNCG5.Rdata")
        load("/Volumes/mourikisa/data/cancerGenesNCG5.Rdata")
        ## Fix gene info table from NCG
        geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
        ## Get a cancer gene with all the associated primary sites and cancer sites
        cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
            group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
                                           cancer_site=paste(unique(cancer_site), collapse=",")) %>%
            ungroup
        geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))

        ## Check if tissue is in NCG tissue
        ncg_sites = geneInfo %>% subset(!is.na(primary_site)) %>% select(primary_site) %>% mutate(primary_site=strsplit(primary_site, ",")) %>% unnest(primary_site) %>% unique() %>% .$primary_site
        if (!(tissue%in%ncg_sites)){
            stop("The tissue provided is not part of NCG5")
        }

        ## get the cancer specific ones
        cs_cgc = geneInfo %>% subset(entrez%in%cgc & grepl(tissue, primary_site)) %>% select(entrez) %>% unique %>% .$entrez
        cs_ncg_cans = geneInfo %>% subset(entrez%in%ncg_cans & grepl(tissue, primary_site)) %>% select(entrez) %>% unique %>% .$entrez
        cs_fps = geneInfo %>% subset(entrez%in%fps & grepl(tissue, primary_site)) %>% select(entrez) %>% unique %>% .$entrez

        d2 = data.frame(info=c("CGCs", "Cancer-specific CGCs", "Sys candidates (cans;FP)", "Cancer-specific (cans;FP)", "Sys candidates (rst)"),
                       values=c(length(cgc), length(cs_cgc), paste0(length(ncg_cans), ";", length(fps)), paste0(length(cs_ncg_cans), ";", length(cs_fps)), length(rst)))
        g2 = tableGrob(d2, rows = NULL, cols=NULL, vp=viewport(layout.pos.row=3,layout.pos.col=4:5))
        grid.draw(g2)



    }

    grid.text("B. Score distribution (overall)", vp = viewport(layout.pos.row = 4, layout.pos.col = 1), gp=gpar(fontface='bold'), just="left")

    ## Get the score distributions
    if(!is.null(cohort)){
        if(grepl("OAC", title)){
            toPlot = rbind(cohort %>% subset(type=="P")%>% mutate(type="All scores"),
                           cohort %>% subset(type=="P" & sys_can)%>% mutate(type="Systems candidates"))
            toPlot = toPlot %>% mutate(cl=ifelse(symbol=="KAT2A", "KAT2A", ifelse(symbol=="KAT2B", "KAT2B", "Rest")))
            tb = toPlot %>% mutate(type=ifelse(type=="All scores", "All", "Sys cans")) %>%
                group_by(type) %>%
                summarise(min=format(min(score), scientific = T, digits = 2), median=format(median(score), scientific = T, digits = 2),
                          mean=format(mean(score), scientific = T, digits = 2), max=format(max(score), scientific = T, digits = 2))
            p = ggplot(toPlot, aes(x=type, y=score)) + geom_boxplot(outlier.colour = NA) + geom_point(aes(color=cl, alpha=.5), position = position_jitter(width = 0.2)) +
                scale_color_manual(values=c("red", "green", "black"), name="Genes", breaks=c("KAT2A","KAT2B")) +
                theme_boss() +
                theme(
                    axis.line.x = element_line(color="black"),
                    axis.line.y = element_line(color="black"),
                    axis.text.x = element_text(size=12),
                    axis.text.y = element_text(size=12),
                    axis.title.x = element_text(size=14)
                ) +
                xlab("")
            g = tableGrob(tb, rows = NULL, vp=viewport(layout.pos.row=5,layout.pos.col=1:2))
            grid.draw(g)
            print(p, vp=viewport(layout.pos.row=6,layout.pos.col=1:2))
        }else{
            toPlot = rbind(cohort %>% subset(type=="P")%>% mutate(type="All scores"),
                          cohort %>% subset(type=="P" & sys_can)%>% mutate(type="Systems candidates"))
            tb = toPlot %>% mutate(type=ifelse(type=="All scores", "All", "Sys cans")) %>%
                group_by(type) %>%
                summarise(min=format(min(score), scientific = T, digits = 2), median=format(median(score), scientific = T, digits = 2),
                          mean=format(mean(score), scientific = T, digits = 2), max=format(max(score), scientific = T, digits = 2))
            p = ggplot(toPlot, aes(x=type, y=score)) + geom_boxplot(outlier.colour = NA) + #geom_point(position = position_jitter(width = 0.2)) +
                theme_boss() +
                theme(
                    axis.line.x = element_line(color="black"),
                    axis.line.y = element_line(color="black"),
                    axis.text.x = element_text(size=12),
                    axis.text.y = element_text(size=12),
                    axis.title.x = element_text(size=14)
                ) +
                xlab("")
            g = tableGrob(tb, rows = NULL, vp=viewport(layout.pos.row=5,layout.pos.col=1:2))
            grid.draw(g)
            print(p, vp=viewport(layout.pos.row=6,layout.pos.col=1:2))
        }
    }

    grid.text("C. Mean score within samples (only sys cans)", vp = viewport(layout.pos.row = 4, layout.pos.col = 4), gp=gpar(fontface='bold'), just="left")

    ## Get the mean score of each sample
    if(!is.null(cohort)){

        toPlot = cohort %>% subset(type=="P" & sys_can) %>% group_by(sample) %>% summarise(mean_score=mean(score))
        toPlot$type = paste0("samples(",nrow(toPlot), ")")
        samples2sys_cans = cohort %>% subset(type=="P" & sys_can) %>% group_by(sample) %>% summarise(n=n())
        toPlot = toPlot %>% left_join(samples2sys_cans)
        tb = cbind(cohort %>% subset(type=="P" & sys_can) %>% group_by(sample) %>% summarise(mean_score=mean(score)) %>% .$mean_score %>% summary(),
                    cohort %>% subset(type=="P" & sys_can) %>% group_by(sample) %>% summarise(n=n()) %>% .$n %>% summary()) %>% data.frame()
        colnames(tb) = c("Mean score", "Sys cans")
        tb = tb %>% t()
        tb = tb[,c(1,3,4,6)]
        tb = tb %>% data.frame() %>% tibble::rownames_to_column()
        colnames(tb) = c("type", "min", "median", "mean", "max")
        tb = tb %>% mutate(min=format(min, scientific = T, digits = 2), median=format(median, scientific = T, digits = 2),
                           mean=format(mean, scientific = T, digits = 2), max=format(max, scientific = T, digits = 2))

        p = ggplot(toPlot, aes(x=type, y=mean_score)) + geom_boxplot(outlier.colour = NA) + geom_point(aes(colour = factor(n)), position = position_jitter(width = 0.2)) +
            scale_color_manual(values = c("red", "orange", "green", "blue"), name="Number of\nsys cans", breaks=c(1,2,3,4)) +
            theme_boss() +
            theme(
                axis.line.x = element_line(color="black"),
                axis.line.y = element_line(color="black"),
                axis.text.x = element_text(size=12),
                axis.text.y = element_text(size=12),
                axis.title.x = element_text(size=14),
                legend.title=element_blank()
            ) + xlab("") + ylab("Mean sys candidate score within each sample")
            scale_colour_discrete(name="Sys\ncandidates")
        g = tableGrob(tb, rows = NULL, vp=viewport(layout.pos.row=5,layout.pos.col=4:5))
        grid.draw(g)
        print(p, vp=viewport(layout.pos.row=6,layout.pos.col=4:5))
    }

    ncg_cans = cohort %>% subset(sys_can & gene_type=="can") %>% select(entrez) %>% unique %>% .$entrez
    fps = cohort %>% subset(sys_can & NCG_FP) %>% .$entrez
    rst = cohort %>% subset(sys_can & gene_type=="rst") %>% select(entrez) %>% unique %>% .$entrez
    grid.text(paste0("D. Recurrence of sys candidates (",length(ncg_cans)+length(rst) ,") across samples"), vp = viewport(layout.pos.row = 7, layout.pos.col = 1), gp=gpar(fontface='bold'), just="left")

    ## Get in how many samples the sys cans are mutated
    if(!is.null(cohort)){
        samples =  cohort %>% subset(sys_can & type=="P") %>% select(sample) %>% unique %>% nrow
        toPlot = cohort %>% subset(type=="P" & sys_can) %>% group_by(symbol) %>% summarise(n=n())
        p = ggplot(toPlot, aes(x=n)) + geom_histogram(binwidth=1) + theme_boss() +
            theme(
                axis.line.x = element_line(color="black"),
                axis.line.y = element_line(color="black"),
                axis.text.x = element_text(size=12),
                axis.text.y = element_text(size=12),
                axis.title.x = element_text(size=14)
            ) + xlab("Number of samples") +
            scale_x_continuous(limits = c(0, samples))
        tb = toPlot %>%  .$n %>% summary() %>% cbind() %>% data.frame()
        colnames(tb) = "samples"
        tb = tb %>% t() %>% data.frame()
        tb = tb[,c(1,3,4,6)]
        colnames(tb) = c("min", "median", "mean", "max")
        g = tableGrob(tb, rows = NULL, vp=viewport(layout.pos.row=8,layout.pos.col=1:2))
        grid.draw(g)
        print(p, vp=viewport(layout.pos.row=9,layout.pos.col=1:2))
    }

    ## Get genes that are recurrent and plot mean score versus variance
    all_genes = cohort %>% subset(type=="P") %>% group_by(entrez, symbol) %>% summarise(n=n(), m=mean(score), r=max(score)-min(score)) %>% ungroup
    rec_genes = all_genes %>% subset(n>1)
    rec_genes = rec_genes %>% mutate(cl=ifelse(entrez%in%c(ncg_cans, rst), TRUE, FALSE))
    grid.text(paste0("E. Gene score across samples (", nrow(rec_genes),"/",nrow(all_genes), " genes)") , vp = viewport(layout.pos.row = 7, layout.pos.col = 4), gp=gpar(fontface='bold'), just="left")

    if(!is.null(cohort)){
        qs = quantile(rec_genes$r, seq(0,1,.05))[c(16,20)]
        tb = rbind(rec_genes %>%  .$r %>% summary() %>% cbind() %>% data.frame(),
                   qs %>% data.frame())
        tb[,1] = format(tb[,1], scientific = T, digits = 2)
        colnames(tb) = "score range"
        tb = tb %>% t() %>% data.frame()
        tb = tb[,c(1,3,4,7, 8, 6)]
        tb = tb %>% data.frame() %>% tibble::rownames_to_column()
        colnames(tb) = c("type", "min", "median", "mean", "75%", "95%", "max")

        p = ggplot(rec_genes, aes(x=r, y=m)) + geom_point(aes(color=cl)) +
#             geom_text(data=subset(rec_genes, r > 10),
#                       aes(x=r,y=m+0.5,label=symbol), size=3) +
            theme_boss() +
            theme(
                axis.line.x = element_line(color="black"),
                axis.line.y = element_line(color="black"),
                axis.text.x = element_text(size=12),
                axis.text.y = element_text(size=12),
                axis.title.x = element_text(size=14)
            ) + xlab("Score range") + ylab("Mean gene score across samples") +
            scale_color_manual(values = c("black", "red"), name="Sys candidates")

        g = tableGrob(tb, rows = NULL, vp=viewport(layout.pos.row=8,layout.pos.col=4:5))
        grid.draw(g)
        print(p, vp=viewport(layout.pos.row=9,layout.pos.col=4:5))
    }


    dev.off()
    message(paste0( "Results saved at: ", path, "/", out_fn,sep=""))
}

## this function adds information from gtex, and takes into account gene
## expression when calculating the rank of each prediction within each sample
getSysCans = function(path=NULL, cancer_type=NULL, tissue=NULL, scores_fn="scores.Rdata",
                      fp_dir="~/Mountpoints/rosalind_lustre/mourikisa/data/NCG_false_positives.txt",
                      cgc_fn="~/Mountpoints/rosalind_lustre/mourikisa/data/518_CGC_annotation_2211.xlsx",
                      previous.studies="~/Mountpoints/rosalind_lustre/mourikisa/data/CGCs_to_be_considered_OAC.tsv",
                      exclude_expr=c("Not Expressed")){

    require(readxl)

    if(is.null(path) | is.null(cancer_type) | is.null(tissue)){
        stop("Missing arguments")
    }

    ## Get the data we need for the definition of the sys cans
    ## Load geneInfo to get the symbols as well
    load("~/Mountpoints/rosalind_lustre/mourikisa/data/geneInfoNCG5.Rdata")
    geneInfoNCG5 = geneInfo %>% select(entrez, symbol) %>% unique

    ## Load the CGC that we will consider drivers
    message("Getting CGC annotation...")
    #wb = loadWorkbook(cgc_fn) ## changed function
    #cgcs = readWorksheet(wb, sheet=1) ## changed funstion
    cgcs = read_excel(cgc_fn, 1)
    ## refine to those that will be retained
    retained_CGC = cgcs %>% subset(RET.DIS=="RET")
    retained_CGC_trans = retained_CGC %>% subset(grepl("Trans", Drivers.to.be.retained))
    ## Get training and validation set
    load(paste0(path, "/training_set_noScale.Rdata"))
    training_ns = training_ns%>%tibble::rownames_to_column()%>%
        separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>%
        mutate(entrez=as.numeric(entrez)) %>% data.frame()
    ## Refine the training set to include only the ones we decided to coonsider drivers
    message("Refining CGC drivers...")
    cgc_drivers = NULL
    for (i in 1:nrow(retained_CGC)){
        en = retained_CGC$entrez[i]
        gene_type= retained_CGC$GENETICS[i]
        symbol = retained_CGC$symbol[i]
        retained_drivers = retained_CGC$Drivers.to.be.retained[i]
        #message(paste(symbol, gene_type, sep = "\t"))
        if (gene_type=="O"){
            if(grepl("Trans", retained_drivers)){
                dd = training_ns %>% subset(entrez==en) %>% subset(CNVGain==1 | no_NTDam_muts!=0 | no_GOF_muts!=0 | BND==1 | INV==1 | INS==1)
            }else{
                dd = training_ns %>% subset(entrez==en) %>% subset(CNVGain==1 | no_NTDam_muts!=0 | no_GOF_muts!=0)
            }
            cgc_drivers = rbind(cgc_drivers, dd)
        }else if (gene_type=="TS"){
            if(grepl("Trans", retained_drivers)){
                dd = training_ns %>% subset(entrez==en) %>% subset(Copy_number==0 | no_TRUNC_muts!=0 | no_NTDam_muts!=0 | (Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)) | BND==1 | INV==1 | INS==1)
            }else{
                dd = training_ns %>% subset(entrez==en) %>% subset(Copy_number==0 | no_TRUNC_muts!=0 | no_NTDam_muts!=0 | (Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)))
            }
            cgc_drivers = rbind(cgc_drivers, dd)
        }else if (gene_type=="O/TS"){
            if(grepl("Trans", retained_drivers)){
                ## All drivers are considered
                dd = training_ns %>% subset(entrez==en)
            }else{
                ## Here we need to exlude all the genes with only SVs
                dd = training_ns %>% subset(entrez==en) %>% subset(!(no_TRUNC_muts==0 & no_NTDam_muts==0 & no_GOF_muts==0 & Copy_number!=0 & CNVGain==0 & (BND!=0 | INS!=0 | INV!=0)))
            }
            cgc_drivers = rbind(cgc_drivers, dd)
        }else{
            stop("Unknown gene type")
        }
    }

    ## Now add in the CGC drivers the CGCs identified from previous studies
    previous = read.delim(previous.studies, header = T, sep="\t")
    ## Get only those that were not considered already
    previous = previous %>% subset(to.be.considered & !already.considered)
    ## Get the entrez for these genes
    previous = previous %>% left_join(geneInfo%>%select(entrez, symbol), by=c("Gene.symbol"="symbol"))

    for(i in 1:nrow(previous)){
        en = previous$entrez[i]
        gene_type= previous$Type[i]
        symbol = previous$Gene.symbol[i]

        if(gene_type=="SVs"){
            dd = training_ns %>% subset(entrez==en) %>% subset(BND==1 | INV==1 | INS==1 | CNVLoss==1 | CNVGain==1)
        }else if(gene_type=="CNVs"){
            dd = training_ns %>% subset(entrez==en) %>% subset(CNVLoss==1 | CNVGain==1)
        }
        cgc_drivers = rbind(cgc_drivers, dd)
    }


    load(paste0(path, "/validation_set_noScale.Rdata"))
    validation_ns = validation_ns%>%tibble::rownames_to_column()%>%
        separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")%>%
        mutate(type="P")%>%mutate(entrez=as.numeric(entrez))%>%data.frame()
    cohort = rbind(cgc_drivers, validation_ns) ## Instead of the whole training, we use only the subset of the refined CGC
    ## Get scores
    load(paste0(path, "/",scores_fn))
    scores = scores[["scores"]]

    cohort = cohort %>% left_join(scores)
    cohort$gene_type[is.na(cohort$gene_type)] = "cgc"
    cohort = cohort %>% left_join(geneInfoNCG5)

    ## Make NAs in CGC score to be the highest score in each patient + 1 so we can sort by score in the output table
    sample2maxscore = cohort %>% group_by(sample) %>% summarise(max_score=max(score, na.rm = T)+1)
    cohort$score[is.na(cohort$score) & cohort$gene_type=="cgc"] = sample2maxscore$max_score[match(cohort$sample, sample2maxscore$sample)][which(is.na(cohort$score))]


    samples = unique(cohort$sample)
    sys_cans = NULL
    if (!is.null(tissue)){
        message("Expression is considered for the selection of sys cans...")
        ## Load expression data (GTeX)
        message("Loading GTeX expression data")
        load("~/Mountpoints/rosalind_lustre/mourikisa/data/expressionNCG5_gtex.Rdata") ## name of data frame is gtex
        gtex = gtex %>% select(entrez, tissue, exp_level) %>% spread(tissue, exp_level)
        ## Checkpoint
        if(!(tissue%in%colnames(gtex))){
            stop("Wrong tissue supplied")
        }
        ## Select tissue of the cancer type
        gtex = gtex[,c("entrez", tissue)]
        colnames(gtex) = c("entrez", "gtex_tissue_expression")
        ## bring in the cohort table the expression
        message("Getting GTEx expression")
        cohort = cohort %>% left_join(gtex)


        total = length(samples)
        pb <- txtProgressBar(min = 0, max = total, style = 3)
        for (i in 1:length(samples)){
            Sys.sleep(0.1)
            ## When you get the candidates here you need to take into account expression in normal tissues and in the cancer sample

            ## Here I get the ranks for the non-CGC genes. The CGCs will be added as patient.rank==NA
            ## I also add 2 ranks, one with all the genes and one without amplifications
            d = cohort %>% subset(sample==samples[i] & gene_type!="cgc" & !gtex_tissue_expression%in%exclude_expr & !is.na(gtex_tissue_expression)) %>%
                arrange(desc(score)) %>% mutate(patient.rank=row_number(-score))
            d_withoutAmp = d %>% subset(sample==samples[i] & gene_type!="cgc" & CNVGain==0 & !gtex_tissue_expression%in%exclude_expr & !is.na(gtex_tissue_expression)) %>%
                arrange(desc(score)) %>% mutate(patient.rank.no.amp=row_number(-score))
            ## Bring in the rank without amplification
            d = d %>% left_join(d_withoutAmp)

            ## Add the CGCs
            d_cgcs = cohort %>% subset(sample==samples[i] & gene_type=="cgc" & !gtex_tissue_expression%in%exclude_expr & !is.na(gtex_tissue_expression))
            d = rbind.fill(d, d_cgcs)

            ## Add all those that are not expressed (I will exclude those, I won't use them anyways) - code left here for testing
            #d_ne = cohort %>% subset(sample==samples[i] & (gtex_tissue_expression=="Not Expressed" | is.na(gtex_tissue_expression)))
            #d = rbind.fill(d, d_ne)

            ## Add them to the big table
            sys_cans = rbind(sys_cans, d)
            setTxtProgressBar(pb, i)
        }
    }

    ## Add false positive annotation
    false_positive_genes <- read.table(fp_dir, header = T, sep = "\t")
    false_positive_genes <- false_positive_genes %>% select(entrez) %>% .$entrez
    sys_cans = sys_cans %>% mutate(NCG_FP=ifelse(entrez%in%false_positive_genes, TRUE, FALSE))

    ## In the final output we need everything (also the CGC that were excluded) so create an export table -- No we don't, I disabled this. If needed in the future maybe turn it back on
    #export = rbind(training_ns, validation_ns)
    #export = export %>% left_join(cohort)
    ## all the rest CGC are now coming in and their gene type will be NA
    #export$gene_type[is.na(export$gene_type)] = "CGC_discarded"


    message(paste0("Score file used: ", scores_fn))

    return(sys_cans)

}


## This function takes in the prediction scores and maps the first 100
## genes to pathways
scores2pathways = function(path, geneset="/Volumes/mourikisa/novel_driver_prediction/pancancer_model_training/geneSets/curated_gene_sets/c2.all.v5.1.entrez.gmt.txt"){

    load(paste0(path, "/scores.Rdata"))
    d = scores[["scores"]] %>% arrange(desc(score)) %>% slice(1:100)

    pathways = make_pathway_list(geneset)
    message(paste0(length(pathways), " pathways loaded..."))
    message(paste0(unlist(pathways) %>% unique() %>% length(), " genes found in those pathways..."))
    len = sapply(pathways, length)
    message("Excluding pathways with \"CANCER\" in their name...")
    pathways = pathways[!grepl("CANCER",names(pathways))]

    genes2map = d %>% select(entrez) %>% unique
    genes2map$path = sapply(as.list(as.character(unique(d$entrez))), function(x,y){
        sel=sapply(y, function(z,w) w%in%z, w = x);
        return( paste(names(y)[sel],collapse='.')) },
        y = pathways)

}

## This function was created to add features in the existing tables of MySQL database
## It pulls the tables from the MySQL database and adds features
## For the moment run it manually in Athena
addFeatures = function(cancer_type="OAC"){

    source("/home/mourikisa/novel_driver_prediction/pancancer_model_training/automation/config.R")

#     ## First feature to be added is the "all somatic SNVs and indels for every gene in the coresponding sample"
#     ## For this I need to load the mutation datasets from the Score_correction directory
#     load("/home/mourikisa/novel_driver_prediction/Score_correction/mutation_datasets.Rdata")
#     nsi2samples2genes = total_muts_nsi %>% subset(!is.na(entrez_19014)) %>% group_by(Cancer_type, Sample, entrez_19014, symbol_19014) %>% summarise(nsi_muts=n()) %>% ungroup %>% data.frame()
#     si2samples2genes = total_muts_si %>% subset(!is.na(entrez_19014)) %>% group_by(Cancer_type, Sample, entrez_19014, symbol_19014) %>% summarise(si_muts=n()) %>% ungroup %>% data.frame()
#     muts2samples2genes = nsi2samples2genes %>% full_join(si2samples2genes)
#     muts2samples2genes$nsi_muts[is.na(muts2samples2genes$nsi_muts)] = 0
#     muts2samples2genes$si_muts[is.na(muts2samples2genes$si_muts)] = 0
#     muts2samples2genes = muts2samples2genes %>% mutate(no_ALL_muts=nsi_muts+si_muts)
#     muts2samples2genes = muts2samples2genes %>% rename(Entrez=entrez_19014) %>% select(-symbol_19014, -nsi_muts, -si_muts)
#     muts2samples2genes$Entrez = as.numeric(muts2samples2genes$Entrez)
#     ## OAC is not in there (Add OAC)
#     load("/Volumes/mourikisa/data/OAC/Rdata/dataset_mutations_damaging.Rdata")
#     ## Watch out for this definition - depending on the version of ANNOVAR the categories may change name
#     ns = c("nonsynonymous","stopgain","frameshift deletion","splicing","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","stoploss","frameshift substitution")
#     oac_nsi2gene2sample = df_mut %>% subset(ExonicFunc.refGene%in%ns)  %>% subset(!is.na(entrez_19014)) %>% group_by(entrez_19014, symbol_19014, sample) %>% summarise(nsi_muts=n()) %>% ungroup %>%
#         mutate(Cancer_type="OAC") %>% rename(Sample=sample, Entrez=entrez_19014) %>% select(Cancer_type, Sample, Entrez, nsi_muts)
#     oac_si2gene2sample = df_mut %>% subset(!(ExonicFunc.refGene%in%ns))  %>% subset(!is.na(entrez_19014)) %>% group_by(entrez_19014, symbol_19014, sample) %>% summarise(si_muts=n()) %>% ungroup %>%
#         mutate(Cancer_type="OAC") %>% rename(Sample=sample, Entrez=entrez_19014) %>% select(Cancer_type, Sample, Entrez, si_muts)
#     oac_muts2samples2genes = oac_nsi2gene2sample %>% full_join(oac_si2gene2sample)
#     oac_muts2samples2genes$nsi_muts[is.na(oac_muts2samples2genes$nsi_muts)] = 0
#     oac_muts2samples2genes$si_muts[is.na(oac_muts2samples2genes$si_muts)] = 0
#     oac_muts2samples2genes = oac_muts2samples2genes %>% mutate(no_ALL_muts=nsi_muts+si_muts)
#     oac_muts2samples2genes = oac_muts2samples2genes %>% select(-nsi_muts, -si_muts)
#
#     d = getTable(tb=cancer_type)
#     d = d %>% left_join(muts2samples2genes)
#     d$no_ALL_muts[is.na(d$no_ALL_muts)] = 0
#     ## Reorder columns
#     d = d[,c(1:3, length(d), 4:(length(d)-1))]


#     ## Add histone compartment covariate in the features
#     gene_cov = read.table("/Volumes/mourikisa/novel_driver_prediction/Score_correction/gene_covariates.txt", header = T) %>% rename(symbol=gene) %>% select(symbol, reptime, hic)
#     ## Get their entrez through geneInfo table in NCG
#     dev_mode()
#     library(ncglib)
#     geneInfo = get_geneinfo(version="NCG5")
#     dev_mode()
#     geneInfo = geneInfo %>% select(entrez, symbol, cancer_type) %>% unique %>% rename(gene_type=cancer_type)
#     gene_cov = gene_cov %>% left_join(geneInfo)
#     ## exclude those for which we dont have entrez
#     gene_cov = gene_cov %>% subset(!is.na(entrez))
#     gene_cov = gene_cov %>% rename(Entrez=entrez) %>% select(Entrez, reptime, hic)
#     write.table(gene_cov, file="/Volumes/mourikisa/novel_driver_prediction/Score_correction/gene_covariates_geneInfo_crossed.tsv", quote = F, row.names = F, sep=)


    gene_cov = read.table(file="/home/mourikisa/novel_driver_prediction/Score_correction/gene_covariates_19014_crossed.tsv", header = T)
    ## Delete the columns if already exist
    d = d %>% select(-hic, -reptime)
    d = d %>% left_join(gene_cov)
    ## For those I have no information for hic, I set it to 0 for LAML and OAC
    ## or the rest I will set it to the median of the distribution - look at the Score_correction directory

    d$hic[is.na(d$hic)] = 25 ## Set it to the median of the overall distribution
    d$reptime[is.na(d$reptime)] = 394 ## Set it to the median of the overall distribution

    ## Correct the mistake with primates
    load("/mnt/lustre/users/k1469280/mourikisa/data/geneProperties.Rdata")
    geneProperties = geneProperties %>% select(Entrez, origin) %>% unique %>% subset(origin=="Primates") %>%
      mutate(primates=1) %>% select(-origin)
    d = d %>% left_join(geneProperties)
    d$primates[is.na(d$primates)] = 0
    ## Make primates factor
    d$primates = as.factor(d$primates)

    ## Reorder columns
    d = d[,c(1:45, length(d), 46:(length(d)-1))]


    ## Add microRNAs
    ## load ncglib
#     con = get_ncg_connection()
#     con %>% tbl_df()
#     mirna = con %>% tbl("NCG51_mirna_interactions") %>% collect(n=Inf)
    load("/mnt/lustre/users/k1469280/mourikisa/data/miRNA_NCG51.Rdata")
    ## Select only the human miRNAs
    mirna = mirna %>% separate(mir, into=c("sp", "crap"), sep="-", extra="merge", remove = F) %>% select(-crap) %>% subset(sp=="hsa")
    mirna = mirna %>% group_by(entrez) %>% summarise(mirna=length(unique(mir))) %>% rename(Entrez=entrez) %>% data.frame()
    ## For testing purposes (to see if there is actually a difference)
    #load("/Volumes/mourikisa/data/geneInfoNCG5.Rdata")
    #mirna = mirna %>% left_join(geneInfo%>%select(entrez, symbol, cancer_type)%>%rename(gene_type=cancer_type), by=c("Entrez"="entrez"))
    d = d %>% left_join(mirna)
    d$mirna[is.na(d$mirna)] = 0


    write.table(d, file="/home/mourikisa/novel_driver_prediction/Rdata/total_table_TCGA_19014_prepared_LAML.tsv", quote = F, sep = "\t", row.names = F, col.names = F)

#     ## then in the MySQLdatabase I run the following
#     rename table LAML to LAML_test;
#     create table PAAD like PAAD_test;
#     alter table PAAD add no_ALL_muts INT(11) after Entrez;
#     alter table PAAD add primates tinyint(1) after mammals;
#     alter table PAAD add reptime INT(11);
#     alter table PAAD add hic INT(11);

#     load data local infile '/home/mourikisa/novel_driver_prediction/Rdata/total_table_TCGA_19014_prepared_PAAD.tsv' into table PAAD;

    ## ------------------------------------------------------------------------
    ## Update WGD using the raw data - the results are always saved in the geneProperties_final.Rdata
    ## WGD data (raw)
    wgd = read.xlsx2("/Users/thmourikis/Desktop/WGD_data.xlsx", 1, startRow = 2)
    ## get geneSymbols
    load("/Users/thmourikis/Desktop/19014_GeneSymbols.Rdata")
    noWGD_19014 = geneSymbols %>% subset(symbol%in%wgd$Gene[wgd$Ohnolog..O.!="O"]) %>%
        select(Entrez,Symbol) %>% unique %>% mutate(WGD=0)
    WGD_19014 = geneSymbols %>% subset(symbol%in%wgd$Gene[wgd$Ohnolog..O.=="O"]) %>%
        select(Entrez,Symbol) %>% unique %>% mutate(WGD=1)
    WGD_data_19014 = rbind(WGD_19014, noWGD_19014) ## It may be the case that there are aliases as Ohnologs and non Ohnologs
    WGD_data_19014 = WGD_data_19014 %>%  arrange(desc(WGD)) %>% subset(!duplicated(Entrez))


    ## Update the Duplicability (I am adding the highest of the two - genic or Genomic because this is what Giovanni found significant)
    ## Get the alldups table with all duplicability information
    alldups = read.table('/Users/fc-8s-imac/ncg-data-update/duplicability_CGC_analysis/output/duplicability_withdomains.csv', header=T) %>%
        subset(mainEntrez==T) %>%
        rename(genegroup=variable)

    ## Get also the geneProperties_complete
    load("/Volumes/mourikisa/data/geneProperties_complete.Rdata")
    ## Remember to take unique in the joining because alldups is redundant for dataset column
    geneProperties = geneProperties %>% rename(Genic_old=Genic) %>% left_join(alldups%>%select(Entrez, Duplicability, Genic, Genomic)%>%unique)
    ## Get duplicability in 60% because this where Giovanni saw difference
    geneProperties = geneProperties %>% mutate(duplicated=ifelse(Duplicability>=60, 1, 0))
    save(geneProperties, file="/Volumes/mourikisa/data/geneProperties_complete.Rdata")

    ## Create geneProperties final
    ## Fic domains cause Giovanni set it to 0
    geneProperties = geneProperties %>% mutate(alldomains=ifelse(inCDD==0, NA, alldomains))
    geneProperties = geneProperties %>% rename(WGD=WGD_new)
    geneProperties = geneProperties %>% select(Entrez, Length.fullrefseq, duplicated, WGD, alldomains,
                              degree, betweenness, hub, central, age, origin, exp.breadth,
                              tot.tissues, hic, mirna)
    save(geneProperties, file="/Volumes/mourikisa/data/geneProperties_final.Rdata")


}


## Prepare tables for each cancer type and put them in the MySQL database
# for (ct in cancer_types){
#     cat(ct, "\n")
#     df = prepareTableML(df=total_table, cancer_type=ct)
#     write.table(df, file=paste("/home/mourikisa/novel_driver_prediction/Rdata/total_table_TCGA_19014_prepared_", ct, ".tsv", sep=""),
#                 row.names = F, col.names = F, sep = "\t", quote = F)
# }
## Load all data from 129_OAC/Rdata/total_table_19014.Rdata
getMLinput <- function(df, geneProperties_dir="/mnt/lustre/users/k1469280/mourikisa/data/geneProperties_final_mmImputed.Rdata"){

    df = df %>% mutate(Cancer_type="OAC") %>% rename(Sample=sample, Entrez=entrez_19014, Copy_number=Total_CN, CNV_type=CNV_type_corrected) %>% select(-symbol_19014, -ploidy, -CNV_entries, -in_19014)

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
## I just run in Athena createAllinput() after sourcing config.R
createAllInput <- function(save_dir="/home/mourikisa/novel_driver_prediction/pancancer_model_training/automation/median_mode_imputation/input_data"){
    cancer_types <- c("ACC", "BLCA", "CESC", "CHOL",
                      "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD",
                      "LUSC", "OV", "PRAD", "READ", "SARC", "SKCM",
                      "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
    message("Loading pancancer table....")
    load("/home/mourikisa/novel_driver_prediction/Rdata/total_table_TCGA_19014.Rdata")
    for (cancer_type in cancer_types){
        message(cancer_type)
        d = getMLinput(ct=cancer_type, df=total_table)
        write.table(d, file=paste0(save_dir, "/", cancer_type, "_ML.tsv"), row.names = F, quote = F, sep = "\t")
    }
}

doFeatureSelection = function(training, columnRange=1:length(training)){
    ## Based on http://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
    ## 1st way
    # The Caret R package provides the findCorrelation which will analyze a correlation matrix of your data’s attributes report on attributes that can be removed.
    # A correlation matrix is created from these attributes and highly correlated attributes are identified, in this case the age attribute is remove as it correlates highly with the pregnant attribute.
    # Generally, you want to remove attributes with an absolute correlation of 0.75 or higher.

    # ensure the results are repeatable
    set.seed(7)
    # load the library
    library(mlbench)
    library(caret)

    # calculate correlation matrix
    num_training = sapply(training%>%select(-type)%>%data.frame(), function(x) as.numeric(x))
    num_training = num_training[,apply(num_training, 2, var, na.rm=TRUE) != 0]
    correlationMatrix <- cor(num_training)
    # summarize the correlation matrix
    print(correlationMatrix)
    # find attributes that are highly corrected (ideally >0.75)
    highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75, names = T)
    # print indexes of highly correlated attributes
    print(highlyCorrelated)


}


## Functions to get the mutations from TCGA data
applyMutationFilters = function(df=NULL){
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
## Watch out! Some of the terms below appear in the ExonicFunc.refGene
## and some in the Func.refGene of Annovar
## Make sure you include all the mutations
exonicFunc_nonsilent = c("frameshift deletion","frameshift insertion",
                         "frameshift substitution","stopgain","stoploss",
                         "nonsynonymous SNV","splicing")
exonicFunc_silent = c("nonframeshift deletion","nonframeshift insertion",
                      "nonframeshift substitution", "synonymous SNV")
## We also need to take intronic mutations
non_coding = c("intronic","UTR3", "UTR5")
nonsilentSNVs = c("stopgain","stoploss",
                  "nonsynonymous SNV","splicing")

getSilentMutations = function(dataDir="/Volumes/mourikisa/data/TCGA/01_03_2015/",
                        resDir="/Volumes/mourikisa/novel_driver_prediction/Score_correction/",
                        cancer_types=NULL){

    total_muts = NULL
    for (cancer_type in cancer_types){
        cat(cancer_type, "\n")
        tumour_path = paste0(dataDir, cancer_type, "/Tumor/Somatic_Mutations/", cancer_type,
                             "_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF_Annovar_dbNSFP_gAnn_SampleFilters_Silent.Rdata")
        load(tumour_path)
        maf = do.call(rbind, somatic_mutations)
        total_muts = rbind(total_muts, maf)
    }
    total_muts = unrowname(total_muts)
    return(total_muts)
}

getNonSilentMutations = function(dataDir="~/athena/data/TCGA/01_03_2015/",
                                  cancer_types=NULL){

    total_muts = NULL
    for (cancer_type in cancer_types){
        cat(cancer_type, "\n")
        tumour_path = paste0(dataDir, cancer_type, "/Tumor/Somatic_Mutations/", cancer_type,
                             "_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF_Annovar_dbNSFP_gAnn_SampleFilters_NonSilent_MutationFilters_OncodriveClust.Rdata")
        load(tumour_path)
        maf = do.call(rbind, somatic_mutations)
        total_muts = rbind(total_muts, maf)
    }
    total_muts = unrowname(total_muts)
    return(total_muts)
}

# ## Load mutation datasets in Score_correction directory
# par(mfrow = c(2, 2))
# with(test,scatter.smooth(x=driver_rate, y=muts_rate, lpars=list(col = "red", lwd = 3, lty = 3),
#                          col=ifelse(cancer_type=="cgc","red", ifelse(cancer_type=="can", "blue",ifelse(cancer_type=="rst","green", "black"))), pch=16,
#                          main="All genes"))
# legend(x="topleft",c("cgc", "can", "rst"),pch=16,col=c("red","blue","green"))
# with(test, text(x=driver_rate, y=muts_rate-0.03, labels=ifelse(driver_rate>0.3 | muts_rate>0.4, symbol, ""), cex = 0.8))
# text(x=0.85, y=0, paste0("cor: ", format(cor.test(test$muts_rate, test$driver_rate)$estimate, digits = 3)))
#
# with(test%>%subset(cancer_type=="cgc"),scatter.smooth(x=driver_rate, y=muts_rate, lpars=list(col = "red", lwd = 3, lty = 3),
#                                                       col=ifelse(cancer_type=="cgc","red", "black"), pch=16,
#                                                       main="CGC", xlim=c(0,0.1), ylim=c(0,0.4)))
# with(test%>%subset(cancer_type=="cgc"), text(x=driver_rate, y=muts_rate-0.03, labels=ifelse(driver_rate>0.05, symbol, ""), cex = 0.8))
# text(x=0.09, y=0, paste0("cor: ", format(cor.test(test$muts_rate[test$cancer_type=="cgc"], test$driver_rate[test$cancer_type=="cgc"])$estimate, digits = 3)))
#
# with(test%>%subset(cancer_type=="can"),scatter.smooth(x=driver_rate, y=muts_rate, lpars=list(col = "red", lwd = 3, lty = 3),
#                                                       col=ifelse(cancer_type=="can","blue", "black"), pch=16,
#                                                       main="Candidates", xlim=c(0,0.1), ylim=c(0,0.4)))
# with(test%>%subset(cancer_type=="can"), text(x=driver_rate, y=muts_rate-0.03, labels=ifelse(driver_rate>0.08 | muts_rate>0.2, symbol, ""), cex = 0.8))
# text(x=0.09, y=0, paste0("cor: ", format(cor.test(test$muts_rate[test$cancer_type=="can"], test$driver_rate[test$cancer_type=="can"])$estimate, digits = 3)))
#
# with(test%>%subset(cancer_type=="rst"),scatter.smooth(x=driver_rate, y=muts_rate, lpars=list(col = "red", lwd = 3, lty = 3),
#                                                       col=ifelse(cancer_type=="rst","green", "black"), pch=16,
#                                                       main="Rest", xlim=c(0,0.1), ylim=c(0,0.4)))
# with(test%>%subset(cancer_type=="rst"), text(x=driver_rate, y=muts_rate-0.03, labels=ifelse(driver_rate>0.07 | muts_rate>0.2, symbol, ""), cex = 0.8))
# text(x=0.09, y=0, paste0("cor: ", format(cor.test(test$muts_rate[test$cancer_type=="rst"], test$driver_rate[test$cancer_type=="rst"])$estimate, digits = 3)))

# ## Plot the synonymous VS the driver rate
# par(mfrow = c(2, 2))
# with(test,scatter.smooth(x=driver_rate, y=syn_muts_rate, lpars=list(col = "red", lwd = 3, lty = 3),
#                          col=ifelse(cancer_type=="cgc","red", ifelse(cancer_type=="can", "blue",ifelse(cancer_type=="rst","green", "black"))), pch=16,
#                          main="All genes"))
# legend(x="topleft",c("cgc", "can", "rst"),pch=16,col=c("red","blue","green"))
# with(test, text(x=driver_rate, y=syn_muts_rate-0.01, labels=ifelse(driver_rate>0.2 | syn_muts_rate>0.10, symbol, ""), cex = 0.8))
# text(x=0.85, y=0.25, paste0("cor: ", format(cor.test(test$syn_muts_rate, test$driver_rate)$estimate, digits = 3)))
#
# with(test%>%subset(cancer_type=="cgc"),scatter.smooth(x=driver_rate, y=syn_muts_rate, lpars=list(col = "red", lwd = 3, lty = 3),
#                                                       col=ifelse(cancer_type=="cgc","red", "black"), pch=16,
#                                                       main="CGC", xlim=c(0,0.2), ylim=c(0,0.15)))
# with(test%>%subset(cancer_type=="cgc"), text(x=driver_rate, y=syn_muts_rate-0.005, labels=ifelse(driver_rate>0.1 | syn_muts_rate>0.02, symbol, ""), cex = 0.8))
# text(x=0.175, y=0.15, paste0("cor: ", format(cor.test(test$syn_muts_rate[test$cancer_type=="cgc"], test$driver_rate[test$cancer_type=="cgc"])$estimate, digits = 3)))
#
#
# with(test%>%subset(cancer_type=="can"),scatter.smooth(x=driver_rate, y=syn_muts_rate, lpars=list(col = "red", lwd = 3, lty = 3),
#                                                       col=ifelse(cancer_type=="can","blue", "black"), pch=16,
#                                                       main="Candidates", xlim=c(0,0.2), ylim=c(0,0.15)))
# with(test%>%subset(cancer_type=="can"), text(x=driver_rate, y=syn_muts_rate-0.005, labels=ifelse(driver_rate>0.075 | syn_muts_rate>0.075, symbol, ""), cex = 0.8))
# text(x=0.175, y=0.15, paste0("cor: ", format(cor.test(test$syn_muts_rate[test$cancer_type=="can"], test$driver_rate[test$cancer_type=="can"])$estimate, digits = 3)))
#
#
# with(test%>%subset(cancer_type=="rst"),scatter.smooth(x=driver_rate, y=syn_muts_rate, lpars=list(col = "red", lwd = 3, lty = 3),
#                                                       col=ifelse(cancer_type=="rst","green", "black"), pch=16,
#                                                       main="Candidates", xlim=c(0,0.2), ylim=c(0,0.15)))
# with(test%>%subset(cancer_type=="rst"), text(x=driver_rate, y=syn_muts_rate-0.005, labels=ifelse(driver_rate>0.08 | syn_muts_rate>0.075, symbol, ""), cex = 0.8))
# text(x=0.175, y=0.15, paste0("cor: ", format(cor.test(test$syn_muts_rate[test$cancer_type=="rst"], test$driver_rate[test$cancer_type=="rst"])$estimate, digits = 3)))

## Get the covariates
# cov = read.table("/Volumes/mourikisa/novel_driver_prediction/Score_correction/gene_covariates.txt", header = T)
# pdf("/Volumes/mourikisa/novel_driver_prediction/Score_correction/replication_time_distribution.pdf")
# print(hist(cov$reptime, breaks = 100, main = "Replication time (1134 NAs excluded)", xlab = "Reptime"))
# dev.off()
# pdf("/Volumes/mourikisa/novel_driver_prediction/Score_correction/hic_distribution.pdf")
# print(hist(cov$hic, breaks = 100, main = "Chromatin compartment (414 NAs excluded)", xlab = "HiC"))
# dev.off()

# ## Check the values of covariates for CGC
# cov = cov %>% left_join(geneInfo, by=c("gene"="symbol")) ## There are about 10 duplicated symbols - not much of a problem at this stage
# par(mfrow = c(2, 2))
# with(cov, hist(reptime, breaks=100, main="All genes"))
# abline(v=100, lty=2, col="red")
# abline(v=1000, lty=2, col="red")
# with(cov,legend(x="topright",legend=cbind(paste0(names(summary(reptime)), ":", summary(reptime)))))
#
# with(cov%>%subset(cancer_type=="cgc"), hist(reptime, breaks=100, main="CGC", col = "red", xlim = c(0,1700)))
# abline(v=100, lty=2, col="red")
# abline(v=1000, lty=2, col="red")
# with(cov%>%subset(cancer_type=="cgc"),legend(x="topright",legend=cbind(paste0(names(summary(reptime)), ":", summary(reptime)))))
#
# with(cov%>%subset(cancer_type=="can"), hist(reptime, breaks=100, main="Candidates", col="blue", xlim = c(0,1700)))
# abline(v=100, lty=2, col="red")
# abline(v=1000, lty=2, col="red")
# with(cov%>%subset(cancer_type=="can"),legend(x="topright",legend=cbind(paste0(names(summary(reptime)), ":", summary(reptime)))))
#
# with(cov%>%subset(cancer_type=="rst"), hist(reptime, breaks=100, main="Rest", col="green", xlim = c(0,1700)))
# abline(v=100, lty=2, col="red")
# abline(v=1000, lty=2, col="red")
# with(cov%>%subset(cancer_type=="rst"),legend(x="topright",legend=cbind(paste0(names(summary(reptime)), ":", summary(reptime)))))
#
# par(mfrow = c(2, 2))
# with(cov, hist(hic, breaks=100, main="All genes"))
# abline(v=-50, lty=2, col="red")
# abline(v=50, lty=2, col="red")
# with(cov,legend(x="topleft",legend=cbind(paste0(names(summary(hic)), ":", summary(hic)))))
#
# with(cov%>%subset(cancer_type=="cgc"), hist(hic, breaks=100, main="CGC", col="red", xlim = c(-100,100)))
# abline(v=-50, lty=2, col="red")
# abline(v=50, lty=2, col="red")
# with(cov%>%subset(cancer_type=="cgc"),legend(x="topleft",legend=cbind(paste0(names(summary(hic)), ":", summary(hic)))))
#
# with(cov%>%subset(cancer_type=="can"), hist(hic, breaks=100, main="Candidates", col="blue", xlim = c(-100,100)))
# abline(v=-50, lty=2, col="red")
# abline(v=50, lty=2, col="red")
# with(cov%>%subset(cancer_type=="can"),legend(x="topleft",legend=cbind(paste0(names(summary(hic)), ":", summary(hic)))))
#
# with(cov%>%subset(cancer_type=="rst"), hist(hic, breaks=100, main="Rest", col="green", xlim = c(-100,100)))
# abline(v=-50, lty=2, col="red")
# abline(v=50, lty=2, col="red")
# with(cov%>%subset(cancer_type=="rst"),legend(x="topleft",legend=cbind(paste0(names(summary(hic)), ":", summary(hic)))))

## Wilcoxon on the distributions of cgc and rest
# gene_cov = read.table("/Volumes/mourikisa/novel_driver_prediction/Score_correction/gene_covariates.txt", header = T) %>% rename(symbol=gene) %>% select(symbol, reptime, hic)
# ## Get their entrez through geneInfo table in NCG
# dev_mode()
# library(ncglib)
# geneInfo = get_geneinfo(version="NCG5")
# dev_mode()
# geneInfo = geneInfo %>% select(entrez, symbol, cancer_type) %>% unique %>% rename(gene_type=cancer_type)
# gene_cov = gene_cov %>% left_join(geneInfo)
# ## exclude those for which we dont have entrez
# gene_cov = gene_cov %>% subset(!is.na(entrez))
# gene_cov = gene_cov %>% rename(Entrez=entrez) %>% select(Entrez, reptime, hic)
# wilcox.test(gene_cov%>%subset(gene_type=="cgc")%>%.$hic, gene_cov%>%subset(gene_type=="rst")%>%.$hic)
# wilcox.test(gene_cov%>%subset(gene_type=="cgc")%>%.$reptime, gene_cov%>%subset(gene_type=="rst")%>%.$reptime)

getNCGannotation = function(){
#     dev_mode()
#     library(ncglib)
#     ncg_connection = get_ncg_connection()
#     geneInfo = get_geneinfo("NCG5")
#     expValidation = ncg_connection %>% tbl("NCG5_experimental_validation") %>% select(entrez, methods) %>% data.frame()
#     cancerGenes = get_cancer_genes(version="NCG5")
#     dupbythresholds = get_dupbythresholds()
#     expression = ncg_connection %>% tbl("NCG5_expression_GTEx2015") %>% collect(n=Inf)
#     dev_mode()

    load("/mnt/lustre/users/k1469280/mourikisa/data/geneInfoNCG5.Rdata")
    load("/mnt/lustre/users/k1469280/mourikisa/data/cancerGenesNCG5.Rdata")
    load("/mnt/lustre/users/k1469280/mourikisa/data/expValidationNCG5.Rdata")
    load("/mnt/lustre/users/k1469280/mourikisa/data/dupByThresholdsNCG5.Rdata")
    load("/mnt/lustre/users/k1469280/mourikisa/data/expressionNCG5_gtex.Rdata")

    ## Fix gene info table from NCG
    geneInfo = geneInfo %>% select(symbol, entrez, description, uniprot, cancer, cancer_type, cancer_dom, cancer_rec, duplicability, origin, degree, complex) %>% unique
    ## Get a cancer gene with all the associated primary sites and cancer sites
    cancerGenes = cancerGenes %>% select(entrez, cancer_site) %>%
        group_by(entrez) %>% summarise(cancer_site=paste(unique(cancer_site), collapse=",")) %>%
        ungroup

    geneInfo = geneInfo %>% left_join(cancerGenes) %>% rename(cancer_info=cancer, gene_type=cancer_type)

    genicDup = dupbythresholds %>% subset(threshold=="60" | threshold=="20") %>% data.frame() %>% spread(threshold, ndups)
    colnames(genicDup) = c("entrez", "genic20", "genic60")
    geneInfo = geneInfo %>% left_join(genicDup)

    return(list(geneInfo=geneInfo, gtex=gtex))

}

## Variance and Dimensionality reduction
getBestmodel = function( test){
    ## without young plus All_muts plus HIC
    ## Change to the directory
    RES_DIR_LOCAL = "/Volumes/mourikisa/novel_driver_prediction/OAC_without_young_plusALLmuts_plusHIC/cgc_518_training"
    cv_stat_summary = read.table(paste0(RES_DIR_LOCAL, "/cv_stat_summary.tsv"), header=T)
    setwd(paste0(RES_DIR_LOCAL, "/excels"))
    best_models = NULL

    linear = read.xlsx(file="old/grid_predictions_linear.xlsx", 1, startRow = 7)
    linear = linear %>% select(-min, -q1, -median, -mean, -q3, -max)
    linear = linear %>% left_join(cv_stat_summary %>% subset(type=="Sensitivity" & set=="test") %>% select(analysis, iterations, min, q1, median, mean, q3, max, var)) %>% data.frame()
    linear = linear%>%rename(Predicted_genes=novel_drivers_genes_platt_rst_0.95)
    ## Print the best model
    best_models = rbind(best_models,
                        linear %>% arrange(desc(mean)) %>%
                            slice(1:5) %>% mutate(var_rank=row_number(var)) %>%
                            subset(var_rank==1) %>% select(analysis, nu, gamma, mean, var, Predicted_genes)
                        )
    p1= ggplot(linear, aes(x=var, y=mean, color=Predicted_genes)) +
        geom_point() + geom_point(data=linear%>%subset(analysis=="linear.0.2"), size=4) +
        geom_point(data=linear%>%subset(analysis=="linear.0.1"), size=5) +
        ggtitle("Linear models") +
        scale_color_gradient(low="blue", high="red") +
        xlab("Sensitivity variance") + ylab("Mean sensitivity")

    radial = read.xlsx(file="old/grid_predictions_radial.xlsx", 1, startRow = 7)
    radial = radial %>% select(-min, -q1, -median, -mean, -q3, -max)
    radial = radial %>% left_join(cv_stat_summary %>% subset(type=="Sensitivity" & set=="test") %>% select(analysis, iterations, min, q1, median, mean, q3, max, var)) %>% data.frame()
    radial = radial%>%rename(Predicted_genes=novel_drivers_genes_platt_rst_0.95)
    ## Print the best model
    best_models = rbind(best_models,
                        radial %>% arrange(desc(mean)) %>%
                            slice(1:5) %>% mutate(var_rank=row_number(var)) %>%
                            subset(var_rank==1) %>% select(analysis, nu, gamma, mean, var, Predicted_genes)
                        )
    p2= ggplot(radial, aes(x=var, y=mean, color=Predicted_genes)) + geom_point() + geom_point(data=radial%>%subset(analysis=="radial.0.2.0.015625"), size=4) +
        geom_point(data=radial%>%subset(analysis=="radial.0.1.0.0078125"), size=5) +
        ggtitle("Radial models") +
        scale_color_gradient(low="blue", high="red") +
        xlab("Sensitivity variance") + ylab("Mean sensitivity")

    polynomial = read.xlsx(file="old/grid_predictions_polynomial.xlsx", 1, startRow = 7)
    polynomial = polynomial %>% select(-min, -q1, -median, -mean, -q3, -max)
    polynomial = polynomial %>% left_join(cv_stat_summary %>% subset(type=="Sensitivity" & set=="test") %>% select(analysis, iterations, min, q1, median, mean, q3, max, var)) %>% data.frame()
    polynomial=polynomial%>%rename(Predicted_genes=novel_drivers_genes_platt_rst_0.95)
    ## Print the best model
    best_models = rbind(best_models,
                        polynomial %>% arrange(desc(mean)) %>%
                            slice(1:5) %>% mutate(var_rank=row_number(var)) %>%
                            subset(var_rank==1) %>% select(analysis, nu, gamma, mean, var, Predicted_genes)
                        )
    p3= ggplot(polynomial, aes(x=var, y=mean, color=Predicted_genes)) + geom_point() +
        geom_point(data=polynomial%>%subset(analysis=="polynomial.0.05.0.125.3"), size=5) +
        ggtitle("Polynomial models") +
        scale_color_gradient(low="blue", high="red") +
        xlab("Sensitivity variance") + ylab("Mean sensitivity")

    sigmoid = read.xlsx(file="old/grid_predictions_sigmoid.xlsx", 1, startRow = 7)
    sigmoid = sigmoid %>% select(-min, -q1, -median, -mean, -q3, -max)
    sigmoid = sigmoid %>% left_join(cv_stat_summary %>% subset(type=="Sensitivity" & set=="test") %>% select(analysis, iterations, min, q1, median, mean, q3, max, var)) %>% data.frame()
    sigmoid=sigmoid%>%rename(Predicted_genes=novel_drivers_genes_platt_rst_0.95)
    ## Print the best model
    best_models = rbind(best_models,
                        sigmoid %>% arrange(desc(mean)) %>%
                            slice(1:5) %>% mutate(var_rank=row_number(var)) %>%
                            subset(var_rank==1) %>% select(analysis, nu, gamma, mean, var, Predicted_genes)
                        )
    p4= ggplot(sigmoid, aes(x=var, y=mean, color=Predicted_genes)) + geom_point() +
        geom_point(data=sigmoid%>%subset(analysis=="sigmoid.0.05.0.0625"), size=4) +
        geom_point(data=sigmoid%>%subset(analysis=="sigmoid.0.05.0.5"), size=5) +
        ggtitle("Sigmoid models") +
        scale_color_gradient(low="blue", high="red") +
        xlab("Sensitivity variance") + ylab("Mean sensitivity")

    cat("Models picked...")
    print(best_models)

    pdf(file="/Volumes/mourikisa/novel_driver_prediction/OAC_without_young_plusALLmuts_plusHIC/cgc_518_training/1489_79_OAC_sensitivity_VS_variance.pdf", height = 20, width = 20)
    grid.arrange(p1,p2,p3,p4,ncol=2)
    dev.off()

    write.xlsx(linear%>%data.frame(), file="grid_predictions_linear.xlsx", row.names = F)
    write.xlsx(radial%>%data.frame(), file="grid_predictions_radial.xlsx", row.names = F)
    write.xlsx(polynomial%>%data.frame(), file="grid_predictions_polynomial.xlsx", row.names = F)
    write.xlsx(sigmoid%>%data.frame(), file="grid_predictions_sigmoid.xlsx", row.names = F)

    ## Get the scores (based on the best models above)
    without_young_plusALL_plusHIC_score = getScore(cancer_type="OAC", path="/Volumes/mourikisa/novel_driver_prediction/OAC_without_young_plusALLmuts_plusHIC/cgc_518_training/gridPrediction/",
                                                   linear_preds = "predictions_linear.0.1.Rdata", linearCVS = 0.61, linearVar = 0.02,
                                                   radial_preds = "predictions_radial.0.1.0.0078125.Rdata", radialCVS = 0.72, radialVar = 0.03,
                                                   polynomial_preds = "predictions_polynomial.0.05.0.125.3.Rdata", polynomialCVS = 0.54, polynomialVar = 0.02,
                                                   sigmoid_preds = "predictions_sigmoid.0.05.0.5.Rdata", sigmoidCVS = 0.93, sigmoidVar = 0.01)

    ## Get geneInfo
    annotation = getNCGannotation()
    TISSUE = "Esophagus"

    ## Get features and scores (for dimensionality reduction)
    load("/Volumes/mourikisa/novel_driver_prediction/OAC_without_young_plusALLmuts_plusHIC/cgc_518_training/training_set.Rdata")
    load("/Volumes/mourikisa/novel_driver_prediction/OAC_without_young_plusALLmuts_plusHIC/cgc_518_training/validation_set.Rdata")

    cohort = rbind(training%>%tibble::rownames_to_column()%>%separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>% data.frame(),
                   validation%>%tibble::rownames_to_column()%>%separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")%>%mutate(type="P")%>%data.frame())
    cohort = cohort %>% mutate(entrez=as.numeric(entrez)) %>% left_join(without_young_plusALL_plusHIC_score[["scores"]])
    cohort$gene_type[is.na(cohort$gene_type)] = "cgc"
    cohort = cohort %>% left_join(annotation[["gtex"]]%>%subset(tissue==TISSUE)%>%select(entrez, exp_level)%>%rename(tissue_exp_level=exp_level))

    test = cohort
    test[,4:37] <- sapply(test[,4:37], as.numeric)
    test = test %>% select(-one_of(nearZeroVar(test[,4:37], name=T)))

    ## PCA
    test.pca <- prcomp(test[1:79,c(8:18,20:31)])
    print(test.pca)
    library(ggbiplot)
    g = ggbiplot(test.pca, obs.scale = 1, var.scale = 1,
                 groups = test[1:79,"type"], ellipse = TRUE,
                 circle = FALSE)
    g <- g + scale_color_discrete(name = '')
    g <- g + theme(legend.direction = 'horizontal',
                   legend.position = 'top')
    print(g)

    ## tSNE
    # load the tsne package
    library(tsne)

    # initialize counter to 0
    x <- 0
    epc <- function(x) {
        x <<- x + 1
        filename <- paste("plot", x, "jpg", sep=".")
        cat("> Plotting TSNE to ", filename, " ")

        # plot to d:\\plot.x.jpg file of 2400x1800 dimension
        jpeg(filename, width=2400, height=1800)

        plot(x, t='n', main="T-SNE")
        text(x, labels=t$label)
        dev.off()
    }

    tsne_data <- tsne(test[,c(8:18,24:31)], max_iter=1000, epoch=100)
    test = cbind(test, tsne_data %>% data.frame())
    test = test %>% left_join(geneInfo)


    ggplot(test%>%subset(is.na(score) | score>28), aes(x=X1, y=X2, shape=type, color=score)) + geom_point(size=2) +scale_shape_manual(values = c(3, 15)) + scale_color_gradient(low = "blue", high = "red")

    p1 = ggplot(test, aes(X1, X2)) + geom_point(aes(color=score, shape=type)) +
        geom_density2d(data=test%>%subset(type=="C"),contour = T) + ## density based on training
        #geom_density2d(data=test%>%subset(type=="P"), contour=T, aes(fill = score)) +
        scale_color_gradient(low = "blue", high = "red") +
        geom_point(data=test%>%subset(type=="C"), size=2, color="black")

    p2 = ggplot(rbind(test%>%subset(type=="C"), test%>%subset(type=="P")%>%arrange(desc(score))%>%slice(1:100)), aes(X1, X2)) + geom_point(aes(color=score, shape=type)) +
        geom_density2d(data=test%>%subset(type=="C"),contour = T) + ## density based on training
        #geom_density2d(data=test%>%subset(type=="P"), contour=T, aes(fill = score)) +
        scale_color_gradient(low = "blue", high = "red") +
        geom_point(data=test%>%subset(type=="C"), size=2, color="black")

    grid.arrange(p1,p2,nrow=2)

    ## Isomap
    library(RDRToolbox)
    ## Only systems-level properties
    #iso_dim2 = Isomap(data=as.matrix(test[,c(8:18,24:31)]), dims=2)
    ## All features
    iso_dim2 = Isomap(data=as.matrix(test[,c(4:31)]), dims=2)

    iso_dim2 = iso_dim2 %>% data.frame()
    colnames(iso_dim2) = c("X1", "X2")

    test = cbind(test, iso_dim2)
    doms_recs = rbind(annotation[["geneInfo"]] %>% subset(cancer_dom==1) %>% select(entrez) %>% mutate(cancer_gene_type="Dominant"),
    annotation[["geneInfo"]] %>% subset(cancer_rec==1) %>% select(entrez) %>% mutate(cancer_gene_type="Recessive"))
    doms_recs = doms_recs %>% arrange(entrez, desc(cancer_gene_type)) %>% subset(!duplicated(entrez)) ## There are 4 genes with both Dominant and recessive annotation
    test = test %>% left_join(doms_recs)
    test = test %>% left_join(annotation[["geneInfo"]]%>%select(entrez, symbol)%>%unique)

    ## Plot the score gradient
    ggplot(test, aes(X1, X2)) + geom_point(aes(shape=type, color=score, text=paste("symbol/sample:", symbol, "/", sample))) +
        scale_color_gradient(low="blue", high="red") +
        geom_density2d(data=test%>%subset(type=="C"),contour = T) + ## density based on training
        geom_point(data=test%>%subset(type=="C"), aes(X1, X2, shape=type, text=paste("symbol/sample:", symbol, "/", sample))) +
        theme_boss() + ggtitle("Isopmap using all features")

    ggplotly()


    ## Make the plot with training and 44 only
    ## Add 44
    genes_44 = read.xlsx(file="/Volumes/mourikisa/novel_driver_prediction/OAC/cgc_518_training/shared_genes_annotation_0.95_intersection.xlsx", 1) %>%
        select(entrez) %>% unique %>% .$entrez
    test = test %>% mutate(in_44=ifelse(entrez%in%genes_44, TRUE, FALSE))
    test = test %>%
        mutate(color=ifelse(gene_type=="cgc" & (cancer_gene_type=="Recessive" | is.na(cancer_gene_type)), "Recessive",
                            ifelse(gene_type=="cgc" & (cancer_gene_type=="Dominant" & !is.na(cancer_gene_type)), "Dominant", NA)))
    test$color[is.na(test$color) & test$in_44] = "Prediction"
    test$color[is.na(test$color)] = "Other"
    p = ggplot(test, aes(X1, X2)) + geom_jitter(data=test%>%subset(color=="Other"), aes(X1, X2, color=color)) +
        geom_jitter(data=test%>%subset(color%in%c("Dominant", "Recessive")), aes(X1, X2, color=color)) +
        geom_jitter(data=test%>%subset(color=="Prediction"), aes(X1, X2, color=color)) +
        scale_colour_manual(values = c("green", "grey", "red", "black")) +
        geom_density2d(data=test%>%subset(type=="C"),contour = T) + ## density based on training
        geom_text(data=test%>%subset(symbol%in%c("RING1", "NCOR1", "SIN3A", "KAT2A", "KAT2B", "ITCH")), aes(X1, X2, label = symbol)) +
        theme_boss()
    pdf("/Volumes/FCShared/Thanos/CRUK/44_OAC_Isomap.pdf", width = 15, height = 10, useDingbats = F)
    print(p)
    dev.off()

    ## Get the minimum score for at least 3 predictions
    max_score = max(test$score[!is.na(test$score)])
    min_score_atleast3 = 10000000
    sys_cans = NULL
    for (s in unique(test$sample)){
        cat(s, "\n")
        d = test %>% subset(sample==s) %>% mutate(score=ifelse(is.na(score), max_score+1, score)) %>% arrange(sample, desc(score)) %>% slice(1:3)
        min_score = min(d$score[!is.na(d$score)])
        if (min_score < min_score_atleast3){
            min_score_atleast3 = min_score
        }
        sys_cans = rbind(sys_cans, d)
    }

    ## At least 3
    test_atleast3 = test %>%
        mutate(color=ifelse(gene_type=="cgc" & (cancer_gene_type=="Recessive" | is.na(cancer_gene_type)), "Recessive",
                            ifelse(gene_type=="cgc" & (cancer_gene_type=="Dominant" & !is.na(cancer_gene_type)), "Dominant", NA)))
    test_atleast3$color[is.na(test_atleast3$color) & test_atleast3$score>=min_score_atleast3] = "Prediction"
    test_atleast3$color[is.na(test_atleast3$color)] = "Other"
    ## Exactly 3
    test_exactly3 = test %>%
        mutate(color=ifelse(gene_type=="cgc" & (cancer_gene_type=="Recessive" | is.na(cancer_gene_type)), "Recessive",
                            ifelse(gene_type=="cgc" & (cancer_gene_type=="Dominant" & !is.na(cancer_gene_type)), "Dominant", NA)))
    sys_cans = sys_cans %>% select(cancer_type, sample, entrez) %>% mutate(check="Prediction")
    test_exactly3 = test_exactly3 %>% left_join(sys_cans)
    test_exactly3$color[is.na(test_exactly3$color) & test_exactly3$check=="Prediction"] = "Prediction"
    test_exactly3$color[is.na(test_exactly3$color)] = "Other"

    p1 = ggplot(test_atleast3, aes(X1, X2)) + geom_point(aes(shape=type, color=color)) +
        geom_point(data=test_atleast3%>%subset(type=="C"), aes(X1, X2, color=color, text=paste("symbol/sample:", symbol, "/", sample))) +
        scale_colour_manual(values = c("green", "grey", "red", "black")) +
        geom_density2d(data=test_atleast3%>%subset(type=="C"),contour = T) + ## density based on training
        theme_boss() +
        ggtitle(paste0("All genes - at least 3 drivers \n No genes:", test_atleast3%>%subset(score>=min_score_atleast3 | is.na(score))%>%select(entrez)%>%unique%>%nrow,
                       "\n No redundant genes:", test_atleast3%>%subset(score>=min_score_atleast3 | is.na(score))%>%nrow))

    p2 = ggplot(test_exactly3, aes(X1, X2)) + geom_point(aes(shape=type, color=color)) +
        geom_point(data=test_exactly3%>%subset(type=="C"), aes(X1, X2, color=color)) +
        scale_colour_manual(values = c("green", "grey", "red", "black")) +
        geom_density2d(data=test_exactly3%>%subset(type=="C"),contour = T) + ## density based on training
        theme_boss() +
        ggtitle(paste0("All genes - exactly 3 drivers \n No genes:", test_exactly3%>%subset(color=="Prediction" | is.na(score))%>%select(entrez)%>%unique%>%nrow,
                       "\n No redundant genes:",  test_exactly3%>%subset(color=="Prediction" | is.na(score))%>%nrow))

    grid.arrange(p1,p2, nrow=1)


    ## ExpT_NE is numeric so 1 factor is 2 now
    test_no_NE = test%>%subset(!(tissue_exp_level=="Not Expressed" & ExpT_NE==2) | is.na(tissue_exp_level))
    min_score_atleast3 = 10000000
    sys_cans_no_NE = NULL
    for (s in unique(test_no_NE$sample)){
        cat(s, "\n")
        d = test_no_NE %>% subset(sample==s) %>% mutate(score=ifelse(is.na(score), max_score+1, score)) %>% arrange(sample, desc(score)) %>% slice(1:3)
        min_score = min(d$score[!is.na(d$score)])
        if (min_score < min_score_atleast3){
            min_score_atleast3 = min_score
        }
        sys_cans_no_NE = rbind(sys_cans_no_NE, d)
    }
    ## At least 3
    test_no_NE_atleast3 = test_no_NE %>%
        mutate(color=ifelse(gene_type=="cgc" & (cancer_gene_type=="Recessive" | is.na(cancer_gene_type)), "Recessive",
                            ifelse(gene_type=="cgc" & (cancer_gene_type=="Dominant" & !is.na(cancer_gene_type)), "Dominant", NA)))
    test_no_NE_atleast3$color[is.na(test_no_NE_atleast3$color) & test_no_NE_atleast3$score>=min_score_atleast3] = "Prediction"
    test_no_NE_atleast3$color[is.na(test_no_NE_atleast3$color)] = "Other"
    ## Exactly 3
    test_no_NE_exactly3 = test_no_NE %>%
        mutate(color=ifelse(gene_type=="cgc" & (cancer_gene_type=="Recessive" | is.na(cancer_gene_type)), "Recessive",
                            ifelse(gene_type=="cgc" & (cancer_gene_type=="Dominant" & !is.na(cancer_gene_type)), "Dominant", NA)))
    sys_cans_no_NE = sys_cans_no_NE %>% select(cancer_type, sample, entrez) %>% mutate(check="Prediction")
    test_no_NE_exactly3 = test_no_NE_exactly3 %>% left_join(sys_cans_no_NE)
    test_no_NE_exactly3$color[is.na(test_no_NE_exactly3$color) & test_no_NE_exactly3$check=="Prediction"] = "Prediction"
    test_no_NE_exactly3$color[is.na(test_no_NE_exactly3$color)] = "Other"


    p3 = ggplot(test_no_NE_atleast3, aes(X1, X2)) + geom_point(aes(shape=type, color=color)) +
        geom_point(data=test_no_NE_atleast3%>%subset(type=="C"), aes(X1, X2, color=color)) +
        scale_colour_manual(values = c("green", "grey", "red", "black")) +
        geom_density2d(data=test_no_NE_atleast3%>%subset(type=="C"),contour = T) + ## density based on training
        theme_boss() +
        ggtitle(paste0("No not expressed genes - at least 3 drivers \n No genes:", test_no_NE_atleast3%>%subset(score>=min_score_atleast3 | is.na(score))%>%select(entrez)%>%unique%>%nrow,
                       "\n No redundant genes:", test_no_NE_atleast3%>%subset(score>=min_score_atleast3 | is.na(score))%>%nrow))

    p4 = ggplot(test_no_NE_exactly3, aes(X1, X2)) + geom_point(aes(shape=type, color=color)) +
        geom_point(data=test_no_NE_exactly3%>%subset(type=="C"), aes(X1, X2, color=color)) +
        scale_colour_manual(values = c("green", "grey", "red", "black")) +
        geom_density2d(data=test_no_NE_exactly3%>%subset(type=="C"),contour = T) + ## density based on training
        theme_boss() +
        ggtitle(paste0("No not expressed genes - exactly 3 drivers \n No genes:", test_no_NE_exactly3%>%subset(color=="Prediction" | is.na(score))%>%select(entrez)%>%unique%>%nrow,
                       "\n No redundant genes:",  test_no_NE_exactly3%>%subset(color=="Prediction" | is.na(score))%>%nrow))

    grid.arrange(p1,p2,p3,p4, nrow=2)

    ## Check from which samples
    test_atleast3 %>% subset(color!="Other") %>% group_by(sample) %>% summarise(all=n()) %>% left_join(test_no_NE %>% subset(color!="Other") %>% group_by(sample) %>% summarise(no_NE=n()))

    ## Get features and scores (for exporting - no scaling in the features)
    load("/Volumes/mourikisa/novel_driver_prediction/OAC_without_young_plusALLmuts_plusHIC/cgc_518_training/training_set_noScale.Rdata")
    load("/Volumes/mourikisa/novel_driver_prediction/OAC_without_young_plusALLmuts_plusHIC/cgc_518_training/validation_set_noScale.Rdata")
    cohort = rbind(training_ns%>%tibble::rownames_to_column()%>%separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>% data.frame(),
                   validation_ns%>%tibble::rownames_to_column()%>%separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")%>%mutate(type="P")%>%data.frame())
    cohort = cohort %>% mutate(entrez=as.numeric(entrez)) %>% left_join(without_young_plusALL_plusHIC_score[["scores"]])
    cohort$gene_type[is.na(cohort$gene_type)] = "cgc"

    ## Molecular properties
    mol_properties = c("no_ALL_muts", "no_TRUNC_muts", "no_NTDam_muts", "no_GOF_muts", "Copy_number", "CNVLoss", "ExpT_ME", "ExpT_HE", "ExpT_LE", "ExpT_NE", "ExpT_NET")

    ## Annotate and export the table
    export = cohort %>% select(cancer_type, sample, entrez, type, gene_type, score, kernels) %>% left_join(annotation[["geneInfo"]])
    export = export %>% left_join(cohort%>%select(cancer_type, sample, entrez, one_of(mol_properties))) %>% left_join(cohort%>%select(-one_of(mol_properties), -type, -gene_type, -score, -score_scaled, -kernels, -kernels_no))
    ## Get 44
    genes_44 = read.xlsx(file="/Volumes/mourikisa/novel_driver_prediction/OAC/cgc_518_training/shared_genes_annotation_0.95_intersection.xlsx", 1) %>%
        select(entrez) %>% unique %>% .$entrez
    export = export %>% mutate(in_44=ifelse(entrez%in%genes_44, TRUE, FALSE))
    export = export %>% left_join(annotation[["gtex"]]%>%subset(tissue==TISSUE)%>%select(entrez, exp_level)%>%rename(tissue_exp_level=exp_level))

    ## Add in_Cambridge annotation column
    reb_genes = read.xlsx("/Volumes/mourikisa/novel_driver_prediction/OAC/31_Cambridge_genes_annotations.xlsx", 1)
    reb_genes = reb_genes %>% select(entrez) %>% unique %>% .$entrez
    export = export %>% mutate(in_Cambridge_31=ifelse(entrez%in%reb_genes, TRUE, FALSE))

    ## Add false positive annotation
    fp_dir="/Volumes/mourikisa/data/NCG_false_positives.txt"
    false_positive_genes <- read.table(fp_dir, header = T, sep = "\t")
    false_positive_genes <- false_positive_genes %>% select(entrez) %>% .$entrez
    export = export %>% mutate(in_false_positives=ifelse(entrez%in%false_positive_genes, TRUE, FALSE))

    ## Add mutation frequency to the table
    ## Load OAC dataset
    load("/Volumes/mourikisa/data/OAC/Rdata/dataset_mutations_damaging.Rdata")
    oac_MutFreq = df_mut %>% subset(!is.na(entrez_19014)) %>% subset(damaging | truncating | gof) %>%
        select(entrez_19014, sample, ReadCount, VariantAlleleCount, ReadCountControl) %>% rename(entrez=entrez_19014) %>%
        mutate(MutFreq=round(VariantAlleleCount/ReadCount, digits=2)) %>% group_by(entrez, sample) %>% summarise(MutFreq=paste(MutFreq, collapse=","))
    ## Add them to export
    export = export %>% left_join(oac_MutFreq)

    ## Get at least 3,...,8 and add columns to export
    ## Get the minimum score for at least 3 predictions
    max_score = max(test$score[!is.na(test$score)])
    min_scores = NULL
    for (i in 3:8){
        cat(i, "\n")
        min_score_atleast = 10000000
        for (s in unique(test$sample)){
            cat(s, "\n")
            d = test %>% subset(sample==s) %>% mutate(score=ifelse(is.na(score), max_score+1, score)) %>% arrange(sample, desc(score)) %>% slice(1:i)
            min_score = min(d$score[!is.na(d$score)])
            if (min_score < min_score_atleast){
                min_score_atleast = min_score
            }
        }
        min_scores = rbind(min_scores, c(i, min_score_atleast))
    }
    min_scores = data.frame(min_scores)
    colnames(min_scores) = c("at_least", "min_score")
    ## Add the columns for at least predictions
    export = export %>% mutate(at_least_3=ifelse(score >= min_scores$min_score[min_scores$at_least==3] & !is.na(score), TRUE, FALSE))
    export = export %>% mutate(at_least_4=ifelse(score >= min_scores$min_score[min_scores$at_least==4] & !is.na(score), TRUE, FALSE))
    export = export %>% mutate(at_least_5=ifelse(score >= min_scores$min_score[min_scores$at_least==5] & !is.na(score), TRUE, FALSE))
    export = export %>% mutate(at_least_6=ifelse(score >= min_scores$min_score[min_scores$at_least==6] & !is.na(score), TRUE, FALSE))
    export = export %>% mutate(at_least_7=ifelse(score >= min_scores$min_score[min_scores$at_least==7] & !is.na(score), TRUE, FALSE))
    export = export %>% mutate(at_least_8=ifelse(score >= min_scores$min_score[min_scores$at_least==8] & !is.na(score), TRUE, FALSE))

    ## Add the coordinates as well
    export = export %>% left_join(test%>%select(entrez, sample, X1, X2, cancer_gene_type))

    #export = export %>% left_join(test%>%select(cancer_type, sample, entrez, X1, X2))
    write.xlsx(export%>%data.frame(), file="/Volumes/FCShared/Thanos/CRUK/1489_79_OAC_scores_plus_annotation_0208.xlsx", row.names = F)

    ## export the minimum values
    write.xlsx(min_scores%>%data.frame(), file="/Volumes/FCShared/Thanos/CRUK/1489_79_OAC_minimum_scores.xlsx")


}

## Calculate the score of each gene
getScore = function(cancer_type=NULL,
                    training="cgc",
                    path=NULL,
                    linear_preds=NULL, linearCVS=NULL, linearVar=NULL,
                    radial_preds=NULL, radialCVS=NULL, radialVar=NULL,
                    polynomial_preds=NULL, polynomialCVS=NULL, polynomialVar=NULL,
                    sigmoid_preds=NULL, sigmoidCVS=NULL, sigmoidVar=NULL){

    if (is.null(cancer_type)| is.null(path)){
        cat("Parameters missing...", "\n")
        return(NULL)
    }

    cat(paste0("Cancer_type: ", cancer_type), "\n")
    cat(paste0("Training: ", training), "\n")
    ## Get the number of predictions before intersection
    all_preds = NULL
    all_scores = NULL

    if(!is.null(linear_preds) & !is.null(linearCVS)){
        load(paste0(path,linear_preds))
        linear = preds
        rm(preds)
        ## Push it to the all_preds
        all_preds = rbind(all_preds, linear %>% subset(label) %>% mutate(kernel="linear", training=training) %>% data.frame())

        ## Kernel-specific score
        linear = linear %>% arrange(sample, dv) %>% group_by(sample) %>% mutate(rank=row_number(-dv)) %>% ungroup
        linear = linear %>% group_by(sample) %>% mutate(N=max(rank), Nlog10=log10(N)) %>% ungroup
        linear = linear %>% mutate(cvs=linearCVS, var=linearVar) %>% mutate(ratio=cvs/var)
        linear = linear %>% mutate(score=-log10(rank/N)*cvs)
        linear = linear %>% mutate(kernel="linear") %>% data.frame()
        ## Old score implementation
        #linear = linear %>% mutate(score=(-10*log10(dv_rank/nrow(linear)))*(linearCVS/linearVar))
        #linear = linear %>% mutate(score=(-10*log10(dv_rank/nrow(linear)))*linearCVS)
        #min_score = min(linear$score)
        #max_score = max(linear$score)
        #linear = linear %>% mutate(score_scaled = (score-min_score)/(max_score-min_score), kernel="linear", training=training)
        all_scores = rbind(all_scores, linear)
    }

    if(!is.null(radial_preds) & !is.null(radialCVS)){
        load(paste0(path,radial_preds))
        radial = preds
        rm(preds)
        ## Push it to the all_preds
        all_preds = rbind(all_preds, radial %>% subset(label) %>% mutate(kernel="radial", training=training) %>% data.frame())

        ## Kernel-specific score
        radial = radial %>% arrange(sample, dv) %>% group_by(sample) %>% mutate(rank=row_number(-dv)) %>% ungroup
        radial = radial %>% group_by(sample) %>% mutate(N=max(rank), Nlog10=log10(N)) %>% ungroup
        radial = radial %>% mutate(cvs=radialCVS, var=radialVar) %>% mutate(ratio=cvs/var)
        radial = radial %>% mutate(score=-log10(rank/N)*cvs)
        radial = radial %>% mutate(kernel="radial") %>% data.frame()
        ## Old score implementation
        #radial = radial %>% mutate(score=(-10*log10(dv_rank/nrow(radial)))*(radialCVS/radialVar))
        #radial = radial %>% mutate(score=(-10*log10(dv_rank/nrow(radial)))*radialCVS)
        #min_score = min(radial$score)
        #max_score = max(radial$score)
        #radial = radial %>% mutate(score_scaled = (score-min_score)/(max_score-min_score), kernel="radial", training=training)
        all_scores = rbind(all_scores, radial)
    }

    if(!is.null(polynomial_preds) & !is.null(polynomialCVS)){
        load(paste0(path,polynomial_preds))
        polynomial = preds
        rm(preds)
        ## Push it to the all_preds
        all_preds = rbind(all_preds, polynomial %>% subset(label) %>% mutate(kernel="polynomial", training=training) %>% data.frame())

        ## Kernel-specific score
        polynomial = polynomial %>% arrange(sample, dv) %>% group_by(sample) %>% mutate(rank=row_number(-dv)) %>% ungroup
        polynomial = polynomial %>% group_by(sample) %>% mutate(N=max(rank), Nlog10=log10(N)) %>% ungroup
        polynomial = polynomial %>% mutate(cvs=polynomialCVS, var=polynomialVar) %>% mutate(ratio=cvs/var)
        polynomial = polynomial %>% mutate(score=-log10(rank/N)*cvs)
        polynomial = polynomial %>% mutate(kernel="polynomial") %>% data.frame()
        #polynomial = polynomial %>% mutate(score=(-10*log10(dv_rank/nrow(polynomial)))*(polynomialCVS/polynomialVar))
        #polynomial = polynomial %>% mutate(score=(-10*log10(dv_rank/nrow(polynomial)))*polynomialCVS)
        #min_score = min(polynomial$score)
        #max_score = max(polynomial$score)
        #polynomial = polynomial %>% mutate(score_scaled = (score-min_score)/(max_score-min_score), kernel="polynomial", training=training)
        all_scores = rbind(all_scores, polynomial)
    }
    if(!is.null(sigmoid_preds) & !is.null(sigmoidCVS)){
        load(paste0(path,sigmoid_preds))
        sigmoid = preds
        rm(preds)
        ## Push it to the all_preds
        all_preds = rbind(all_preds, sigmoid %>% subset(label) %>% mutate(kernel="sigmoid", training=training) %>% data.frame())

        ## Kernel-specific score
        sigmoid = sigmoid %>% arrange(sample, dv) %>% group_by(sample) %>% mutate(rank=row_number(-dv)) %>% ungroup
        sigmoid = sigmoid %>% group_by(sample) %>% mutate(N=max(rank), Nlog10=log10(N)) %>% ungroup
        sigmoid = sigmoid %>% mutate(cvs=sigmoidCVS, var=sigmoidVar) %>% mutate(ratio=cvs/var)
        sigmoid = sigmoid %>% mutate(score=-log10(rank/N)*cvs)
        sigmoid = sigmoid %>% mutate(kernel="sigmoid") %>% data.frame()
        #sigmoid = sigmoid %>% mutate(score=(-10*log10(dv_rank/nrow(sigmoid)))*(sigmoidCVS/sigmoidVar))
        #sigmoid = sigmoid %>% mutate(score=(-10*log10(dv_rank/nrow(sigmoid)))*sigmoidCVS)
        #min_score = min(sigmoid$score)
        #max_score = max(sigmoid$score)
        #sigmoid = sigmoid %>% mutate(score_scaled = (score-min_score)/(max_score-min_score), kernel="sigmoid", training=training)
        all_scores = rbind(all_scores, sigmoid)
    }

    ## Final score for each gene across kernels
    scores = all_scores %>% group_by(cancer_type, sample, entrez, gene_type, N, Nlog10) %>%
        summarise(sum_of_scores=sum(score), kernels_used=length(unique(kernel))) %>%
        mutate(score=(1/(kernels_used*Nlog10))*sum_of_scores) %>% ungroup %>% data.frame()



    ## Add from how many kernels they have predicted as TRUE
    preds2kernels = all_preds %>% subset(label) %>% group_by(cancer_type, sample, entrez) %>% summarise(kernels_predicted=paste(unique(kernel), collapse=","), kernels_predicted_no=length(unique(kernel))) %>% ungroup
    scores = scores %>% left_join(preds2kernels)
    scores$kernels_predicted_no[is.na(scores$kernels_predicted_no)] = 0

    return(list(all_scores=all_scores, scores=scores))

}

## A function to get the value of objective function in libsvm
## Usage from http://stats.stackexchange.com/questions/25387/problem-with-e1071-libsvm
# a <- svm(x, y, scale = FALSE, type = 'C-classification', kernel = 'linear', cost = 50000)
# w <- t(a$coefs) %*% a$SV;
# b <- -a$rho;
# obj_func_str1 <- get_obj_func_info(w, b, 50000, x, y)
# obj_func_str2 <- get_obj_func_info(w, b - 5, 50000, x, y)
get_obj_func_info <- function(w, b, c_par, x, y) {
    xi <- rep(0, nrow(x))

    for (i in 1:nrow(x)) {
        xi[i] <- 1 - as.numeric(as.character(y[i]))*(sum(w*x[i,]) + b)
        if (xi[i] < 0) xi[i] <- 0
    }

    return(list(obj_func_value = 0.5*sqrt(sum(w * w)) + c_par*sum(xi),
                sum_xi = sum(xi), xi = xi))
}

mapply_pb = function(FUN, X, Y,  ...){
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
  res <- mapply(wrapper, X, Y, ...)
  close(pb)
  res
}


make_pathway_list = function(pathMsigDbFile) {
  inputFile <- pathMsigDbFile
  con  <- file(inputFile, open = "r")
  c = 1
  pathway.list <- vector(mode="list",length=0)
  print("Loading gene sets....")
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    myVector <- do.call("rbind",strsplit(oneLine, "\t"))
    t = vector(mode="list",length=1)
    t[[1]] = myVector[3:length(myVector)]
    names(t) = myVector[1]
    pathway.list = c(pathway.list,t)
    c = c+1
  }

  close(con)
  return(pathway.list)
}

base_breaks_x = function(x, br, la){
  d <- data.frame(y=-Inf, yend=-Inf, x=0, xend=100)
  list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE), scale_x_reverse(breaks=br, labels=la))
}

base_breaks_y = function(m, br, la ){
  d <- data.frame(x=Inf, xend=Inf, y=0, yend=m)
  list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE), scale_y_continuous(breaks=br, labels=la))
}

map_to_geneset <- function ( x , tissue,
                             geneset="geneSets/c2.cp.kegg.v5.1.entrez.gmt.txt",
                             fout="KEGG_intersection_best_models.xlsx",
                             fpdf="KEGG_intersection_best_models.pdf",
                             install_NCG=T,
                             NCG_version="NCG5",
                             NCG_path=NULL,
                             uninformative_paths=c('CANCER','LEUKEMIA')) {


  kegg = make_pathway_list(geneset)
  len = sapply(kegg, length)

  sel.genes = unique(x[,c('entrez','gene_type')])
  sel.genes$path = sapply(as.list(as.character(unique(x$entrez))), function(x,y){
    sel=sapply(y, function(z,w) w%in%z, w = x);
    return( paste(names(y)[sel],collapse='.')) },
    y = kegg)

  dev_mode()
  if(install_NCG & !is.null(NCG_path)) install(NCG_path)
  library(ncglib)
  geneInfo=as.data.frame(get_geneinfo(version = NCG_version), stringsAsFactors=F)
  cancer_genes = as.data.frame(get_cancer_genes(version = NCG_version))
  dev_mode()

  sel.genes$symbol = geneInfo[match(sel.genes$entrez, geneInfo$entrez), "symbol"]

  ## refine a data frame to return for downstream analysis - part 1
  exp1 = sel.genes %>% select(entrez, symbol, gene_type, path) %>%
      mutate(type="prediction", primary_site=NA)  %>% data.frame()

  y = sapply(as.list(as.character(unique(cancer_genes$entrez))), function(x,y){
      sel=sapply(y, function(z,w) w%in%z, w = x);
      return( paste(names(y)[sel],collapse='.')) },
      y = kegg)

  # cancer genes
  can.genes = data.frame(entrez=unique(cancer_genes$entrez), path = y)
  can.genes$cancer_type = geneInfo[match(can.genes$entrez, geneInfo$entrez),'cancer_type']

  allGenes = geneInfo %>% subset(duplicability!=-1 | cancer_type=="can" | cancer_type=="cgc") %>% select(entrez, symbol)
  y = sapply(as.list(as.character(unique(allGenes$entrez))), function(x,y){
    sel=sapply(y, function(z,w) w%in%z, w = x);
    return( paste(names(y)[sel],collapse='.')) },
    y = kegg)

  ## Get all the 19014 genes
  allGenes = data.frame(entrez=unique(allGenes$entrez), path = y)
  allGenes$cancer_type = geneInfo[match(allGenes$entrez, geneInfo$entrez),'cancer_type']


  ## refine a data frame to return for downstream analysis - part 2

#   exp2 = can.genes %>% rename(gene_type=cancer_type) %>%
#       left_join(cancer_genes%>%select(entrez, primary_site), by=c("entrez")) %>%
#       group_by(entrez, path, gene_type) %>% summarise(primary_site=paste0(unique(primary_site), collapse=",")) %>%
#       left_join(geneInfo%>%select(entrez,symbol), by=c("entrez")) %>%
#       mutate(type="NCG5") %>%
#       select(entrez, symbol, gene_type, path, type, primary_site)

  exp2 = allGenes %>% rename(gene_type=cancer_type) %>%
      left_join(cancer_genes%>%select(entrez, primary_site), by=c("entrez")) %>%
      group_by(entrez, path, gene_type) %>% summarise(primary_site=paste0(unique(primary_site), collapse=",")) %>%
      left_join(geneInfo%>%select(entrez,symbol), by=c("entrez")) %>%
      mutate(type=ifelse(gene_type=="rst", "rst", "NCG5")) %>%
      select(entrez, symbol, gene_type, path, type, primary_site) %>% data.frame()

  exp = rbind(exp1,exp2)

  unique_path = unique(unlist(strsplit(can.genes$path, "\\.")))
  m = matrix(0, nr = nrow(can.genes), nc = length(unique_path), dimnames = list(can.genes$entrez, unique_path))
  for(i in unique_path) m[grep(i, can.genes$path), i ] =1
  df = data.frame(value=sort(apply(m, 2, sum)))
  df$path=rownames(df)
  df = unrowname(df)
  df = df[,c("path","value")]
  df$path_size = len[df$path]
  df$perc = with(df, (value/path_size)*100 )
  df$category = "NCG_cancer_genes_1571"
  res = df

  unique_path = unique(unlist(strsplit(subset(can.genes, cancer_type=='cgc')$path, "\\.")))
  m = matrix(0, nr = nrow(subset(can.genes, cancer_type=='cgc')), nc = length(unique_path), dimnames = list(subset(can.genes, cancer_type=='cgc')$entrez, unique_path))
  for(i in unique_path) m[grep(i, subset(can.genes, cancer_type=='cgc')$path), i ] =1
  df = data.frame(value=sort(apply(m, 2, sum)))
  df$path=rownames(df)
  df = unrowname(df)
  df = df[,c("path","value")]
  df$path_size = len[df$path]
  df$perc = with(df, (value/path_size)*100 )
  df$category = "NCG_cancer_genes_518"
  res = rbind(res,df)

  can.genes = data.frame(entrez=unique(subset(cancer_genes, primary_site==tissue)$entrez),
                         path = sapply(as.list(as.character(unique(subset(cancer_genes, primary_site==tissue)$entrez))), function(x,y){
    sel=sapply(y, function(z,w) w%in%z, w = x);
    return( paste(names(y)[sel],collapse='.')) },
    y = kegg))

  can.genes$cancer_type = geneInfo[match(can.genes$entrez, geneInfo$entrez),'cancer_type']

  unique_path = unique(unlist(strsplit(can.genes$path, "\\.")))
  m = matrix(0, nr = nrow(can.genes), nc = length(unique_path), dimnames = list(can.genes$entrez, unique_path))
  for(i in unique_path) m[grep(i, can.genes$path), i ] =1
  df = data.frame(value=sort(apply(m, 2, sum)))
  df$path=rownames(df)
  df = unrowname(df)
  df = df[,c("path","value")]
  df$path_size = len[df$path]
  df$perc = with(df, (value/path_size)*100 )
  df$category = paste0("NCG_cancer_genes_",nrow(can.genes),"_",tissue)
  res = rbind(res, df)


#   df$perc = 100*(df$value/nrow(can.genes))
#   df$label = paste0(df$path," ( ", df$value," / ", len[df$path], " )" )
#   df$path = factor(df$path, levels = names(sort(apply(m, 2, sum))))
#   df$label = factor(df$label, levels = unique(df$label[match(levels(df$path), df$path)]))

#   pdf(file=pdfname,h=pdf_height,w=pdf_width)
#
#   print(ggplot(df, aes(x=label,y=perc))+geom_bar(stat='identity')+coord_flip()+theme_boss()+ylab("")+xlab("")+
#           ggtitle(paste0(nrow(can.genes)," cancer genes of ", tissue)))


  # all novel prediction

  unique_path = unique(unlist(strsplit(sel.genes$path, "\\.")))
  m = matrix(0, nr = nrow(sel.genes), nc = length(unique_path), dimnames = list(sel.genes$entrez, unique_path))
  for(i in unique_path) m[grep(i, sel.genes$path), i ] =1

  df = data.frame(value=sort(apply(m, 2, sum)))
  df$path=rownames(df)
  df = unrowname(df)
  df = df[,c("path","value")]
  df$path_size = len[df$path]
  df$perc = with(df, (value/path_size)*100 )
  df$category = paste0("Novel_prediction_ALL_",nrow(sel.genes))
  res = rbind(res, df)

  # print(ggplot(df, aes(x=label,y=perc))+geom_bar(stat='identity')+coord_flip()+theme_boss()+ylab("")+xlab("")+ggtitle(paste0(nrow(sel.genes)," novel genes (all)")))

  # only novel rst

  unique_path = unique(unlist(strsplit(subset(sel.genes, gene_type=="rst")$path, "\\.")))
  m = matrix(0, nr = nrow(subset(sel.genes, gene_type=="rst")), nc = length(unique_path), dimnames = list(subset(sel.genes, gene_type=="rst")$entrez, unique_path))
  for(i in unique_path) m[grep(i, subset(sel.genes, gene_type=="rst")$path), i ] =1

  df = data.frame(value=sort(apply(m, 2, sum)))
  df$path=rownames(df)
  df = unrowname(df)
  df = df[,c("path","value")]
  df$path_size = len[df$path]
  df$perc = with(df, (value/path_size)*100 )
  df$category = paste0("Novel_prediction_RST_",nrow(subset(sel.genes, gene_type=="rst")))
  res = rbind(res, df)

  # print(ggplot(df, aes(x=label,y=perc))+geom_bar(stat='identity')+coord_flip()+theme_boss()+ylab("")+xlab("")+ggtitle(paste0(nrow(subset(sel.genes, gene_type=="rst"))," novel genes (rst)")))


  # dev.off()
  write.xlsx(subset(res, category=="NCG_cancer_genes_1571"), "NCG_cancer_genes_1571", file=fout, row.names = F, showNA = F)
  for (i in unique(res$category)[!unique(res$category)%in%"NCG_cancer_genes_1571"]){
    write.xlsx(subset(res, category==i), i, file=fout, row.names = F, showNA = F, append=T)
  }

  ## Rank-based measure of association (Spearman's rho)
  l = list()
  u = unique(res$category)
  for(i in 1:length(u)){
      l[[i]] = subset(res, category==u[i])$path
      l[[i]] = l[[i]][ with(subset(res, category==u[i]), order(value, perc, decreasing = T) ) ]
  }

  get_hard_core = function(x,y,xtitle,ytitle){
         # check if x is longer than y

      a = as.numeric(factor(x, levels=x))
      b = as.numeric(factor(y, levels=y))

      names(a) = x
      names(b) = y

      dfxy = data.frame(a=a, b=b[names(a)])
      pattern = paste(unique(uninformative_paths), collapse="|")
      dfxy = dfxy %>% tibble::rownames_to_column() %>% mutate(informative=ifelse(grepl(pattern, rowname), FALSE, TRUE)) %>% data.frame()
      dfxy$informative = factor(dfxy$informative, levels = c(TRUE, FALSE))
      #scatter.smooth(a, b[names(a)]) # geom_smooth
      options(scipen=2)
      pv = cor.test( dfxy$a, dfxy$b, m="spearman")$p.value
      rho = cor.test( dfxy$a, dfxy$b, m="spearman")$estimate

      ## Keep colours steady irrespectively from which level of the factor you find first
      myColors = c("red", "black")
      names(myColors) = levels(dfxy$informative)
      p = ggplot(dfxy, aes(x=a, y=b)) + geom_point(aes(colour = informative), size=rel(1.5)) + geom_smooth() +
          #geom_text(aes(label=ifelse(a>10 & b<10 & informative==TRUE,as.character(rowname),'')),hjust=0,vjust=0,size=rel(3)) +
          geom_text(aes(label=ifelse(informative==TRUE,as.character(rowname),'')),hjust=0,vjust=0,size=rel(3)) +
          #geom_text(aes(label=as.character(rowname)),hjust=0,vjust=0,size=rel(2.5)) +
          scale_colour_manual(values=myColors, name="Informative\npathways") +
          theme_boss() +
          theme(
              axis.line = element_line(colour = "black")
          ) +
          xlab(paste0(xtitle, " KEGG pathway ranks")) +
          ylab(paste0(ytitle, " KEGG pathway ranks")) +
          ggtitle(paste0("p-value=", pv, "; rho=", rho))


#       +
#           base_breaks_y(ymax, c(0,ymed,ymax), as.character(c(0,lmed,lmax)))+
#           base_breaks_x(gg$x, c(0,35,80,100), c("0","35","80","100"))
      return(p)
  }

  names(l) = u
  p1 = get_hard_core(l[[1]],l[[2]], names(l)[1], names(l)[2])
  p2 = get_hard_core(l[[1]],l[[3]], names(l)[1], names(l)[3])
  p3 = get_hard_core(l[[3]],l[[4]], names(l)[3], names(l)[4])
  p4 = get_hard_core(l[[3]],l[[5]], names(l)[3], names(l)[5])
  pdf(file=paste(fpdf, sep=""), h= 20, w=24)
  grid.arrange(p1, p2, p3, p4, ncol=2)
  dev.off()

  ## Attention if a gene is a cgc but not in the tissue examined (e.g blood)
  ## it exists twice in the table returned, once for the predictions and
  ## once for the NCG5 genes (in the second case we give primary site...etc)
  return(exp)
}


# rad = get_stat_novel_predictions("/Volumes/ceredam/novel_driver_prediction/ESCA/gridPrediction/",
#                                  "predictions_radial",
#                                  "Results/SVM_radial.novel_by_models.pdf")

get_stat_novel_predictions = function(path, kern="linear", gene_type, res_dir, plot="barplot"){

    specify_decimal <- function(x, k) format(round(x, k), nsmall=k)
    if (!(gene_type%in%c("all", "rst", "kcg"))){
        gene_type = "all"
        cat("----> get_stat_novel_predictions: unknown gene_type, it was set to \"all\"", "n")
    }

    if (kern=="linear"){pattern="predictions_linear"}
    if (kern=="radial"){pattern="predictions_radial"}
    if (kern=="polynomial"){pattern="predictions_polynomial"}
    if (kern=="sigmoid"){pattern="predictions_sigmoid"}

    if (gene_type=="rst"){genes="rst"}
    if (gene_type=="all"){genes="all"}
    if (gene_type=="kcg"){genes="NCG5"}

    ff = list.files(path, pattern, full.names = T)
    mm= list.files(path , pattern, full.names = F)
    ff = ff[-grep("xlsx", ff)]
    mm = mm[-grep("xlsx", mm)]
    mm = gsub(".Rdata","", mm)
    mm = gsub("predictions_","", mm)


    novel = list()
    for(i in 1:length(ff)){
      print(i)
      load(ff[i])
      if (gene_type=="rst"){
        rst = subset(as.data.frame(preds), gene_type=="rst" & label & prob_platt>.95)
        novel[[i]] = unique(rst$entrez)
      }else if (gene_type=="all"){
        all = subset(as.data.frame(preds), label & prob_platt>.95)
        novel[[i]] = unique(all$entrez)
      }else if (gene_type=="kcg"){
        kcg = subset(as.data.frame(preds), (gene_type=="cgc" | gene_type=="can") & label & prob_platt>.95)
        novel[[i]] = unique(kcg$entrez)
      }
    }

    names(novel)=mm


    if (plot=="barplot"){

        ## No predictions - clean novel list from 0s
        ix = which(lapply(novel, function(x) length(x))==0)
        cat(paste0(length(ix), " models had no predictions and excluded..."), "\n")
        novel[ix] = NULL

        print(length(unique(unlist(novel))))
        ## Filter mm also for models with 0 predictions
        mm = mm[!(mm%in%names(ix))]

        novel = mapply(function(x,y){ cbind('entrez'=x, 'model'=y) }, novel, as.list(mm), SIMPLIFY = T)
        df = do.call(rbind,novel)
        df = as.data.frame(df)
        df$previous = NA
        for(i in 2:length(mm) ){
            df$previous[ which(df$model == mm[i] ) ] = df$entrez[ which(df$model == mm[i] ) ]%in%df$entrez[ which(df$model == mm[i-1] ) ]
        }
        sdf = ddply(df, .(model), summarise, tot = length(entrez), previous = sum(previous))
        sdf$previous[1] = 0
        sdf$private = sdf$tot - sdf$previous
        sdf$p.private = (sdf$private/sdf$tot )*100
        sdf$p.previous = (sdf$previous/sdf$tot )*100
        sdf$label=paste0(sdf$tot,"\n(",sdf$private,")")
        mdf = melt(sdf, id.vars = "model", measure.vars = c("private","previous"))
        if (kern=="linear"){
            pdf(file=paste(res_dir,"/",gene_type, "_shared_genes_between_models_", kern, ".pdf", sep=""), h=6, w=10)
        }
        if (kern=="radial"){
            pdf(file=paste(res_dir,"/",gene_type, "_shared_genes_between_models_", kern, ".pdf", sep=""), h=6, w=40)
        }
        if (kern=="polynomial"){
            pdf(file=paste(res_dir,"/",gene_type, "_shared_genes_between_models_", kern, ".pdf", sep=""), h=20, w=150)
        }
        if (kern=="sigmoid"){
            pdf(file=paste(res_dir,"/",gene_type, "_shared_genes_between_models_", kern, ".pdf", sep=""), h=6, w=40)
        }
        p = ggplot(mdf, aes(x=model,y=value))+
            geom_bar(stat="identity", col="black", aes( fill=variable))+theme_boss()+scale_fill_manual(values=c("green","grey"))+theme(axis.text.x = element_text(angle=90, size=rel(1.5)))+
            geom_text(data=sdf, aes(x=model,y=tot, label=label),position = position_dodge(width=1), vjust=-0.15, size=rel(4))+ylab("Unique genes (n)")+xlab("")+
            ggtitle( paste0("Novel driver genes (",genes,",Platt) per model\nUnique genes across models = ", length(unique(df$entrez))) )+
            scale_y_continuous(limits=c(0, max(sdf$tot+500)), breaks=c(0,round(max(sdf$tot),0)/2,max(sdf$tot)), labels=as.character(c(0,round(max(sdf$tot)/2,0),max(sdf$tot))))
        print(p)
        dev.off()
        write.table(df, file=paste(res_dir,"/", gene_type, "_shared_genes_between_models_", kern, ".tsv",sep=""), quote=F, row.names = F, sep = "\t")
        write.table(sdf, file=paste(res_dir,"/", gene_type, "_shared_genes_between_models_per_gene_", kern, ".tsv",sep=""), quote=F, row.names = F, sep = "\t")
    }else if (plot=="heatmap"){

#        ## Get all the combinations of the models
#         mcombs = unique(data.frame(do.call(rbind, combn(mm, 2, simplify = FALSE))))
#
#         m1 = as.list(mcombs[,1])
#         m2 = as.list(mcombs[,2])
#
#         x = paste0(gsub("radial.","",mcombs[,1]),"_vs_",gsub("radial.","", mcombs[,2]) )
#
#         names(m1) = x
#         names(m2) = x
#
#         n1 = lapply(m1, function(x,y) y[[x]], y=novel)
#         n2 = lapply(m2, function(x,y) y[[x]], y=novel)
#
#         gimme_money_bad_guy = function(x,y){
#           length(setdiff(y,x)) / length(y)
#         }
#
#         res = mapply_pb( gimme_money_bad_guy, n1, n2, SIMPLIFY = T, USE.NAMES=F)
#         res = unlist(res)
#         y  = as.data.frame(do.call("rbind", strsplit(x, "\\_vs_"))); colnames(y)=c("previous","selected")
#         y$value = res
#         y$value[ which(!is.finite(y$value))] = NA
#         myPalette <- colorRampPalette(c("white",(brewer.pal(7, "YlOrRd"))), space="Lab")
#         r = c(0,max(y$value, na.rm=T))
#         ggplot(df, aes(x=previous, y=selected, fill=value)) +          geom_tile() + coord_equal() +
#           # geom_text(aes(label = as.character(100-as.numeric(df$m1INTERm2Perc2)))) +
#           scale_fill_gradientn(colours = myPalette(100), name='Percentage', limits = r, breaks =r, labels=as.character(r))
#                                ,na.value ="grey80") +
#           #theme_boss() +
#           theme(
#             legend.position = "bottom",
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),
#             panel.background = element_blank(),
#             axis.title.x = element_text(colour=NA),
#             axis.title.y = element_blank(),
#             axis.line = element_blank(),
#             axis.text.x=element_text(face="bold", angle=90, size=10),
#             axis.text.y=element_text(face="bold", size=10)
#           )
#
#
        # novel = NULL
#         for(i in 1:nrow(mcombs)){
#             cat(c(i, " ", unname(unlist(mcombs[i,]))), "\n")
#             ## Load the predictions
#             d1 = ff[grep(paste0(mcombs[i,1], ".Rdata"), ff)]
#             d2 = ff[grep(paste0(mcombs[i,2], ".Rdata"), ff)]
#             load(d1)
#             d1 = preds
#             load(d2)
#             d2 = preds
#
#             ## Rst predictions
#             rst1 = d1 %>% subset(gene_type=="rst" & label & prob_platt>.95) %>% select(entrez) %>% unique %>% .$entrez
#             rst2 = d2 %>% subset(gene_type=="rst" & label & prob_platt>.95) %>% select(entrez) %>% unique %>% .$entrez
#             m1INTERm2 = length(intersect(rst1, rst2))
#             m1SETDIFFm2 = length(setdiff(rst1, rst2))
#             m2SETDIFFm1 = length(setdiff(rst2, rst1))
#             novel = rbind(novel, data.frame(m1=mcombs[i,1], m2=mcombs[i,2],
#                                             m1Genes=length(rst1), m2Genes=length(rst2),
#                                             gene_type="rst", m1INTERm2=m1INTERm2,
#                                             m1SETDIFFm2=m1SETDIFFm2,
#                                             m2SETDIFFm1=m2SETDIFFm1))
#
#
#             ## All (cancer genes and rst predictions)
#             all1 = d1 %>% subset(label & prob_platt>.95) %>% select(entrez) %>% unique %>% .$entrez
#             all2 = d2 %>% subset(label & prob_platt>.95) %>% select(entrez) %>% unique %>% .$entrez
#             m1INTERm2 = length(intersect(all1, all2))
#             m1SETDIFFm2 = length(setdiff(all1, all2))
#             m2SETDIFFm1 = length(setdiff(all2, all1))
#             novel = rbind(novel, data.frame(m1=mcombs[i,1], m2=mcombs[i,2],
#                                             m1Genes=length(all1), m2Genes=length(all2),
#                                             gene_type="all", m1INTERm2=m1INTERm2,
#                                             m1SETDIFFm2=m1SETDIFFm2,
#                                             m2SETDIFFm1=m2SETDIFFm1))
#
#             ## Predictions in NCG5
#             kcg1 = d1 %>% subset((gene_type=="cgc" | gene_type=="can") & label & prob_platt>.95) %>% select(entrez) %>% unique %>% .$entrez
#             kcg2 = d2 %>% subset((gene_type=="cgc" | gene_type=="can") & label & prob_platt>.95) %>% select(entrez) %>% unique %>% .$entrez
#             m1INTERm2 = length(intersect(kcg1, kcg2))
#             m1SETDIFFm2 = length(setdiff(kcg1, kcg2))
#             m2SETDIFFm1 = length(setdiff(kcg2, kcg1))
#             novel = rbind(novel, data.frame(m1=mcombs[i,1], m2=mcombs[i,2],
#                                             m1Genes=length(kcg1), m2Genes=length(kcg2),
#                                             gene_type="NCG5", m1INTERm2=m1INTERm2,
#                                             m1SETDIFFm2=m1SETDIFFm2,
#                                             m2SETDIFFm1=m2SETDIFFm1))
#         }

        ## Get parameters in the df
        #novel %>% separate(m1, into=c())

#         novel %>% mutate(m1INTERm2Perc2) %>%
#             mutate(printm1=paste0(m1," (", m1Genes, ")"), printm2=paste0(m2," (", m2Genes, ")"),
#                    m1INTERm2Perc2=as.numeric(specify_decimal((m1INTERm2/m2Genes)*100, 2)),
#                    m1INTERm2Perc1=as.numeric(specify_decimal((m1INTERm2/m1Genes)*100, 2)))
#
#
#         for (t in unique(novel$gene_type)){
#             if (kern=="linear"){
#                 pdf(file=paste(res_dir,"/",t, "_shared_genes_between_models_", kern, ".pdf", sep=""), h=6, w=10)
#             }
#             if (kern=="radial"){
#                 pdf(file=paste(res_dir,"/",t, "_shared_genes_between_models_", kern, ".pdf", sep=""), h=6, w=40)
#             }
#
#             df = novel %>% subset(gene_type==t)
#             myPalette <- colorRampPalette(c("white",(brewer.pal(7, "YlOrRd"))), space="Lab")
#             r = c(0,100)
#             ggplot(df, aes(x=as.factor(printm1), y=as.factor(printm2))) +
#                 geom_tile(aes(fill = 100-as.numeric(m1INTERm2Perc2))) + coord_equal() +
#                 geom_text(aes(label = as.character(100-as.numeric(df$m1INTERm2Perc2)))) +
#                 scale_fill_gradientn(colours = myPalette(100), name='Percentage', limits = r, breaks =r, labels=as.character(r),na.value ="white") +
#                 #theme_boss() +
#                 theme(
#                     legend.position = "bottom",
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(),
#                     panel.background = element_blank(),
#                     axis.title.x = element_text(colour=NA),
#                     axis.title.y = element_blank(),
#                     axis.line = element_blank(),
#                     axis.text.x=element_text(face="bold", angle=90, size=10),
#                     axis.text.y=element_text(face="bold", size=10)
#                 )
#
#         }

    }





}

stats_genes = function(cancer_type){

  load("/mnt/lustre/users/k1469280/mourikisa/data/geneInfoNCG5.Rdata")
  load("/mnt/lustre/users/k1469280/mourikisa/data/dupByThresholds.Rdata")

  geneInfo = as.data.frame(geneInfo, stringsAsFactors=F)
  dupByThresholds = as.data.frame(dupByThresholds, stringsAsFactors=F)

  d = getTable(tb=cancer_type)

  lof = rbind(lof_trunc = subset(d,no_TRUNC_muts>0, select =c('Cancer_type', 'Sample','Entrez')) %>%mutate(type="trunc"),
              lof_nondam = subset(d, no_NTDam_muts>0, select =c('Cancer_type', 'Sample','Entrez')) %>%mutate(type="ntdam"),
              lof_homdel = subset(d,(CNVLoss==1 & Copy_number==0), select =c('Cancer_type', 'Sample','Entrez')) %>% mutate(type="homdel"),
              lof_hetdelTruncORNTdam = subset(d, (CNVLoss==1 & Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)),select =c('Cancer_type', 'Sample','Entrez')) %>%
              mutate(type="hetdelTruncORNTdam"))

  lof.sample = ddply(lof, .(Cancer_type,Sample), summarise, n_genes=length(unique(Entrez)))

  #lof.genes  = ddply(lof, .(Cancer_type,Entrez), summarise, n_samples=length(unique(Sample)))
  #lof.genes = cbind(lof.genes, geneInfo[ match(lof.genes$Entrez, geneInfo$entrez), c('symbol', 'cancer_type', 'complex')]) # , 'duplicability'
  #colnames(lof.genes)[5] = "gene_category"
#   colnames(lof.genes)[7] = "genic_dup_60"

  lof.genes = lof %>% group_by(Cancer_type, Entrez) %>%
    summarise(n_samples=length(unique(Sample)),
              driver_alteration_types=paste(unique(type), collapse=","))
  lof.genes = lof.genes %>%
    left_join(geneInfo%>%select(entrez, symbol, cancer_type, complex), by=c("Entrez"="entrez")) %>%
    rename(gene_category=cancer_type)
  lof.genes$genic_dup_60 =  with( subset(dupByThresholds,threshold=="60"), ndups[ match(lof.genes$Entrez, entrez ) ] )
  lof.genes$genic_dup_20 =  with( subset(dupByThresholds,threshold=="20"), ndups[ match(lof.genes$Entrez, entrez ) ] )

  ix = which(is.na(lof.genes$genic_dup_20)); if(length(ix)>0) lof.genes$genic_dup_20[ix] = 0
  ix = which(is.na(lof.genes$genic_dup_60)); if(length(ix)>0) lof.genes$genic_dup_60[ix] = 0

  list(LOF_samples = lof.sample, LOF_genes = lof.genes)
}


## Create training and summary statistics
createDescribeTraining = function(cancer_type, tissue, res_dir, gains=T){
    if (!file.exists(res_dir)){
      dir.create(file.path(res_dir))
    }
    ## Get data from the database
    d = getTable(tb=cancer_type)
    ## Mark training observations before you scale
    d = markTrainNovelty(df=d, tissue=tissue, gains = gains)
    ## Save training and validation set without scaling (Feature refinement as well)
    training_ns = d %>% subset(!is.na(type))
    # training_ns = cleanFeatures(df=training_ns, ct=cancer_type)
    training_ns = cleanFeatures(df=training_ns, paste0(res_dir,"/config.txt"))
    training_ns = training_ns %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)
    if(length(grep("training_set_noScale.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(training_ns, file=paste(res_dir,"/training_set_noScale.Rdata", sep=""))
    }
    validation_ns = d %>% subset(is.na(type)) %>% select(-type)
#     validation_ns = cleanFeatures(df=validation_ns, ct=cancer_type)
    validation_ns = cleanFeatures(df=validation_ns, paste0(res_dir,"/config.txt"))
    validation_ns = validation_ns %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)
    if(length(grep("validation_set_noScale.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(validation_ns, file=paste(res_dir,"/validation_set_noScale.Rdata", sep=""))
    }

    ## Scale dataset
    d = getScaledTable(df = d)
    ## Clean and prepare it for the SVM
    # d = cleanFeatures(df=d, ct=cancer_type)
    d = cleanFeatures(df=d, paste0(res_dir,"/config.txt"))
    ## Feature refinement
    d = d %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)
    ## Separate training and prediction
    training = d %>% subset(!is.na(type))
    ## if training is not there save it
    if(length(grep("training_set.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(training, file=paste(res_dir,"/training_set.Rdata", sep=""))
    }
    validation = d %>% subset(is.na(type)) %>% select(-type)
    ## Save validation set
    if(length(grep("validation_set.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(validation, file=paste(res_dir,"/validation_set.Rdata", sep=""))
    }

    ## Get Correlation matrix
    varTypes = sapply(training, function(x) class(x)) %>% data.frame() %>% tibble::rownames_to_column()
    colnames(varTypes) = c("varType", "type")
    fvars = varTypes %>% subset(type=="factor") %>% .$varType
    t = training
    for (v in fvars){
        t[,v] = as.numeric(t[,v])
    }
    corMatrix = round(cor(t),2)
    upper = corMatrix
    upper[upper.tri(corMatrix)] = ""
    upper = as.data.frame(upper)
    write.xlsx(upper,file=paste0(res_dir, "/corMatrix.xlsx"), showNA = FALSE)

    ## Produce statistics for training
    stats = rbind(ddply(training_ns,.(type),numcolwise(quantile,na.rm = TRUE)) %>% mutate(statistic=c("min", "q1", "median", "q3", "max")),
          ddply(training_ns,.(type),numcolwise(mean,na.rm = TRUE)) %>% mutate(statistic="mean"))
    write.table(stats, file = paste(res_dir, "/training_summary_contFeatures_noScale.tsv", sep=""), row.names = F, quote = FALSE, sep="\t")
    stats = stats %>% gather(feature, value, -type, -statistic) %>% spread(statistic, value)
#     pdf(file=paste(res_dir, "/training_numerical_summary_noScale.pdf" , sep=""), h= 16, w=20)
#     p = ggplot(stats, aes(x=factor(feature), ymin = min, lower = q1, middle = median, upper = q3, ymax = max)) +geom_boxplot(stat = "identity") + scale_y_log10()
#     print(p)
#     dev.off()

    stats = rbind(ddply(training,.(type),numcolwise(quantile,na.rm = TRUE)) %>% mutate(statistic=c("min", "q1", "median", "q3", "max")),
                  ddply(training_ns,.(type),numcolwise(mean,na.rm = TRUE)) %>% mutate(statistic="mean"))
    write.table(stats, file = paste(res_dir, "/training_summary_contFeatures_Scale.tsv", sep=""), row.names = F, quote = FALSE, sep="\t")
    stats = stats %>% gather(feature, value, -type, -statistic) %>% spread(statistic, value)
#     pdf(file=paste(res_dir, "/training_numerical_summary_Scale.pdf" , sep=""), h= 16, w=20)
#     p = ggplot(stats, aes(x=factor(feature), ymin = min, lower = q1, middle = median, upper = q3, ymax = max)) + geom_boxplot(stat = "identity") + scale_y_log10()
#     print(p)
#     dev.off()
    ## This will be the same to bot scaled and unscaled training because categorical stay as they are independently of the scaling
    if(gains){
    stats = rbind(training_ns %>% count(type, Genic) %>% rename(categories=Genic, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="Genic"),
          training %>% count(type, memberofcomplex) %>% rename(categories=memberofcomplex, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="memberofcomplex"),
          training %>% count(type, WGD) %>% rename(categories=WGD, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="WGD"),
          training %>% count(type, hub) %>% rename(categories=hub, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="hub"),
          training %>% count(type, central) %>% rename(categories=central, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="central"),
          training %>% count(type, CNVGain) %>% rename(categories=CNVGain, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="CNVGain"),
          training %>% count(type, CNVLoss) %>% rename(categories=CNVLoss, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="CNVLoss"),
          training %>% count(type, ExpT_ME) %>% rename(categories=ExpT_ME, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_ME"),
          training %>% count(type, ExpT_HE) %>% rename(categories=ExpT_HE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_HE"),
          training %>% count(type, ExpT_LE) %>% rename(categories=ExpT_LE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_LE"),
          training %>% count(type, ExpT_NE) %>% rename(categories=ExpT_NE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_NE"),
          training %>% count(type, ExpT_NET) %>% rename(categories=ExpT_NET, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_NET"),
          training %>% count(type, old) %>% rename(categories=old, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="old"),
          training %>% count(type, young) %>% rename(categories=young, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="young"),
          training %>% count(type, luca) %>% rename(categories=luca, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="luca"),
          training %>% count(type, eukaryotes) %>% rename(categories=eukaryotes, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="eukaryotes"),
          training %>% count(type, metazoans) %>% rename(categories=metazoans, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="metazoans"),
          training %>% count(type, vertebrates) %>% rename(categories=vertebrates, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="vertebrates"),
          training %>% count(type, opisthokonts) %>% rename(categories=opisthokonts, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="opisthokonts"),
          training %>% count(type, mammals) %>% rename(categories=mammals, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="mammals"),
          training %>% count(type, selective) %>% rename(categories=selective, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="selective"),
          training %>% count(type, always_expressed) %>% rename(categories=always_expressed, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="always_expressed"),
          training %>% count(type, middle) %>% rename(categories=middle, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="middle"),
          training %>% count(type, one_tissue) %>% rename(categories=one_tissue, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="one_tissue"),
          training %>% count(type, never_expressed) %>% rename(categories=never_expressed, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="never_expressed")
    )
    }else{
      stats = rbind(training_ns %>% count(type, Genic) %>% rename(categories=Genic, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="Genic"),
                    training %>% count(type, memberofcomplex) %>% rename(categories=memberofcomplex, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="memberofcomplex"),
                    training %>% count(type, WGD) %>% rename(categories=WGD, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="WGD"),
                    training %>% count(type, hub) %>% rename(categories=hub, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="hub"),
                    training %>% count(type, central) %>% rename(categories=central, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="central"),
                    training %>% count(type, CNVLoss) %>% rename(categories=CNVLoss, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="CNVLoss"),
                    training %>% count(type, ExpT_ME) %>% rename(categories=ExpT_ME, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_ME"),
                    training %>% count(type, ExpT_HE) %>% rename(categories=ExpT_HE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_HE"),
                    training %>% count(type, ExpT_LE) %>% rename(categories=ExpT_LE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_LE"),
                    training %>% count(type, ExpT_NE) %>% rename(categories=ExpT_NE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_NE"),
                    training %>% count(type, ExpT_NET) %>% rename(categories=ExpT_NET, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_NET"),
                    training %>% count(type, old) %>% rename(categories=old, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="old"),
                    training %>% count(type, young) %>% rename(categories=young, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="young"),
                    training %>% count(type, luca) %>% rename(categories=luca, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="luca"),
                    training %>% count(type, eukaryotes) %>% rename(categories=eukaryotes, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="eukaryotes"),
                    training %>% count(type, metazoans) %>% rename(categories=metazoans, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="metazoans"),
                    training %>% count(type, vertebrates) %>% rename(categories=vertebrates, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="vertebrates"),
                    training %>% count(type, opisthokonts) %>% rename(categories=opisthokonts, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="opisthokonts"),
                    training %>% count(type, mammals) %>% rename(categories=mammals, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="mammals"),
                    training %>% count(type, selective) %>% rename(categories=selective, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="selective"),
                    training %>% count(type, always_expressed) %>% rename(categories=always_expressed, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="always_expressed"),
                    training %>% count(type, middle) %>% rename(categories=middle, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="middle"),
                    training %>% count(type, one_tissue) %>% rename(categories=one_tissue, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="one_tissue"),
                    training %>% count(type, never_expressed) %>% rename(categories=never_expressed, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="never_expressed")
      )

    }


    write.table(stats, file = paste(res_dir, "/training_summary_catFeatures.tsv", sep=""), row.names = F, quote = FALSE, sep="\t")
#     pdf(file=paste(res_dir, "/training_categorical_summary.pdf" , sep=""), h= 16, w=20)
#     p = ggplot(stats, aes(x=factor(feature), y=perc*100, fill=factor(categories))) + geom_bar(stat="identity") +
#         theme_boss() + theme(
#             axis.text.x = element_text(angle=90, size=12),
#             axis.text.y = element_text(size=12),
#             axis.title.x = element_blank(),
#             axis.title.y = element_text(size = 12),
#             legend.title = element_blank()
#         ) + scale_fill_brewer(palette="Set1") +
#         ylab("Number of occurences (%)")
#     print(p)
#     dev.off()

}

createDescribeTrainingVogel = function(cancer_type, res_dir, gains=T){

    if (!file.exists(res_dir)){
        dir.create(file.path(res_dir))
    }
    ## Get data from the database
    d = getTable(tb=cancer_type)
    ## Mark training observations before you scale
    d = markTrainVogel(df=d, gains = gains)
    ## Save training and validation set without scaling (Feature refinement as well)
    training_ns = d %>% subset(!is.na(type))
    # training_ns = cleanFeatures(df=training_ns, ct=cancer_type)
    training_ns = cleanFeatures(df=training_ns, paste0(res_dir,"/config.txt"))
    training_ns = training_ns %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)
    if(length(grep("training_set_noScale.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(training_ns, file=paste(res_dir,"/training_set_noScale.Rdata", sep=""))
    }
    validation_ns = d %>% subset(is.na(type)) %>% select(-type)
    #     validation_ns = cleanFeatures(df=validation_ns, ct=cancer_type)
    validation_ns = cleanFeatures(df=validation_ns, paste0(res_dir,"/config.txt"))
    validation_ns = validation_ns %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)
    if(length(grep("validation_set_noScale.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(validation_ns, file=paste(res_dir,"/validation_set_noScale.Rdata", sep=""))
    }

    ## Scale dataset
    d = getScaledTable(df = d)
    ## Clean and prepare it for the SVM
    # d = cleanFeatures(df=d, ct=cancer_type)
    d = cleanFeatures(df=d, paste0(res_dir,"/config.txt"))
    ## Feature refinement
    d = d %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)
    ## Separate training and prediction
    training = d %>% subset(!is.na(type))
    ## if training is not there save it
    if(length(grep("training_set.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(training, file=paste(res_dir,"/training_set.Rdata", sep=""))
    }
    validation = d %>% subset(is.na(type)) %>% select(-type)
    ## Save validation set
    if(length(grep("validation_set.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(validation, file=paste(res_dir,"/validation_set.Rdata", sep=""))
    }

    ## Get Correlation matrix
    varTypes = sapply(training, function(x) class(x)) %>% data.frame() %>% tibble::rownames_to_column()
    colnames(varTypes) = c("varType", "type")
    fvars = varTypes %>% subset(type=="factor") %>% .$varType
    t = training
    for (v in fvars){
        t[,v] = as.numeric(t[,v])
    }
    corMatrix = round(cor(t),2)
    upper = corMatrix
    upper[upper.tri(corMatrix)] = ""
    upper = as.data.frame(upper)
    write.xlsx(upper,file=paste0(res_dir, "/corMatrix.xlsx"), showNA = FALSE)

    ## Produce statistics for training
    stats = rbind(ddply(training_ns,.(type),numcolwise(quantile,na.rm = TRUE)) %>% mutate(statistic=c("min", "q1", "median", "q3", "max")),
                  ddply(training_ns,.(type),numcolwise(mean,na.rm = TRUE)) %>% mutate(statistic="mean"))
    write.table(stats, file = paste(res_dir, "/training_summary_contFeatures_noScale.tsv", sep=""), row.names = F, quote = FALSE, sep="\t")
    stats = stats %>% gather(feature, value, -type, -statistic) %>% spread(statistic, value)
    #     pdf(file=paste(res_dir, "/training_numerical_summary_noScale.pdf" , sep=""), h= 16, w=20)
    #     p = ggplot(stats, aes(x=factor(feature), ymin = min, lower = q1, middle = median, upper = q3, ymax = max)) +geom_boxplot(stat = "identity") + scale_y_log10()
    #     print(p)
    #     dev.off()

    stats = rbind(ddply(training,.(type),numcolwise(quantile,na.rm = TRUE)) %>% mutate(statistic=c("min", "q1", "median", "q3", "max")),
                  ddply(training_ns,.(type),numcolwise(mean,na.rm = TRUE)) %>% mutate(statistic="mean"))
    write.table(stats, file = paste(res_dir, "/training_summary_contFeatures_Scale.tsv", sep=""), row.names = F, quote = FALSE, sep="\t")
    stats = stats %>% gather(feature, value, -type, -statistic) %>% spread(statistic, value)
    #     pdf(file=paste(res_dir, "/training_numerical_summary_Scale.pdf" , sep=""), h= 16, w=20)
    #     p = ggplot(stats, aes(x=factor(feature), ymin = min, lower = q1, middle = median, upper = q3, ymax = max)) + geom_boxplot(stat = "identity") + scale_y_log10()
    #     print(p)
    #     dev.off()
    ## This will be the same to bot scaled and unscaled training because categorical stay as they are independently of the scaling
    if(gains){
        stats = rbind(training_ns %>% count(type, Genic) %>% rename(categories=Genic, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="Genic"),
                      training %>% count(type, memberofcomplex) %>% rename(categories=memberofcomplex, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="memberofcomplex"),
                      training %>% count(type, WGD) %>% rename(categories=WGD, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="WGD"),
                      training %>% count(type, hub) %>% rename(categories=hub, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="hub"),
                      training %>% count(type, central) %>% rename(categories=central, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="central"),
                      training %>% count(type, CNVGain) %>% rename(categories=CNVGain, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="CNVGain"),
                      training %>% count(type, CNVLoss) %>% rename(categories=CNVLoss, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="CNVLoss"),
                      training %>% count(type, ExpT_ME) %>% rename(categories=ExpT_ME, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_ME"),
                      training %>% count(type, ExpT_HE) %>% rename(categories=ExpT_HE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_HE"),
                      training %>% count(type, ExpT_LE) %>% rename(categories=ExpT_LE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_LE"),
                      training %>% count(type, ExpT_NE) %>% rename(categories=ExpT_NE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_NE"),
                      training %>% count(type, ExpT_NET) %>% rename(categories=ExpT_NET, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_NET"),
                      training %>% count(type, old) %>% rename(categories=old, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="old"),
                      training %>% count(type, young) %>% rename(categories=young, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="young"),
                      training %>% count(type, luca) %>% rename(categories=luca, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="luca"),
                      training %>% count(type, eukaryotes) %>% rename(categories=eukaryotes, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="eukaryotes"),
                      training %>% count(type, metazoans) %>% rename(categories=metazoans, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="metazoans"),
                      training %>% count(type, vertebrates) %>% rename(categories=vertebrates, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="vertebrates"),
                      training %>% count(type, opisthokonts) %>% rename(categories=opisthokonts, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="opisthokonts"),
                      training %>% count(type, mammals) %>% rename(categories=mammals, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="mammals"),
                      training %>% count(type, selective) %>% rename(categories=selective, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="selective"),
                      training %>% count(type, always_expressed) %>% rename(categories=always_expressed, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="always_expressed"),
                      training %>% count(type, middle) %>% rename(categories=middle, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="middle"),
                      training %>% count(type, one_tissue) %>% rename(categories=one_tissue, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="one_tissue"),
                      training %>% count(type, never_expressed) %>% rename(categories=never_expressed, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="never_expressed")
        )
    }else{
        stats = rbind(training_ns %>% count(type, Genic) %>% rename(categories=Genic, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="Genic"),
                      training %>% count(type, memberofcomplex) %>% rename(categories=memberofcomplex, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="memberofcomplex"),
                      training %>% count(type, WGD) %>% rename(categories=WGD, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="WGD"),
                      training %>% count(type, hub) %>% rename(categories=hub, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="hub"),
                      training %>% count(type, central) %>% rename(categories=central, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="central"),
                      training %>% count(type, CNVLoss) %>% rename(categories=CNVLoss, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="CNVLoss"),
                      training %>% count(type, ExpT_ME) %>% rename(categories=ExpT_ME, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_ME"),
                      training %>% count(type, ExpT_HE) %>% rename(categories=ExpT_HE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_HE"),
                      training %>% count(type, ExpT_LE) %>% rename(categories=ExpT_LE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_LE"),
                      training %>% count(type, ExpT_NE) %>% rename(categories=ExpT_NE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_NE"),
                      training %>% count(type, ExpT_NET) %>% rename(categories=ExpT_NET, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_NET"),
                      training %>% count(type, old) %>% rename(categories=old, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="old"),
                      training %>% count(type, young) %>% rename(categories=young, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="young"),
                      training %>% count(type, luca) %>% rename(categories=luca, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="luca"),
                      training %>% count(type, eukaryotes) %>% rename(categories=eukaryotes, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="eukaryotes"),
                      training %>% count(type, metazoans) %>% rename(categories=metazoans, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="metazoans"),
                      training %>% count(type, vertebrates) %>% rename(categories=vertebrates, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="vertebrates"),
                      training %>% count(type, opisthokonts) %>% rename(categories=opisthokonts, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="opisthokonts"),
                      training %>% count(type, mammals) %>% rename(categories=mammals, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="mammals"),
                      training %>% count(type, selective) %>% rename(categories=selective, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="selective"),
                      training %>% count(type, always_expressed) %>% rename(categories=always_expressed, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="always_expressed"),
                      training %>% count(type, middle) %>% rename(categories=middle, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="middle"),
                      training %>% count(type, one_tissue) %>% rename(categories=one_tissue, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="one_tissue"),
                      training %>% count(type, never_expressed) %>% rename(categories=never_expressed, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="never_expressed")
        )

    }


    write.table(stats, file = paste(res_dir, "/training_summary_catFeatures.tsv", sep=""), row.names = F, quote = FALSE, sep="\t")

}


createDescribeTrainingCGC = function(cancer_type, res_dir, gains=T, mysql=F, input_fn=NULL, means_sds=NULL){

    if (!file.exists(res_dir)){
        dir.create(file.path(res_dir))
    }
    ## Get data from the database
    if(mysql){
        d = getTable(tb=cancer_type)
    }else if(mysql==FALSE & !is.null(input_fn)){
        d = read.table(paste(input_fn), sep="\t", header = T)
        ## We need to make it factors
        cat("Converting features to factors...", "\n")
        fcols <- c("duplicated",
                   "WGD", "hub", "central", "CNVGain", "CNVLoss",
                   "ExpT_ME", "ExpT_HE", "ExpT_LE", "ExpT_NE",
                   "ExpT_NET", "old", "young", "luca", "eukaryotes",
                   "metazoans", "vertebrates", "opisthokonts",
                   "mammals", "primates", "selective", "always.expressed",
                   "middle", "one.tissue", "never.expressed")
        cols <- which(colnames(d) %in% fcols)
        for(i in cols){
            d[,i] = factor(d[,i], levels = c(0,1))
        }
    }else{
        stop("Input file argument is not present")
    }

    ## Mark training observations before you scale
    d = markTrainCGC(df=d, gains = gains)
    cat("Getting negative observations...", "\n")
    negatives_ns = d[["n"]]
    d = d[["tp"]]
    ## Save training and validation set without scaling (Feature refinement as well)
    training_ns = d %>% subset(!is.na(type))
    # training_ns = cleanFeatures(df=training_ns, ct=cancer_type)
    training_ns = cleanFeatures(df=training_ns, paste0(res_dir,"/config.txt"))
    #training_ns = training_ns %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)
    if(length(grep("training_set_noScale.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(training_ns, file=paste(res_dir,"/training_set_noScale.Rdata", sep=""))
    }
    validation_ns = d %>% subset(is.na(type)) %>% select(-type)
    #     validation_ns = cleanFeatures(df=validation_ns, ct=cancer_type)
    validation_ns = cleanFeatures(df=validation_ns, paste0(res_dir,"/config.txt"))
    #validation_ns = validation_ns %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)
    if(length(grep("validation_set_noScale.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(validation_ns, file=paste(res_dir,"/validation_set_noScale.Rdata", sep=""))
    }
    negatives_ns = cleanFeatures(df=negatives_ns, paste0(res_dir,"/config.txt"))
    #negatives_ns = negatives_ns %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)
    if(length(grep("negative_set_noScale.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(negatives_ns, file=paste(res_dir,"/negative_set_noScale.Rdata", sep=""))
    }


    ## Scale dataset
    ##load(means_sds)
    d = getScaledTable(df = d, means_sds = means_sds)
    mean_sd = d[["mean_sd"]]
    ## Save means and sds
    if(length(grep("mean_sd_of_training_validation_features.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(mean_sd, file=paste(res_dir,"/mean_sd_of_training_validation_features.Rdata", sep=""))
    }
    d = d[["df_scaled"]]
    ## Clean and prepare it for the SVM
    # d = cleanFeatures(df=d, ct=cancer_type)
    d = cleanFeatures(df=d, paste0(res_dir,"/config.txt"))
    ## Feature refinement
    #d = d %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)
    ## Separate training and prediction
    training = d %>% subset(!is.na(type))
    ## if training is not there save it
    if(length(grep("training_set.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(training, file=paste(res_dir,"/training_set.Rdata", sep=""))
    }
    validation = d %>% subset(is.na(type)) %>% select(-type)
    ## Save validation set
    if(length(grep("validation_set.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(validation, file=paste(res_dir,"/validation_set.Rdata", sep=""))
    }

    ## Scale the negatives using the mean and sd of the TP set
    negatives = negatives_ns
    for(i in colnames(negatives)){
        if(i %in% colnames(mean_sd)){
            scale_mean = mean_sd["means",grep(i, colnames(mean_sd))]
            scale_sd = mean_sd["sds",grep(i, colnames(mean_sd))]
            negatives[,i] = (negatives[,i] - as.numeric(scale_mean))/as.numeric(scale_sd)
        }
    }
    ## Save negative set
    if(length(grep("negative_set.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(negatives, file=paste(res_dir,"/negative_set.Rdata", sep=""))
    }

    ## Get Correlation matrix
    varTypes = sapply(training, function(x) class(x)) %>% data.frame() %>% tibble::rownames_to_column()
    colnames(varTypes) = c("varType", "type")
    fvars = varTypes %>% subset(type=="factor") %>% .$varType
    t = training
    for (v in fvars){
        t[,v] = as.numeric(t[,v])
    }
    corMatrix = round(cor(t),2)
    upper = corMatrix
    upper[upper.tri(corMatrix)] = ""
    upper = as.data.frame(upper)
    write.xlsx(upper,file=paste0(res_dir, "/corMatrix.xlsx"), showNA = FALSE)

    ## Produce statistics for training
    stats = rbind(ddply(training_ns,.(type),numcolwise(quantile,na.rm = TRUE)) %>% mutate(statistic=c("min", "q1", "median", "q3", "max")),
                  ddply(training_ns,.(type),numcolwise(mean,na.rm = TRUE)) %>% mutate(statistic="mean"))
    write.table(stats, file = paste(res_dir, "/training_summary_contFeatures_noScale.tsv", sep=""), row.names = F, quote = FALSE, sep="\t")
    stats = stats %>% gather(feature, value, -type, -statistic) %>% spread(statistic, value)
    #     pdf(file=paste(res_dir, "/training_numerical_summary_noScale.pdf" , sep=""), h= 16, w=20)
    #     p = ggplot(stats, aes(x=factor(feature), ymin = min, lower = q1, middle = median, upper = q3, ymax = max)) +geom_boxplot(stat = "identity") + scale_y_log10()
    #     print(p)
    #     dev.off()

    stats = rbind(ddply(training,.(type),numcolwise(quantile,na.rm = TRUE)) %>% mutate(statistic=c("min", "q1", "median", "q3", "max")),
                  ddply(training_ns,.(type),numcolwise(mean,na.rm = TRUE)) %>% mutate(statistic="mean"))
    write.table(stats, file = paste(res_dir, "/training_summary_contFeatures_Scale.tsv", sep=""), row.names = F, quote = FALSE, sep="\t")

    ## get counts for NAs (if any)
    nas = sapply(training_ns, function(x) (sum(is.na(x))/length(x))*100) %>% data.frame()
    colnames(nas) = "percentage_NAs"
    write.xlsx(nas, file = paste(res_dir, "/NAs_percentages.xlsx", sep=""), row.names = T, sheetName = "Training_ns")
    nas = sapply(training, function(x) (sum(is.na(x))/length(x))*100) %>% data.frame()
    colnames(nas) = "percentage_NAs"
    write.xlsx(nas, file = paste(res_dir, "/NAs_percentages.xlsx", sep=""), row.names = T, append = T, sheetName = "Training")
    nas = sapply(validation_ns, function(x) (sum(is.na(x))/length(x))*100) %>% data.frame()
    colnames(nas) = "percentage_NAs"
    write.xlsx(nas, file = paste(res_dir, "/NAs_percentages.xlsx", sep=""), row.names = T, append = T, sheetName = "Validation_ns")
    nas = sapply(validation, function(x) (sum(is.na(x))/length(x))*100) %>% data.frame()
    colnames(nas) = "percentage_NAs"
    write.xlsx(nas, file = paste(res_dir, "/NAs_percentages.xlsx", sep=""), row.names = T, append = T, sheetName = "Validation")
}

createDescribeTrainingTwoClasses = function(cancer_type, tissue, res_dir, gains=T){
    if (!file.exists(res_dir)){
        dir.create(file.path(res_dir))
    }
    ## Get data from the database
    d = getTable(tb=cancer_type)
    ## Mark training observations before you scale
    d = markTrainTwoClasses(df=d, tissue=tissue, gains = gains)
    ## Save training and validation set without scaling (Feature refinement as well)
    training_ns = d %>% subset(!is.na(type))
    # training_ns = cleanFeatures(df=training_ns, ct=cancer_type)
    training_ns = cleanFeatures(df=training_ns, paste0(res_dir,"/config.txt"))
    training_ns = training_ns %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)
    if(length(grep("training_set_noScale.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(training_ns, file=paste(res_dir,"/training_set_noScale.Rdata", sep=""))
    }
    validation_ns = d %>% subset(is.na(type)) %>% select(-type)
    #     validation_ns = cleanFeatures(df=validation_ns, ct=cancer_type)
    validation_ns = cleanFeatures(df=validation_ns, paste0(res_dir,"/config.txt"))
    validation_ns = validation_ns %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)
    if(length(grep("validation_set_noScale.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(validation_ns, file=paste(res_dir,"/validation_set_noScale.Rdata", sep=""))
    }

    ## Scale dataset
    d = getScaledTable(df = d)
    ## Clean and prepare it for the SVM
    # d = cleanFeatures(df=d, ct=cancer_type)
    d = cleanFeatures(df=d, paste0(res_dir,"/config.txt"))
    ## Feature refinement
    d = d %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)
    ## Separate training and prediction
    training = d %>% subset(!is.na(type))
    ## if training is not there save it
    if(length(grep("training_set.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(training, file=paste(res_dir,"/training_set.Rdata", sep=""))
    }
    validation = d %>% subset(is.na(type)) %>% select(-type)
    ## Save validation set
    if(length(grep("validation_set.Rdata",list.files(res_dir, recursive = F, full.names=F)))==0){
        save(validation, file=paste(res_dir,"/validation_set.Rdata", sep=""))
    }

    ## Produce statistics for training
    number_of_categories = training_ns %>% select(type) %>% unique %>% nrow
    stats = rbind(ddply(training_ns,.(type),numcolwise(quantile,na.rm = TRUE)) %>% mutate(statistic=rep(c("min", "q1", "median", "q3", "max"), number_of_categories)),
                  ddply(training_ns,.(type),numcolwise(mean,na.rm = TRUE)) %>% mutate(statistic=rep(c("mean"), number_of_categories)))
    write.table(stats, file = paste(res_dir, "/training_summary_contFeatures_noScale.tsv", sep=""), row.names = F, quote = FALSE, sep="\t")
    stats = stats %>% gather(feature, value, -type, -statistic) %>% spread(statistic, value)
    #     pdf(file=paste(res_dir, "/training_numerical_summary_noScale.pdf" , sep=""), h= 16, w=20)
    #     p = ggplot(stats, aes(x=factor(feature), ymin = min, lower = q1, middle = median, upper = q3, ymax = max)) +geom_boxplot(stat = "identity") + scale_y_log10()
    #     print(p)
    #     dev.off()

    stats = rbind(ddply(training,.(type),numcolwise(quantile,na.rm = TRUE)) %>% mutate(statistic=rep(c("min", "q1", "median", "q3", "max"), number_of_categories)),
                  ddply(training_ns,.(type),numcolwise(mean,na.rm = TRUE)) %>% mutate(statistic=rep(c("mean"), number_of_categories)))
    write.table(stats, file = paste(res_dir, "/training_summary_contFeatures_Scale.tsv", sep=""), row.names = F, quote = FALSE, sep="\t")
    stats = stats %>% gather(feature, value, -type, -statistic) %>% spread(statistic, value)
    #     pdf(file=paste(res_dir, "/training_numerical_summary_Scale.pdf" , sep=""), h= 16, w=20)
    #     p = ggplot(stats, aes(x=factor(feature), ymin = min, lower = q1, middle = median, upper = q3, ymax = max)) + geom_boxplot(stat = "identity") + scale_y_log10()
    #     print(p)
    #     dev.off()
    ## This will be the same to bot scaled and unscaled training because categorical stay as they are independently of the scaling
    if(gains){
        stats = rbind(training_ns %>% count(type, Genic) %>% rename(categories=Genic, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="Genic"),
                      training %>% count(type, memberofcomplex) %>% rename(categories=memberofcomplex, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="memberofcomplex"),
                      training %>% count(type, WGD) %>% rename(categories=WGD, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="WGD"),
                      training %>% count(type, hub) %>% rename(categories=hub, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="hub"),
                      training %>% count(type, central) %>% rename(categories=central, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="central"),
                      training %>% count(type, CNVGain) %>% rename(categories=CNVGain, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="CNVGain"),
                      training %>% count(type, CNVLoss) %>% rename(categories=CNVLoss, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="CNVLoss"),
                      training %>% count(type, ExpT_ME) %>% rename(categories=ExpT_ME, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_ME"),
                      training %>% count(type, ExpT_HE) %>% rename(categories=ExpT_HE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_HE"),
                      training %>% count(type, ExpT_LE) %>% rename(categories=ExpT_LE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_LE"),
                      training %>% count(type, ExpT_NE) %>% rename(categories=ExpT_NE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_NE"),
                      training %>% count(type, ExpT_NET) %>% rename(categories=ExpT_NET, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_NET"),
                      training %>% count(type, old) %>% rename(categories=old, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="old"),
                      training %>% count(type, young) %>% rename(categories=young, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="young"),
                      training %>% count(type, luca) %>% rename(categories=luca, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="luca"),
                      training %>% count(type, eukaryotes) %>% rename(categories=eukaryotes, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="eukaryotes"),
                      training %>% count(type, metazoans) %>% rename(categories=metazoans, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="metazoans"),
                      training %>% count(type, vertebrates) %>% rename(categories=vertebrates, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="vertebrates"),
                      training %>% count(type, opisthokonts) %>% rename(categories=opisthokonts, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="opisthokonts"),
                      training %>% count(type, mammals) %>% rename(categories=mammals, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="mammals"),
                      training %>% count(type, selective) %>% rename(categories=selective, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="selective"),
                      training %>% count(type, always_expressed) %>% rename(categories=always_expressed, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="always_expressed"),
                      training %>% count(type, middle) %>% rename(categories=middle, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="middle"),
                      training %>% count(type, one_tissue) %>% rename(categories=one_tissue, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="one_tissue"),
                      training %>% count(type, never_expressed) %>% rename(categories=never_expressed, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="never_expressed")
        )
    }else{
        stats = rbind(training_ns %>% count(type, Genic) %>% rename(categories=Genic, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="Genic"),
                      training %>% count(type, memberofcomplex) %>% rename(categories=memberofcomplex, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="memberofcomplex"),
                      training %>% count(type, WGD) %>% rename(categories=WGD, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="WGD"),
                      training %>% count(type, hub) %>% rename(categories=hub, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="hub"),
                      training %>% count(type, central) %>% rename(categories=central, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="central"),
                      training %>% count(type, CNVLoss) %>% rename(categories=CNVLoss, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="CNVLoss"),
                      training %>% count(type, ExpT_ME) %>% rename(categories=ExpT_ME, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_ME"),
                      training %>% count(type, ExpT_HE) %>% rename(categories=ExpT_HE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_HE"),
                      training %>% count(type, ExpT_LE) %>% rename(categories=ExpT_LE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_LE"),
                      training %>% count(type, ExpT_NE) %>% rename(categories=ExpT_NE, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_NE"),
                      training %>% count(type, ExpT_NET) %>% rename(categories=ExpT_NET, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="ExpT_NET"),
                      training %>% count(type, old) %>% rename(categories=old, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="old"),
                      training %>% count(type, young) %>% rename(categories=young, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="young"),
                      training %>% count(type, luca) %>% rename(categories=luca, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="luca"),
                      training %>% count(type, eukaryotes) %>% rename(categories=eukaryotes, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="eukaryotes"),
                      training %>% count(type, metazoans) %>% rename(categories=metazoans, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="metazoans"),
                      training %>% count(type, vertebrates) %>% rename(categories=vertebrates, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="vertebrates"),
                      training %>% count(type, opisthokonts) %>% rename(categories=opisthokonts, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="opisthokonts"),
                      training %>% count(type, mammals) %>% rename(categories=mammals, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="mammals"),
                      training %>% count(type, selective) %>% rename(categories=selective, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="selective"),
                      training %>% count(type, always_expressed) %>% rename(categories=always_expressed, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="always_expressed"),
                      training %>% count(type, middle) %>% rename(categories=middle, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="middle"),
                      training %>% count(type, one_tissue) %>% rename(categories=one_tissue, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="one_tissue"),
                      training %>% count(type, never_expressed) %>% rename(categories=never_expressed, occurrences=n) %>% mutate(perc=occurrences/sum(occurrences)) %>% mutate(feature="never_expressed")
        )

    }


    write.table(stats, file = paste(res_dir, "/training_summary_catFeatures.tsv", sep=""), row.names = F, quote = FALSE, sep="\t")
    #     pdf(file=paste(res_dir, "/training_categorical_summary.pdf" , sep=""), h= 16, w=20)
    #     p = ggplot(stats, aes(x=factor(feature), y=perc*100, fill=factor(categories))) + geom_bar(stat="identity") +
    #         theme_boss() + theme(
    #             axis.text.x = element_text(angle=90, size=12),
    #             axis.text.y = element_text(size=12),
    #             axis.title.x = element_blank(),
    #             axis.title.y = element_text(size = 12),
    #             legend.title = element_blank()
    #         ) + scale_fill_brewer(palette="Set1") +
    #         ylab("Number of occurences (%)")
    #     print(p)
    #     dev.off()

}

## Take the unscaled training and validation, combine them and infer mean and sd for each feature
## Then scale new incoming samples and predict using the best model from the corresponding cancer type
## Probably subset for what you consider driver first
## Example usage:
# scalePredictNewSamples(cancer_type_path = "/Volumes/mourikisa/novel_driver_prediction/ESCA/vogel_114_training",
#                        newSamples = oac, kern=c("linear", "radial", "polynomial","sigmoid"), mynu = c(0.15, 0.15, 0.05, 0.15), mygamma = c(0,0.0078125,0.015625,0.03125),
#                        mydegree = c(0,0,3,0), gains = F, fp_dir="/Volumes/mourikisa/data/NCG_false_positives.txt",
#                        res_dir = "/Volumes/mourikisa/novel_driver_prediction/OAC/prediction_using_best_ESCA/vogel_114_training", install_NCG = F)
scalePredictNewSamples = function(cancer_type_path, newSamples, kern=c("linear", "radial", "polynomial", "sigmoid"), mynu=c(0,0,0,0),
                                  mygamma=c(0,0,0,0), mydegree=c(0,0,0,0), gains=T, fp_dir="/mnt/lustre/users/k1469280/mourikisa/data/NCG_false_positives.txt", res_dir,
                                  install_NCG=F){

    models = data.frame(kern=kern, mynu=mynu, mygamma=mygamma, mydegree=mydegree)

    dev_mode()
    if(install_NCG) install(NCG_path)
    library(ncglib)
    geneInfo=get_geneinfo(version="NCG5")
    cancerGenes = get_cancer_genes(version="NCG5")
    dev_mode()

    ## Fix gene info table from NCG
    geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
    ## Get a cancer gene with all the associated primary sites and cancer sites
    cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
        group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
                                       cancer_site=paste(unique(cancer_site), collapse=",")) %>%
        ungroup
    geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))

    ## Get false positive genes
    false_positive_genes <- read.table(fp_dir, header = T, sep = "\t")
    false_positive_genes <- false_positive_genes %>% select(entrez, symbol)
    ## For now lets include only driver alterations - based on our definition
    if(gains){
        newSamples = newSamples %>% subset(no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | (CNVGain==1 & Copy_number>=4) |
                                                             (CNVLoss==1 & Copy_number==0) | ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) & !(Entrez%in%false_positive_genes$entrez))
    }else{
        newSamples = newSamples %>% subset(no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 |
                                                             (CNVLoss==1 & Copy_number==0) | ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) & !(Entrez%in%false_positive_genes$entrez))
    }

    ## training_ns
    load(paste0(cancer_type_path, "/training_set_noScale.Rdata"))
    ## validation_ns
    load(paste0(cancer_type_path, "/validation_set_noScale.Rdata"))
    ns = rbind(training_ns%>%select(-type), validation_ns)
    ## training
    load(paste0(cancer_type_path, "/training_set.Rdata"))
    ## validation
    load(paste0(cancer_type_path, "/validation_set.Rdata"))

    ## Get mean and sd of the corresponding cancer type
    stats = as.data.frame(lapply(ns, mean)) %>% gather(feature, mean) %>% subset(!is.na(mean)) %>% left_join(
        as.data.frame(lapply(ns, sd)) %>% gather(feature, sd) %>% subset(!is.na(sd)), by=c("feature")
    )

    ## Scale using the mean and sd of the corresponding cancer type
    x = newSamples%>%select(one_of(as.character(stats$feature)))
    x = as.list(x)
    y = as.list(stats$mean)
    names(y) = stats$feature
    z = as.list(stats$sd)
    names(z) = stats$feature

    ## Add factor columns etc and reorder the columns using the order of the unscaled data frame
    newSamples_scaled = cbind(data.frame(mapply( function(x,y,z) (x-y)/z, x, y, z, SIMPLIFY = F )), newSamples[,!(names(newSamples) %in% stats$feature)])
    newSamples_scaled = newSamples_scaled%>%select(one_of(names(newSamples)))

    ## Prepare new samples for prediction
    newSamples_scaled = cleanFeatures(df=newSamples_scaled, paste0(cancer_type_path,"/config.txt"))
    ## Feature refinement
    newSamples_scaled = newSamples_scaled %>% select(-no_NSI_muts, -private_domains, -commondomains, -High, -Low, -Medium, -NotExpressed)

    ## What features should be excluded
    config = read.table(paste0(cancer_type_path,"/config.txt"), header = F)
    config = config[10:nrow(config), 1]
    s = c("cancer_type", "sample", "entrez", "no_TRUNC_muts", "no_NTDam_muts",
          "no_GOF_muts", "Copy_number", "CNVLoss", "CNVGain")
    s = s[!(s%in%config)]

    newSamples = newSamples %>% rename(cancer_type=Cancer_type,
                                       sample=Sample, entrez=Entrez) %>%
        select(one_of(s))

    prediction_table = NULL
    training_size = nrow(training)
    training_genes = training%>%tibble::rownames_to_column()%>%separate(rowname, into=c("Cancer_type", "sample", "entrez"), sep="\\.")%>%select(entrez)%>%unique%>%.$entrez
    training_samples = training%>%tibble::rownames_to_column()%>%separate(rowname, into=c("Cancer_type", "sample", "entrez"), sep="\\.")%>%select(sample)%>%unique%>%.$sample
    prediction_size = nrow(newSamples_scaled)
    prediction_genes = newSamples_scaled%>%tibble::rownames_to_column()%>%separate(rowname, into=c("Cancer_type", "sample", "entrez"), sep="\\.")%>%select(entrez)%>%unique%>%.$entrez
    prediction_samples = newSamples_scaled%>%tibble::rownames_to_column()%>%separate(rowname, into=c("Cancer_type", "sample", "entrez"), sep="\\.")%>%select(sample)%>%unique%>%.$sample

    for (i in 1:nrow(models)){
        ## Get the model
        if (models$kern[i]=="linear"){
            svm.model <- svm(type ~ ., data = training, kernel=models$kern[i], type="one-classification", scale=FALSE, nu=models$mynu[i])
        }else if (models$kern[i]=="radial"){
            svm.model <- svm(type ~ ., data = training, kernel=models$kern[i], type="one-classification", scale=FALSE, gamma=models$mygamma[i], nu=models$mynu[i])
        }else if (models$kern[i]=="polynomial"){
            svm.model <- svm(type ~ ., data = training, kernel=models$kern[i], type="one-classification", scale=FALSE, gamma=models$mygamma[i], nu=models$mynu[i], degree=models$mydegree[i])
        }else if (models$kern[i]=="sigmoid"){
            svm.model <- svm(type ~ ., data = training, kernel=models$kern[i], type="one-classification", scale=FALSE, gamma=models$mygamma[i], nu=models$mynu[i])
        }

        ## Asses training sensitivity
        trainset = training %>% select(type)
        trainset$fitted = svm.model$fitted
        trainset = trainset %>% select(type, fitted)
        trainset = trainset %>% mutate(f=ifelse(fitted==TRUE, "C", "N"))
        trainset$type = factor(as.character(trainset$type), levels = c("C", "N"))
        trainset$f = factor(as.character(trainset$f), levels = c("C", "N"))
        ## Confusion matrix
        cM <- confusionMatrix(trainset$f,trainset$type, positive = c("C"))
        training_sensitivity = cM$byClass["Sensitivity"]

        ## Predict on the new samples using the model
        predictions = predict(svm.model, newSamples_scaled, decision.values = TRUE)

        ## Add probability from Platt
        predictions_platt = plattScaling(svm.model, predictions)
        colnames(predictions_platt) = c("dv", "label", "prob_platt")
        preds = predictions_platt %>% tibble::rownames_to_column()

        ## Take information from gene info
        preds = preds %>%
            separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>%
            mutate(entrez=as.numeric(entrez)) %>%
            left_join(geneInfo%>%rename(gene_type=cancer_type), by=c("entrez"))

        ## Join also with systems level properties from the unscaled data frame
        preds = preds %>% left_join(newSamples, by=c("cancer_type", "sample", "entrez"))

        ## Save preds
        if (models$kern[i]=="linear"){
            analysis = paste(models$kern[i], ".", models$mynu[i], sep="")
            save(svm.model, file=paste(res_dir,"/svm.model_",analysis, ".Rdata", sep=""))
            save(preds, file=paste(res_dir,"/predictions_",analysis, ".Rdata", sep=""))
        }else if (models$kern[i]=="radial"){
            analysis = paste(models$kern[i], ".", models$mynu[i], ".", models$mygamma[i], sep="")
            save(svm.model, file=paste(res_dir,"/svm.model_",analysis, ".Rdata", sep=""))
            save(preds, file=paste(res_dir, "/predictions_",analysis, ".Rdata", sep=""))
        }else if (models$kern[i]=="polynomial"){
            analysis = paste(models$kern[i], ".", models$mynu[i], ".", models$mygamma[i], ".", models$mydegree[i], sep="")
            save(svm.model, file=paste(res_dir,"/svm.model_",analysis, ".Rdata", sep=""))
            save(preds, file=paste(res_dir, "/predictions_",analysis, ".Rdata", sep=""))
        }else if (models$kern[i]=="sigmoid"){
            analysis = paste(models$kern[i], ".", models$mynu[i], ".", models$mygamma[i], sep="")
            save(svm.model, file=paste(res_dir,"/svm.model_",analysis, ".Rdata", sep=""))
            save(preds, file=paste(res_dir, "/predictions_",analysis, ".Rdata", sep=""))
        }

        ## Get some stats for the predictions
        df = data.frame(analysis=analysis, nu=models$mynu[i],
                        gamma=models$mygamma[i], degree=models$mydegree[i],
                        novel_drivers_entries=preds%>%subset(label==TRUE)%>%nrow,
                        novel_drivers_genes=preds%>%subset(label==TRUE)%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_genes_inNCG5 = preds%>%subset(label==TRUE & (gene_type=="cgc" | gene_type=="can"))%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_genes_rst = preds%>%subset(label==TRUE & gene_type=="rst")%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_entries_platt=preds%>%subset(label==TRUE & prob_platt >0.95)%>%nrow,
                        novel_drivers_genes_platt=preds%>%subset(label==TRUE & prob_platt >0.95)%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_in_samples_mutated = preds%>%subset(label==TRUE & prob_platt >0.95)%>%select(sample)%>%unique%>%nrow,
                        novel_drivers_genes_platt_inNCG5=preds%>%subset(label==TRUE & prob_platt >0.95 & (gene_type=="cgc" | gene_type=="can"))%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_genes_platt_rst=preds%>%subset(label==TRUE & prob_platt >0.95 & gene_type=="rst")%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_genes_platt_rst_in_samples_mutated=preds%>%subset(label==TRUE & prob_platt >0.95 & gene_type=="rst")%>%select(sample)%>%unique%>%nrow
        )

        ## Read the cv_stats to get test sensitivity
        cv_stats = read.table(file=paste0(cancer_type_path, "/cv_stat_summary.tsv"), header = T, sep="\t")
        cv_stats = cv_stats %>% subset(kernel==models$kern[i] & type=="Sensitivity" & set=="test") %>% select(analysis, min, q1, median, mean, q3, max)
        df = df %>% left_join(cv_stats, by=c("analysis"))
        prediction_table = rbind(prediction_table, df)
    }


    ## Excel configuration
    require(xlsx)
    wb<-createWorkbook(type="xlsx")

    TEXT_STYLE           <- CellStyle(wb) + Font(wb,  heightInPoints=14)
    TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
        Alignment(rotation=45, horizontal="ALIGN_CENTER") +
        Border(color="black", position=c("TOP", "LEFT","BOTTOM", "RIGHT"), pen=rep("BORDER_THIN",4))


    write_text_in_cell = function(sheet, rowIndex, text, textStyle){
        rows <-createRow( sheet, rowIndex=rowIndex )
        sheetText <-createCell( rows, colIndex=1 )
        setCellValue( sheetText[[1,1]], text )
        setCellStyle( sheetText[[1,1]], textStyle )
    }

    sheet <- createSheet(wb, sheetName = "predictions")

    text = with( prediction_table, paste0("TRAINING SET. Number_of_genes = ", length(training_genes), " ; Samples = ", length(training_samples), " ; Total_entries = ", training_size))
    write_text_in_cell(sheet, 1, text, TEXT_STYLE)

    text = with( prediction_table, paste0("TEST SET.     Number_of_genes = ", length(prediction_genes), " ; Samples = ", length(prediction_samples), " ; Total_entries = ", prediction_size))
    write_text_in_cell(sheet, 2, text, TEXT_STYLE)

    to_remove = c("training_genes","training_samples",'training_set_size','prediction_genes','prediction_samples','prediction_set_size')

    rownames(prediction_table) = NULL

    addDataFrame(prediction_table[, !colnames(prediction_table)%in%to_remove ], sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, row.names = F)

    saveWorkbook(wb,file=paste(res_dir,"/predictions.xlsx", sep=""))

}

## Get the table for each cancer type from the MySQL database
getTable = function(tb="LAML",user="db_user",pass="db_pass",db="tcga_31cts_prepared",...){

    cat(paste("Retrieving data from ",db," MySQL database...",sep=""), "\n")
    ## Get in which node you are
    node = Sys.info()[["nodename"]]
    temp2 <- gregexpr("[0-9]+", node)
    node = as.numeric(unique(unlist(regmatches(node, temp2))))
    if (node==42 | node==43){
        mydb = dbConnect(MySQL(), user=user,
                         password=pass,
                         dbname=db,
                         host='ibnode041')
    }else if (node==41){
        mydb = dbConnect(MySQL(), user='root',
                         dbname=db)
    }

    #tables = dbListTables(mydb)
    rs = dbSendQuery(mydb, paste("select * from ", tb, sep=""))
    data = fetch(rs, n=-1)

    ## Fix factors and NAs in vogel
    ## Keep in mind that factor levels may be different from the pancancer data if all levels are not represented in this cancer type
    fcols <- c("Genic", "memberofcomplex",
               "WGD", "hub", "central", "CNVGain", "CNVLoss",
               "ExpT_ME", "ExpT_HE", "ExpT_LE", "ExpT_NE",
               "ExpT_NET", "old", "young", "luca", "eukaryotes",
               "metazoans", "vertebrates", "opisthokonts",
               "mammals", "primates","selective", "always_expressed",
               "middle", "one_tissue", "never_expressed")
    cols <- which(colnames(data) %in% fcols)
    for(i in cols){
        ## Check unique values to see the levels
        ## If there is only one level explicilty define at least 2
        if (length(unique(data[,i]))<2){
            data[,i] = factor(data[,i], levels = c("0","1"))
        }else{
            data[,i] = factor(data[,i])
        }

    }

    return(data)
}

## Mark the training set because after scaling you cannot filter (e.g no_TRUNC_muts!=0)
markTrainNovelty <- function(df, tissue=NULL, fp_dir="/mnt/lustre/users/k1469280/mourikisa/data/NCG_false_positives.txt", gains=T,...){

    if (is.null(tissue)){
        stop("Please provide tissue to get training set")
    }

  ## To take trainig set based on the candidate genes for each cancer type
  ## Load geneInfo.Rdata and cancerGenes.Rdata from /mnt/lustre/users/k1469280/mourikisa/data in athena
  load("/mnt/lustre/users/k1469280/mourikisa/data/geneInfoNCG5.Rdata")
  load("/mnt/lustre/users/k1469280/mourikisa/data/cancerGenesNCG5.Rdata")

  ## Get information for cgc/cans and primary sites
  ## Fix gene info table from NCG
  geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
  ## Get a cancer gene with all the associated primary sites and cancer sites
  cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
    group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
                                   cancer_site=paste(unique(cancer_site), collapse=",")) %>%
    ungroup
  geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez")) %>%
      rename(gene_type=cancer_type)

  ## Get false positive genes
  false_positive_genes <- read.table(fp_dir, header = T, sep = "\t")
  false_positive_genes <- false_positive_genes %>% select(entrez, symbol)

  ## Get primary site information for the df
  df = df %>% left_join(geneInfo, by=c("Entrez"="entrez"))

  ## Training set extraction
  cat("Marking genes for training...", "\n")
  #df = df %>% mutate(TP=ifelse((vogel=="Vog.Oncogene" | vogel=="Vog.TS") & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | CNVGain==1 | CNVLoss==1 | ExpT_NET==1), T, F))
  #if( df%>%grepl(tissue, primary_site)%>%nrow==0 ) stop("Tissue not present")

  if(gains){
    df = df %>% mutate(TP=ifelse((grepl(tissue, primary_site)) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | (CNVGain==1 & Copy_number>=4) |
                                    (CNVLoss==1 & Copy_number==0) | ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0))) & !(Entrez%in%false_positive_genes$entrez), T, F))
  } else {
    df = df %>% mutate(TP=ifelse((grepl(tissue, primary_site)) &
                                   (no_TRUNC_muts!=0 |
                                    no_NTDam_muts!=0 |
                                    no_GOF_muts!=0  |
                                   (CNVLoss==1 & Copy_number==0) |
                                  ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0))) & !(Entrez%in%false_positive_genes$entrez), T, F))

  }
  #df = df %>% mutate(TN=ifelse((Entrez%in%false_positive_genes$entrez) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | CNVGain==1 | CNVLoss==1 | ExpT_NET==1), T, F))
  #df = df %>% mutate(TN2=ifelse((vogel=="Vog.Oncogene" | vogel=="Vog.TS") & (no_TRUNC_muts==0 & no_NTDam_muts==0 & no_GOF_muts==0 & CNVGain==0 & CNVLoss==0 & ExpT_NET==0), T, F))
  cat("Adding type annotation for training...", "\n")
  #df = df %>% mutate(type=ifelse(TP==T & vogel=="Vog.Oncogene", "O", ifelse(TP==T & vogel=="Vog.TS", "T",
  #                                                                          ifelse(TN==T, "N", NA))))

  ## Some of the cancer specific genes are also in the true negatives list (in NCG we only show a warning)
  #df = df %>% mutate(type=ifelse((TP==T & TN==F), "C", ifelse((TN==T & TP==F) | (TN==T & TP==T), "N", NA)))
  df = df %>% mutate(type=ifelse(TP==T, "C", NA))

  df = df %>% select(-TP)
  ## For cancer-specific training set
  df = df %>% select(-gene_type, -primary_site, -cancer_site)
  df$type = as.factor(df$type)
  ## Exclude genes without any alterations from training & prediction set
  if(gains){
    df = df %>% subset(no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | (CNVGain==1 & Copy_number>=4) |
                         (CNVLoss==1 & Copy_number==0) | ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) & !(Entrez%in%false_positive_genes$entrez))
  }else{
    df = df %>% subset(no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 |
                       (CNVLoss==1 & Copy_number==0) | ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) & !(Entrez%in%false_positive_genes$entrez))
  }
  return(df)
}

markTrainVogel <- function(df, fp_dir="/mnt/lustre/users/k1469280/mourikisa/data/NCG_false_positives.txt", gains=T,...){

    ## To take trainig set based on the candidate genes for each cancer type
    ## Load geneInfo.Rdata and cancerGenes.Rdata from /mnt/lustre/users/k1469280/mourikisa/data in athena
    load("/mnt/lustre/users/k1469280/mourikisa/data/geneInfoNCG5.Rdata")
    load("/mnt/lustre/users/k1469280/mourikisa/data/cancerGenesNCG5.Rdata")

    ## Get information for cgc/cans and primary sites
    ## Fix gene info table from NCG
    geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
    ## Get a cancer gene with all the associated primary sites and cancer sites
    cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
        group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
                                       cancer_site=paste(unique(cancer_site), collapse=",")) %>%
        ungroup
    geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez")) %>%
        rename(gene_type=cancer_type)

    ## Get false positive genes
    false_positive_genes <- read.table(fp_dir, header = T, sep = "\t")
    false_positive_genes <- false_positive_genes %>% select(entrez, symbol)

    ## Get primary site information for the df
    df = df %>% left_join(geneInfo, by=c("Entrez"="entrez"))

    ## Training set extraction
    cat("Marking genes for training...", "\n")
    #df = df %>% mutate(TP=ifelse((vogel=="Vog.Oncogene" | vogel=="Vog.TS") & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | CNVGain==1 | CNVLoss==1 | ExpT_NET==1), T, F))
    #if( df%>%grepl(tissue, primary_site)%>%nrow==0 ) stop("Tissue not present")

    if(gains){
        df = df %>% mutate(TP=ifelse((vogel=="Vog.Oncogene" | vogel=="Vog.TS") & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | (CNVGain==1 & Copy_number>=4) |
                                                                          (CNVLoss==1 & Copy_number==0) | ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0))) & !(Entrez%in%false_positive_genes$entrez), T, F))
    } else {
        df = df %>% mutate(TP=ifelse((vogel=="Vog.Oncogene" | vogel=="Vog.TS") &
                                         (no_TRUNC_muts!=0 |
                                              no_NTDam_muts!=0 |
                                              no_GOF_muts!=0  |
                                              (CNVLoss==1 & Copy_number==0) |
                                              ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0))) & !(Entrez%in%false_positive_genes$entrez), T, F))

    }
    #df = df %>% mutate(TN=ifelse((Entrez%in%false_positive_genes$entrez) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | CNVGain==1 | CNVLoss==1 | ExpT_NET==1), T, F))
    #df = df %>% mutate(TN2=ifelse((vogel=="Vog.Oncogene" | vogel=="Vog.TS") & (no_TRUNC_muts==0 & no_NTDam_muts==0 & no_GOF_muts==0 & CNVGain==0 & CNVLoss==0 & ExpT_NET==0), T, F))
    cat("Adding type annotation for training...", "\n")
    #df = df %>% mutate(type=ifelse(TP==T & vogel=="Vog.Oncogene", "O", ifelse(TP==T & vogel=="Vog.TS", "T",
    #                                                                          ifelse(TN==T, "N", NA))))

    ## Some of the cancer specific genes are also in the true negatives list (in NCG we only show a warning)
    #df = df %>% mutate(type=ifelse((TP==T & TN==F), "C", ifelse((TN==T & TP==F) | (TN==T & TP==T), "N", NA)))
    df = df %>% mutate(type=ifelse(TP==T, "C", NA))

    df = df %>% select(-TP)
    ## For cancer-specific training set
    df = df %>% select(-gene_type, -primary_site, -cancer_site)
    df$type = as.factor(df$type)
    ## Exclude genes without any alterations from training & prediction set
    if(gains){
        df = df %>% subset(no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | (CNVGain==1 & Copy_number>=4) |
                               (CNVLoss==1 & Copy_number==0) | ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) & !(Entrez%in%false_positive_genes$entrez))
    }else{
        df = df %>% subset(no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 |
                               (CNVLoss==1 & Copy_number==0) | ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) & !(Entrez%in%false_positive_genes$entrez))
    }
    return(df)
}

markTrainCGC <- function(df, fp_dir="/mnt/lustre/users/k1469280/mourikisa/data/NCG_false_positives.txt",
                         geneInfo_fn="/mnt/lustre/users/k1469280/mourikisa/data/geneInfoNCG5.Rdata",
                         cancerGenes_fn="/mnt/lustre/users/k1469280/mourikisa/data/cancerGenesNCG5.Rdata",
                         gains=T, ...){

    ## To take trainig set based on the candidate genes for each cancer type
    ## Load geneInfo.Rdata and cancerGenes.Rdata from /mnt/lustre/users/k1469280/mourikisa/data in athena
    load(geneInfo_fn)
    load(cancerGenes_fn)

    ## Get information for cgc/cans and primary sites
    ## Fix gene info table from NCG
    geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
    ## Get a cancer gene with all the associated primary sites and cancer sites
    cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
        group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
                                       cancer_site=paste(unique(cancer_site), collapse=",")) %>%
        ungroup
    geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez")) %>%
        rename(gene_type=cancer_type)

    ## Get false positive genes
    false_positive_genes <- read.table(fp_dir, header = T, sep = "\t")
    false_positive_genes <- false_positive_genes %>% select(entrez, symbol)

    ## Get primary site information for the df
    df = df %>% left_join(geneInfo, by=c("Entrez"="entrez"))

    ## Training set extraction
    cat("Marking genes for training...", "\n")
    #df = df %>% mutate(TP=ifelse((vogel=="Vog.Oncogene" | vogel=="Vog.TS") & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | CNVGain==1 | CNVLoss==1 | ExpT_NET==1), T, F))
    #if( df%>%grepl(tissue, primary_site)%>%nrow==0 ) stop("Tissue not present")

    if(gains){
        df = df %>% mutate(TP=ifelse((gene_type=="cgc") & (no_TRUNC_muts!=0 |
                                                               no_NTDam_muts!=0 |
                                                               no_GOF_muts!=0 |
                                                               (CNVGain==1) |
                                                               (Copy_number==0) |
                                                               ((Copy_number==1) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) |
                                                               BND!=0 |
                                                               INS!=0 |
                                                               INV!=0) & !(Entrez%in%false_positive_genes$entrez), T, F))
    } else {
        df = df %>% mutate(TP=ifelse((gene_type=="cgc") &
                                         (no_TRUNC_muts!=0 |
                                              no_NTDam_muts!=0 |
                                              no_GOF_muts!=0 |
                                              (Copy_number==0) |
                                              ((Copy_number==1) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) |
                                              BND!=0 |
                                              INS!=0 |
                                              INV!=0) & !(Entrez%in%false_positive_genes$entrez), T, F))

    }
    #df = df %>% mutate(TN=ifelse((Entrez%in%false_positive_genes$entrez) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | CNVGain==1 | CNVLoss==1 | ExpT_NET==1), T, F))
    #df = df %>% mutate(TN2=ifelse((vogel=="Vog.Oncogene" | vogel=="Vog.TS") & (no_TRUNC_muts==0 & no_NTDam_muts==0 & no_GOF_muts==0 & CNVGain==0 & CNVLoss==0 & ExpT_NET==0), T, F))
    cat("Adding type annotation for training...", "\n")
    #df = df %>% mutate(type=ifelse(TP==T & vogel=="Vog.Oncogene", "O", ifelse(TP==T & vogel=="Vog.TS", "T",
    #                                                                          ifelse(TN==T, "N", NA))))

    ## Some of the cancer specific genes are also in the true negatives list (in NCG we only show a warning)
    #df = df %>% mutate(type=ifelse((TP==T & TN==F), "C", ifelse((TN==T & TP==F) | (TN==T & TP==T), "N", NA)))
    df = df %>% mutate(type=ifelse(TP==T, "C", NA))

    df = df %>% select(-TP)
    ## For cancer-specific training set
    df = df %>% select(-gene_type, -primary_site, -cancer_site)
    df$type = as.factor(df$type)

    ## Get negatives
    negatives = df %>% subset((no_TRUNC_muts==0 &
                                   no_NTDam_muts==0 &
                                   no_GOF_muts==0 &
                                   CNVGain==0 &
                            CNVLoss==0) | Entrez%in%false_positive_genes$entrez)

    ## Exclude genes without any alterations from training & prediction set  (I now added the exclusion of the NCG FPs as an option by default they are not excluded)
    if(gains){
        # df = df %>% subset((no_TRUNC_muts!=0 |
        #                        no_NTDam_muts!=0 |
        #                        no_GOF_muts!=0 |
        #                        (CNVGain==1) |
        #                        (Copy_number==0) |
        #                        ((Copy_number==1) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) |
        #                        BND!=0 |
        #                        INS!=0 |
        #                        INV!=0) & !(Entrez%in%false_positive_genes$entrez))
        df = df %>% subset((no_TRUNC_muts!=0 |
                                no_NTDam_muts!=0 |
                                no_GOF_muts!=0 |
                                (CNVGain==1) |
                                (Copy_number==0) |
                                ((Copy_number==1) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) |
                                BND!=0 |
                                INS!=0 |
                                INV!=0) )
    }else{
        # df = df %>% subset((no_TRUNC_muts!=0 |
        #                        no_NTDam_muts!=0 |
        #                        no_GOF_muts!=0 |
        #                        (Copy_number==0) |
        #                        ((Copy_number==1) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) |
        #                        BND!=0 |
        #                        INS!=0 |
        #                        INV!=0) & !(Entrez%in%false_positive_genes$entrez))
        df = df %>% subset((no_TRUNC_muts!=0 |
                                no_NTDam_muts!=0 |
                                no_GOF_muts!=0 |
                                (Copy_number==0) |
                                ((Copy_number==1) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) |
                                BND!=0 |
                                INS!=0 |
                                INV!=0) )
    }
    return(list(tp=df, n=negatives))
}

markTrainTwoClasses <- function(df, tissue=NULL, fp_dir="/mnt/lustre/users/k1469280/mourikisa/data/NCG_false_positives.txt", gains=T,...){

    if (is.null(tissue)){
        stop("Please provide tissue to get training set")
    }

    ## To take trainig set based on the candidate genes for each cancer type
    ## Load geneInfo.Rdata and cancerGenes.Rdata from /mnt/lustre/users/k1469280/mourikisa/data in athena
    load("/mnt/lustre/users/k1469280/mourikisa/data/geneInfoNCG5.Rdata")
    load("/mnt/lustre/users/k1469280/mourikisa/data/cancerGenesNCG5.Rdata")

    ## Get information for cgc/cans and primary sites
    ## Fix gene info table from NCG
    geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
    ## Get a cancer gene with all the associated primary sites and cancer sites
    cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
        group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
                                       cancer_site=paste(unique(cancer_site), collapse=",")) %>%
        ungroup
    geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))

    ## Get false positive genes
    false_positive_genes <- read.table(fp_dir, header = T, sep = "\t")
    false_positive_genes <- false_positive_genes %>% select(entrez, symbol)

    ## Get primary site information for the df
    df = df %>% left_join(geneInfo, by=c("Entrez"="entrez"))

    ## Training set extraction
    cat("Marking genes for training...", "\n")
    #df = df %>% mutate(TP=ifelse((vogel=="Vog.Oncogene" | vogel=="Vog.TS") & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | CNVGain==1 | CNVLoss==1 | ExpT_NET==1), T, F))
    #if( df%>%grepl(tissue, primary_site)%>%nrow==0 ) stop("Tissue not present")

    if(gains){
        df = df %>% mutate(TP=ifelse((grepl(tissue, primary_site)) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | (CNVGain==1 & Copy_number>=4) |
                                                                          (CNVLoss==1 & Copy_number==0) | ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0))) & !(Entrez%in%false_positive_genes$entrez), T, F))
    } else {
        df = df %>% mutate(TP=ifelse((grepl(tissue, primary_site)) &
                                         (no_TRUNC_muts!=0 |
                                              no_NTDam_muts!=0 |
                                              no_GOF_muts!=0  |
                                              (CNVLoss==1 & Copy_number==0) |
                                              ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0))) & !(Entrez%in%false_positive_genes$entrez), T, F))

    }
    #df = df %>% mutate(TN=ifelse((Entrez%in%false_positive_genes$entrez) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | CNVGain==1 | CNVLoss==1 | ExpT_NET==1), T, F))
    #df = df %>% mutate(TN2=ifelse((vogel=="Vog.Oncogene" | vogel=="Vog.TS") & (no_TRUNC_muts==0 & no_NTDam_muts==0 & no_GOF_muts==0 & CNVGain==0 & CNVLoss==0 & ExpT_NET==0), T, F))
    cat("Adding type annotation for training...", "\n")
    #df = df %>% mutate(type=ifelse(TP==T & vogel=="Vog.Oncogene", "O", ifelse(TP==T & vogel=="Vog.TS", "T",
    #                                                                          ifelse(TN==T, "N", NA))))

    ## Some of the cancer specific genes are also in the true negatives list (in NCG we only show a warning)
    #df = df %>% mutate(type=ifelse((TP==T & TN==F), "C", ifelse((TN==T & TP==F) | (TN==T & TP==T), "N", NA)))
    df = df %>% mutate(type=ifelse(TP==T, "C", NA))

    df = df %>% select(-TP)
    ## For cancer-specific training set
    df = df %>% select(-cancer_type, -primary_site, -cancer_site)
    df$type = as.factor(df$type)
    ## Exclude genes without any alterations from training & prediction set
    if(gains){
        df = df %>% subset(no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | (CNVGain==1 & Copy_number>=4) |
                               (CNVLoss==1 & Copy_number==0) | ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) & !(Entrez%in%false_positive_genes$entrez))
    }else{
        df = df %>% subset(no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 |
                               (CNVLoss==1 & Copy_number==0) | ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) & !(Entrez%in%false_positive_genes$entrez))
    }
    return(df)
}

markTrainCclass <- function(df, fp_dir="/mnt/lustre/users/k1469280/mourikisa/data/NCG_false_positives.txt", ...){

    ## To take trainig set based on the candidate genes for each cancer type
    ## Load geneInfo.Rdata and cancerGenes.Rdata from /mnt/lustre/users/k1469280/mourikisa/data in athena
    load("/mnt/lustre/users/k1469280/mourikisa/data/geneInfoNCG5_2.Rdata")
    load("/mnt/lustre/users/k1469280/mourikisa/data/cancerGenesNCG5_2.Rdata")

    ## Get information for cgc/cans and primary sites
    ## Fix gene info table from NCG
    geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
    ## Get a cancer gene with all the associated primary sites and cancer sites
    cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
        group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
                                       cancer_site=paste(unique(cancer_site), collapse=",")) %>%
        ungroup
    geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))

    ## Get false positive genes
    false_positive_genes <- read.table(fp_dir, header = T, sep = "\t")
    false_positive_genes <- false_positive_genes %>% select(entrez, symbol)

    ## Get primary site information for the df
    df = df %>% left_join(geneInfo, by=c("Entrez"="entrez"))

    ## Training set extraction
    cat("Marking genes for training...", "\n")
    #df = df %>% mutate(TP=ifelse((vogel=="Vog.Oncogene" | vogel=="Vog.TS") & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | CNVGain==1 | CNVLoss==1 | ExpT_NET==1), T, F))
    df = df %>% mutate(TP=ifelse((grepl("blood", primary_site)) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | (CNVGain==1 & Copy_number>=4) |
                                                                       (CNVLoss==1 & Copy_number==0) | ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0))), T, F))
    #df = df %>% mutate(TN1=ifelse(Entrez%in%false_positive_genes$entrez & , T, F))
    df = df %>% mutate(TN=ifelse(grepl("blood", primary_site) & TP==F, T, F))
    cat("Adding type annotation for training...", "\n")
    #df = df %>% mutate(type=ifelse(TP==T & vogel=="Vog.Oncogene", "O", ifelse(TP==T & vogel=="Vog.TS", "T",
    #                                                                          ifelse(TN==T, "N", NA))))

    ## Some of the cancer specific genes are also in the true negatives list (in NCG we only show a warning)
    df = df %>% mutate(type=ifelse((TP==T & TN==F), "C", ifelse(TN==T & TP==F, "N", NA)))

    df = df %>% select(-TP, -TN)
    ## For cancer-specific training set
    df = df %>% select(-cancer_type, -primary_site, -cancer_site)
    df$type = as.factor(df$type)
    ## Exclude genes without any alterations from training & prediction set
    df = df %>% subset((no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | (CNVGain==1 & Copy_number>=4) |
                           (CNVLoss==1 & Copy_number==0) | ((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0))) | !is.na(type))
    return(df)
}

## Scaling
getScaledTable <- function(df, excludeColumn=c( "Cancer_type", "Sample", "Entrez","type", "vogel"), type="cs",means_sds=NULL, ...){ #

    ## There is no is.nan for data frames!
    is.nan.data.frame <- function(x){do.call(cbind, lapply(x, is.nan))}

    cat("--> getScaledTable: All factors are excluded from scaling...\n")

    factor_cols <- names(sapply(df, function(x) class(x))[sapply(df, function(x) class(x))=="factor"])

    df_f <- df[,names(df)%in%c(excludeColumn,factor_cols)]

    df <-df[,!(names(df)%in%c(excludeColumn,factor_cols))]

    if(is.null(means_sds)){
        ## Get mean and sd for scaling new observations
        mean_sd = sapply(df, function(cl) list(means=mean(cl,na.rm=TRUE), sds=sd(cl,na.rm=TRUE)))

        cat("Scaling...", "\n")
        if (type=="cs"){
            df <- data.frame(lapply(df, function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))
            df <- cbind(df, df_f)
        }else if(type=="range"){
            df <- data.frame(lapply(df, function(x) rescale(x, to=c(-1,1))))
            df <- cbind(df, df_f)
        }
    }else{
        mean_sd = NULL
        cat("Scaling...", "\n")
        if (type=="cs"){
            df <- data.frame(lapply(seq_along(df), function(y,n,i) {
                (y[i]-unlist(means_sds["means",n[i]]))/unlist(means_sds["sds",n[i]])
            }, y=df, n=names(df)))
            df <- cbind(df, df_f)
        }
    }

    return(list(df_scaled=df, mean_sd=mean_sd))
}

write.config = function(x=NULL, res_dir, fname="config.txt"){
  if (!file.exists(res_dir)){
    dir.create(file.path(res_dir))
  }
  y = c('Cancer_type', 'Sample', 'Entrez', 'no_NSI_muts','no_NTDamFunction_muts', 'no_NTDamCons_muts', 'no_NTDamSC_muts','TPM', 'vogel', 'inCDD')
  if (length(x)>0) y = c(y,x)
  write.table(y, file=paste0(res_dir,'/',fname), sep="\t", row.names = F, col.names = F, quote=F)
}

## Clean features
cleanFeatures = function(df, config=NULL){
  # cat("Cleaning features for SVM model...\n")
  if(!is.null(config)){
    if(file.exists(config)){
      to_remove = read.delim(config, sep="\t", header=F)[,1]
      if(length(to_remove)>0){
         df = df %>% mutate(rowname=paste(Cancer_type, Sample, Entrez, sep="."))
        # df = df %>% mutate(rowname=paste(Sample, Entrez, sep="."))

        df = df[,which(!colnames(df)%in%to_remove)]

      }
    }
  }else{
    cat("--> cleanFeatures: it seems you don't have a cleaning config file...i'm going with defalut cleaning\n")
    df = df %>% mutate(rowname=paste(Cancer_type, Sample, Entrez, sep=".")) %>%
      select(-Cancer_type, -Sample, -Entrez, -no_NSI_muts, -no_NTDamFunction_muts, -no_NTDamCons_muts, -no_NTDamSC_muts,
             -TPM, -vogel, -inCDD)
  }

    rnames= df$rowname
    df = df %>% select(-rowname)
    df = data.frame(df)
    row.names(df) = rnames
    return(df)
}

## Check training distributions
checkTraining = function(df, saveDir){
    factor_cols <- names(sapply(df, function(x) class(x))[sapply(df, function(x) class(x))=="factor"])
    pdf(file=paste(saveDir,"/training_feature_dists.pdf",sep=""),height = unit(12, "inches"), width = unit(20, "inches"), useDingbats = F)
    total_stats = NULL
    for (c in colnames(df)[!(colnames(df)%in%c(factor_cols, "type"))]){
        print(c)
        a = df[,c(c,"type")]
        stats = data.frame(
        mi = summary(a[,1])[[1]],
        fQ = summary(a[,1])[[2]],
        md = summary(a[,1])[[3]],
        mn = summary(a[,1])[[4]],
        tQ = summary(a[,1])[[5]],
        mx = summary(a[,1])[[6]],
        row.names = c
        )
        total_stats = rbind(total_stats, stats)
        p = ggplot(a, aes_string(x=colnames(a)[2], y=colnames(a)[1],fill = colnames(a)[2])) +
                geom_boxplot() +
                theme(
                    axis.title.x = element_blank(),
                    axis.text.x = element_text(size = 12, color="black"),
                    axis.text.y = element_text(size = 12)
                ) +
                ylab(paste(colnames(a)[1],"(n)",sep="")) +
                ggtitle(colnames(a)[1])
        print(p)
    }
    write.table(total_stats, file=paste(saveDir,"/training_feature_dists.tsv",sep=""), quote = F, row.names = T, sep = "\t")
    for (c in colnames(df)[colnames(df)%in%c(factor_cols, "type")]){
        a = df[,c(c,"type")]
        a = table(a) %>% data.frame()
        p = ggplot(a, aes_string(x=colnames(a)[2], y="Freq",fill = colnames(a)[1])) +
            geom_bar(stat = "identity", position = "dodge") +
            theme(
                axis.title.x = element_blank(),
                axis.text.x = element_text(size = 12, color="black"),
                axis.text.y = element_text(size = 12)
            ) +
            ylab(paste(colnames(a)[1],"(n)",sep="")) +
            ggtitle(colnames(a)[1])
        print(p)
    }
    dev.off()
}

## Get a table with the commands to run in athena
## Starting range nu: seq(0.05, 0.9, 0.05)
## Starting range gamma: 2^seq(-7,4)
runCommands <- function(ct, kern, mynu, mygamma, mydegree, save_dir, fname="jobs.txt", ncore=1, res_dir = NULL, source_dir=NULL){
    if(is.null(res_dir)) res_dir = paste0(save_dir,'/',ct)
    jobName = paste0(kern, "_", ct)
    sub <- paste0("qsub -N ", jobName," -pe smp ",ncore," -l h_vmem=12G /mnt/lustre/users/k1469280/mourikisa/data/OAC/submitNoveltyDetection.sh")
    params <- expand.grid(mynu, mygamma, mydegree)
    params$p <- paste(params$Var1, params$Var2, params$Var3, sep=" ")
    df <- NULL
    for (p in unique(params$p)){
      cmd <- paste(sub, ct, kern, p, res_dir, source_dir, sep=" ")
      df <- rbind(df, cmd)
    }
    write.table(df, sep="\t", quote=F, row.names=F, col.names=F, file=paste0(save_dir,"/",fname))
}

## Implementation of Nelder-Mead optimization to find parameters giving best cross-validation
nmOptimizeLinear = function(training, startingPoints=3){
  x = training %>% subset(type=="C") %>% dplyr::select(-type) %>% data.matrix()
  y = rep(1, nrow(x))

  ## one-class SVM
  ocSVM <- function(x, y, xtest, ytest, mynu) {
    svm.model <- svm(x,y,kernel="linear", type="one-classification", scale=FALSE, nu=mynu)
    predict(svm.model, xtest, decision.values=T)
  }
  cv <- cv.setup(x, y, score=accuracy, num_folds = 3, num_iter = 3)
  ## Choose 100 random starting points of nu
  starting_nu = runif(startingPoints, min=0, max=0.9)
  optima = list()
  pb <- txtProgressBar(min = 0, max = startingPoints, style = 3)
  for (i in 1:startingPoints){
    Sys.sleep(0.1)
    res <- cv.nelder_mead(cv, ocSVM, mynu = starting_nu[i], num_evals = 10, maximize=TRUE)
    optima[paste(i,"_optimum", sep="")] = res$optimum
    optima[paste(i,"_starting_nu", sep="")] = starting_nu[i]
    optima[paste(i,"_solution_nu", sep="")] = res$solution
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(optima)
}

nmOptimizeRadial = function(training, startingPoints=3){
  x = training %>% subset(type=="C") %>% dplyr::select(-type) %>% data.matrix()
  y = rep(1, nrow(x))

  ## one-class SVM
  ocSVM <- function(x, y, xtest, ytest, mynu, mygamma) {
    svm.model <- svm(x,y,kernel="radial", type="one-classification", scale=FALSE, nu=mynu, gamma=mygamma)
    predict(svm.model, xtest, decision.values=T)
  }
  cv <- cv.setup(x, y, score=accuracy, num_folds = 3, num_iter = 3)
  ## Choose random starting points of nu
  starting_nu = runif(startingPoints, min=0, max=0.9)
  starting_gamma = 2^seq(-6,3)
  stps = expand.grid(starting_nu, starting_gamma)
  colnames(stps) = c("nu", "gamma")
  optima = list()
  pb <- txtProgressBar(min = 0, max = nrow(stps), style = 3)
  for (i in 1:nrow(stps)){
    Sys.sleep(0.1)
    res <- cv.nelder_mead(cv, ocSVM, mynu = stps$nu[i], mygamma=stps$gamma[i], num_evals = 10, maximize=TRUE)
    optima[paste(i,"_optimum", sep="")] = res$optimum
    optima[paste(i,"_starting_nu", sep="")] = stps$nu[i]
    optima[paste(i,"_starting_gamma", sep="")] = stps$gamma[i]
    optima[paste(i,"_solution_nu", sep="")] = res$solution$mynu
    optima[paste(i,"_solution_gamma", sep="")] = res$solution$mygamma
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(optima)
}

## After the new implementation in Novelty detection for each model I create the
## stats, so this function simply concatenates everything
gatherMetrics = function(res_dir){
  total_cv_stats = NULL
  modelDirs = list.dirs(res_dir, recursive=F)
  for (md in modelDirs){
    if (length(grep("cv_stats.tsv", list.files(md, recursive=F)))!=0){
      df = read.table(paste(md, "/cv_stats.tsv", sep=""), header = T, sep="\t")
      if (!("degree" %in% colnames(df))){
          df = df %>% mutate(degree=0) %>%
              select(analysis,kernel,nu,gamma,degree,iteration,type,trainSize,class,set,value)
      }
      total_cv_stats = rbind(total_cv_stats, df)
    }else{
      cat(paste(md, " has no cv_stats.tsv", sep=""), "\n")
    }
  }

  return(total_cv_stats)
}


## Given a directory gather results of cross-validation
getCVStats <- function(res_dir, svmmode="C-classification", classes=2){

    cat("Creating the folder list...", "\n")
    cv_dirs <- list.dirs(res_dir, recursive = F)
    all_fns <- NULL
    ## Make the list of the directories to be gathered
    for (cv_dir in cv_dirs){
        setwd(cv_dir)
        iter_dirs <- list.dirs(cv_dir, recursive = F)
        for (iter_dir in iter_dirs){
            cat(iter_dir, "\n")
            setwd(iter_dir)
            ## Check if all 4 files are there (for cases where the cross-validation still running)
            fns <- c("trainset.Rdata", "testset.Rdata", "svmModel.Rdata", "prediction.Rdata")
            if (sum(file.exists(fns), na.rm=TRUE)==4){
                ## Load the trainset, testset, model & prediciton
                load("testset.Rdata")
                load("prediction.Rdata")

                ## Check point for prediction hault
                if(length(prediction)!=nrow(testset)){
                    cat(paste(iter_dir," skipped! Not equal number of predictions and test set observations!", sep=""), "\n")
                    next
                }else{
                    all_fns <- c(all_fns, iter_dir)
                }

            }else{
                cat(paste(iter_dir, " skipped! Not all the files are present or the structur eof directory is wrong!", sep=""), "\n")
                next
            }
        }
    }

    cat("initiating data mining from cross-validation results...", "\n")
    ## Run it in parallel
    total_cv_stats <- foreach(i=all_fns, .combine=rbind) %dopar% {
        cv_analysis <- tail(unlist(strsplit(i, split = "/")), n=2)[1]
        setwd(i)
        cat(i, "\n")
        iter <- as.numeric(tail(unlist(strsplit(i, split = "\\.")), n=1))
        cv_stats <- NULL
        load("trainset.Rdata")
        load("testset.Rdata")
        load("svmModel.Rdata")
        load("prediction.Rdata")

        noTrain = trainset %>% nrow
        ## Concatenate the predictions to the true values
        testset$prediction = prediction
        testset = testset %>% select(type, prediction)
        ## For C-classification/one-classification
        if (svmmode=="one-classification"){
            testset = testset %>% mutate(p=ifelse(prediction==TRUE, "C", "N"))
        }else if (svmmode=="C-classification"){
            testset$p = testset$prediction
        }

        ## Performance metrics

        ## ---------------------- Test set -------------------------------------
        ## Confusion matrix
        if (classes==3){
            cM <- confusionMatrix(testset$p,testset$type, positive = c("T", "O"))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="accuracy", trainSize=noTrain,
                                                   class="overall", set="test",
                                                   value=cM$overall[["Accuracy"]]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Sensitivity", trainSize=noTrain,
                                                   class="N", set="test",
                                                   value=cM$byClass["Class: N","Sensitivity"]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Sensitivity", trainSize=noTrain,
                                                   class="T", set="test",
                                                   value=cM$byClass["Class: T","Sensitivity"]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Sensitivity", trainSize=noTrain,
                                                   class="O", set="test",
                                                   value=cM$byClass["Class: O","Sensitivity"]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Specificity", trainSize=noTrain,
                                                   class="N", set="test",
                                                   value=cM$byClass["Class: N","Specificity"]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Specificity", trainSize=noTrain,
                                                   class="T", set="test",
                                                   value=cM$byClass["Class: T","Specificity"]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Specificity", trainSize=noTrain,
                                                   class="O", set="test",
                                                   value=cM$byClass["Class: O","Specificity"]))
        }else if (classes==2){
            cM <- confusionMatrix(testset$p,testset$type, positive = c("C"))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="accuracy", trainSize=noTrain,
                                                   class="overall", set = "test",
                                                   value=cM$overall[["Accuracy"]]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Sensitivity", trainSize=noTrain,
                                                   class="overall", set = "test",
                                                   value=cM$byClass["Sensitivity"]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Specificity", trainSize=noTrain,
                                                   class="overall", set = "test",
                                                   value=cM$byClass["Specificity"]))
        }

        ## Area under the curve
        #auc <- multiclass.roc(as.numeric(pred$type), as.numeric(pred$prediction), levels=c(1,2,3))$auc[1]


        ## ---------------------- Train set -------------------------------------
        ## Confusion matrix
        trainset$fitted = svm.model$fitted
        trainset = trainset %>% select(type, fitted)
        ## For C-classification/one-classification
        if (svmmode=="one-classification"){
            trainset = trainset %>% mutate(f=ifelse(fitted==TRUE, "C", "N"))
        }else if (svmmode=="C-classification"){
            trainset$f = trainset$fitted
        }

        ## Confusion matrix
        if (classes==3){
            cM <- confusionMatrix(trainset$f,trainset$type, positive = c("T", "O"))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="accuracy", trainSize=noTrain,
                                                   class="overall", set="train",
                                                   value=cM$overall[["Accuracy"]]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Sensitivity", trainSize=noTrain,
                                                   class="N", set="train",
                                                   value=cM$byClass["Class: N","Sensitivity"]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Sensitivity", trainSize=noTrain,
                                                   class="T", set="train",
                                                   value=cM$byClass["Class: T","Sensitivity"]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Sensitivity", trainSize=noTrain,
                                                   class="O", set="train",
                                                   value=cM$byClass["Class: O","Sensitivity"]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Specificity", trainSize=noTrain,
                                                   class="N", set="train",
                                                   value=cM$byClass["Class: N","Specificity"]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Specificity", trainSize=noTrain,
                                                   class="T", set="train",
                                                   value=cM$byClass["Class: T","Specificity"]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Specificity", trainSize=noTrain,
                                                   class="O", set="train",
                                                   value=cM$byClass["Class: O","Specificity"]))
        }else if (classes==2){
            cM <- confusionMatrix(trainset$f,trainset$type, positive = c("C"))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="accuracy", trainSize=noTrain,
                                                   class="overall", set = "train",
                                                   value=cM$overall[["Accuracy"]]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Sensitivity", trainSize=noTrain,
                                                   class="overall", set = "train",
                                                   value=cM$byClass["Sensitivity"]))
            cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, iteration=iter,
                                                   type="Specificity", trainSize=noTrain,
                                                   class="overall", set = "train",
                                                   value=cM$byClass["Specificity"]))
        }
        return(cv_stats)
    }
}

plotCVStats = function(total_cv_stats, saveDir){

    #iterations per analsyis
    an2iter <- total_cv_stats %>% select(analysis, iteration) %>% unique %>% group_by(analysis) %>% summarise(n=n()) %>% ungroup %>% rename(iterations=n)

    cv_stats_summary <- total_cv_stats %>% subset(!is.na(value)) %>% group_by(analysis, kernel, nu, gamma, type, set) %>% summarize(min=min(value),
                                                                                q1=quantile(value)[2],
                                                                                median=median(value),
                                                                                mean=mean(value),
                                                                                q3=quantile(value)[4],
                                                                                max=max(value)) %>% ungroup

    cv_stats_summary <- cv_stats_summary %>% left_join(an2iter, by="analysis")

    ## Plot overall sensitivity and specificity across different parameters
    #overall <- cv_stats_summary %>% subset(Cancer_type=="overall")
    overall <- cv_stats_summary %>% subset(type=="Sensitivity" | type=="Specificity")
    overall <- overall %>% select(analysis, kernel, nu, gamma, type, set, median) %>% spread(type, median)

    ## Plot test stats
    for (k in unique(overall$kernel)){

        df = overall %>% subset(kernel==k & set=="test")
        an2sy = df %>% select(analysis) %>% unique
        an2sy$symbol <- 1:nrow(an2sy)
        df = df %>% left_join(an2sy, by=c("analysis"))
        pdf(file=paste(saveDir,"/",k,"_CV_Test_ROCspace.pdf",sep=""),height = unit(12, "inches"), width = unit(20, "inches"), useDingbats = F)
        p <- ggplot(df,aes(x=1-Specificity,y=Sensitivity,group=analysis,shape=analysis)) +
        geom_text(aes(label=symbol))+
        #geom_point(size=5) +
        #geom_jitter(position = position_jitter(), size=5) +
        scale_x_continuous(limit=c(0,1)) +
        scale_y_continuous(limit=c(0,1)) +
        geom_abline(intercept = 0, slope =1, lty="dashed", color="red") +
        ggtitle(paste("ROC space (",k,")", sep = "")) +
        xlab("FPR (1-Specificity)") +
        ylab("TPR (Sensitivity)") +
        theme(
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14),
            axis.text.x=element_text(colour="black", size=12),
            axis.text.y=element_text(colour = "black", size=12),
            axis.line = element_line(colour = "black"),
            plot.title=element_text(size = 16),
            legend.text=element_text(size=16),
            legend.key.size = unit(3, "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()
        ) +
        scale_shape_manual(values = df$symbol)
        print(p)
        dev.off()
        write.table(df, file=paste(saveDir,"/",k,"_CV_Test_ROCspace.tsv",sep=""), quote = F, row.names=F, sep="\t")
    }

    ## Plot train stats
    for (k in unique(overall$kernel)){
      df = overall %>% subset(kernel==k & set=="train")
      an2sy = df %>% select(analysis) %>% unique
      an2sy$symbol <- 1:nrow(an2sy)
      df = df %>% left_join(an2sy, by=c("analysis"))
      pdf(file=paste(saveDir,"/",k,"_CV_Train_Sensitivity.pdf",sep=""),height = unit(12, "inches"), width = unit(20, "inches"), useDingbats = F)
      p = ggplot(df, aes(x=factor(analysis), y=Sensitivity)) +
        geom_text(aes(label=symbol)) +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=90, size=12),
          axis.text.y = element_text(size=12)
        ) + ggtitle(paste("Training sensitivity (",k,")", sep = ""))
      print(p)
      dev.off()
      write.table(df, file=paste(saveDir,"/",k,"_CV_Train_Sensitivity.tsv",sep=""), quote = F, row.names=F, sep="\t")
    }
}

################################################
# Platt scaling given a prediction
################################################
## For the moment the sigmoid is trained on the training values of the iteration selected
## and assess the predictions made in the validation set using this model
plattScaling = function(model, predictions){
  ## Impementation was taken from pseudocode in:
  ## Probabilistic Outputs for Support Vector Machines and Comparisons to Regularized Likelihood Methods, Platt, 1999
  ## Construct the training set for Platt's sigmoid

  decision_values = model$decision.values %>% as.vector()
  labels = model$fitted %>% as.vector()
  df = data.frame(dv=decision_values, lbs = labels)

  ## the prediction at the moment is a predict object coming from svm predict
  dv = attr(predictions, "decision.values") %>% as.vector()
  preds = data.frame(dv = dv)
  preds$label = as.vector(predictions)
  rownames(preds) = names(predictions)

  ## Input parameters
  out = df$dv
  target = df$lbs
  prior1 = df%>%subset(lbs==TRUE)%>%nrow
  prior0 = df%>%subset(lbs==FALSE)%>%nrow
  ## Output
  ## A, B = parameters of sigmoid

  A = 0
  B = log((prior0+1)/(prior1+1))
  hiTarget = (prior1+1)/(prior1+2)
  loTarget = 1/(prior0+2)
  lambda = 1e-3
  olderr = 1e300


  pp = rep((prior1+1)/(prior0+prior1+2), nrow(df))
  count = 0

  for (it in 1:100){
    a=0
    b=0
    c=0
    d=0
    e=0

    ## First compute the Hessian & gradient of error function
    ## with respect to A & B
    for (i in 1:length(pp)){
      if (target[i]==TRUE){
        t = hiTarget
      }else{
        t = loTarget
      }

      d1 = pp[i]-t
      d2 = pp[i]*(1-pp[i])
      a = a + out[i]*out[i]*d2
      b = b + d2
      c = c + out[i]*d2
      d = d + out[i]*d1
      e = e + d1
    }
    ## If gradient is really tiny, then stop
    if ((abs(d) < 1e-9) & (abs(e) < 1e-9)){break}

    oldA = A
    oldB = B
    err = 0

    ## Loop until goodness of fit increases
    while(TRUE){
      det = (a+lambda)*(b+lambda)-c*c
      if (det==0){ ## If determinant if Hessian is 0
        ## Increase stabilizer
        lambda = lambda*10
        next
      }
      A = oldA + ((b+lambda)*d-c*e)/det
      B = oldB + ((a+lambda)*e-c*d)/det

      ## Now compute the goodness of fit
      err = 0
      for (i in 1:length(pp)){
        p = 1/(1+exp(out[i]*A+B))
        pp[i] = p
        ## At this step, make sure log(0) returns -200
        if (p<=1.383897e-87){
          err = err - t*(-200)+(1-t)*log(1-p)
        }else if (p==1){
          err = err - t*log(p)+(1-t)*(-200)
        }else{
          err = err - t*log(p)+(1-t)*log(1-p)
        }
        if(err==-Inf) browser()
      }
      if (err < olderr*(1+1e-7)){
        lambda = lambda*0.1
        break
      }
      ## Error did not decrease: increase stabilizer by factor of 10
      ## and try again
      lambda = lambda* 10
      if (lambda >= 1e6){ ## something is broken. Give up
        break
      }
    }
    diff = err-olderr
    scale = 0.5*(err+olderr+1)
    if (diff > -1e-3*scale & diff < 1e-7*scale){
      count = count + 1
    }else{
      count = 0
    }
    olderr = err
    if (count==3){
      break
    }


  }

  applyPlatt = function(x){
    p = 1/(1+exp(x*A+B))
  }

  preds$prob = apply(preds, 1, function(x) applyPlatt(x[1]))
  return(preds)
}

##########
# Binning
##########
binning = function(model, predictions){

  ## Make a dataframe with decision values and labes from the training object
  df = data.frame(label=model$fitted, dv=as.vector(model$decision.values))
  ## Number of bins
  bins = 10
  dv = sort(df$dv, decreasing=FALSE)
  ## Find the bounds of each bin in the training vector
  ## Carefull this function returns the elements of the vector encoded as ranges
  dv_encoded = cut(dv, breaks=bins)
  bin_bounds = unique(cut(dv, breaks=bins))
  probs = NULL
  for (i in 1:length(bin_bounds)){
    bin = bin_bounds[i]
    lowerLimit = as.numeric(gsub("\\(", "", unlist(strsplit(as.character(bin), ","))[1]))
    upperLimit = as.numeric(gsub("\\]", "", unlist(strsplit(as.character(bin), ","))[2]))

    labelMatrix = df %>% subset(dv >= lowerLimit & dv < upperLimit) %>% count(label)
    if (length(unique(labelMatrix$label))==2){
      prob = labelMatrix$n[labelMatrix$label==TRUE]/sum(labelMatrix$n)
    }else if (length(unique(labelMatrix$label))==1){
      if (unique(labelMatrix$label)==FALSE){
        prob = 0
      }else if (unique(labelMatrix$label)==TRUE){
        prob = 1
      }
    }
    probs = c(probs, prob)
  }
  probsDF = data.frame(bin=bin_bounds, prob=probs)
  probsDF = probsDF %>% separate(bin, into = c("lowLimit", "upLimit"), sep=",", remove = FALSE) %>%
    mutate(lowLimit=as.numeric(gsub("\\(", "", lowLimit)),
           upLimit=as.numeric(gsub("\\]", "", upLimit)))

  ## Get the prediction decision values
  preds = attr(predictions, "decision.values")[,1] %>% data.frame()
  preds$label = as.vector(predictions)
  colnames(preds) = c("dv", "label")

  getBinProb = function(x){
    prob = subset(probsDF, x >= lowLimit & x < upLimit)$prob[1]
    ## if prob = NA it means that the value is outside of the bins
    ## Assign it to the closest bin
    ## Get the maximum and the minimum
    maximum = max(probsDF$upLimit)
    minimum = min(probsDF$lowLimit)
    if (is.na(prob)){
      if (x>maximum){
        prob = 1
      }else if (x<minimum){
        prob = 0
      }
    }
    return(prob)
  }
  preds$prob = sapply(preds$dv, function(x) getBinProb(x))
  return(preds)
}

## After inspection of the ROC space, we decide the best model based on the CV
## I give the parameters to the function, train the model on the whole training
## training set and predicts on the validation set
## the function also implements probability estimation using Platt's algorithm
## and binning
## For Platt training set at the moment I use the decision values of the training
## You need to open the tunnel to NCG5 brfore running the function
bestModelPredict  =function(kern, mynu, mygamma, trainingPath, validationPath, validation_nsPath){

  dev_mode()
  install("/Users/fc-8s-imac/ncg-data-update/ncglib/")
  library(ncglib)
  geneInfo=get_geneinfo()
  cancerGenes = get_cancer_genes()
  dev_mode()

  ## Load training and validation sets
  load(trainingPath)
  load(validationPath)
  load(validation_nsPath)
  ## Fix gene info table from NCG
  geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
  ## Get a cancer gene with all the associated primary sites and cancer sites
  cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
    group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
                                   cancer_site=paste(unique(cancer_site), collapse=",")) %>%
    ungroup
  geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))

  if (kern=="linear"){
    mygamma = NULL
  }

  ## Training also contains "N" class
  training = training%>%subset(type=="C")

  ## Train on the whole set
  if (kern=="linear"){
    svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, nu=mynu)
  }else if (kern=="radial"){
    svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, gamma = mygamma, nu=mynu)
  }

  ## make prediction
  predictions = predict(svm.model, validation, decision.values = TRUE)
  ## add probability from Platt and Binning
  predictions_platt = plattScaling(svm.model, predictions)
  colnames(predictions_platt) = c("dv", "label", "prob_platt")
  #predictions_binning = binning(svm.model, predictions)
  #colnames(predictions_binning) = c("dv", "label", "prob_binning")
  preds = predictions_platt %>% tibble::rownames_to_column() #%>%
  #left_join(predictions_binning%>%add_rownames(), by=c("rowname","dv", "label"))

  ## Take information from gene info
  preds = preds %>%
    separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>%
    mutate(entrez=as.numeric(entrez)) %>%
    left_join(geneInfo%>%rename(gene_type=cancer_type), by=c("entrez"))

  ## Join also with systems level properties from validaton set
  validation_ns = validation_ns %>% tibble::rownames_to_column() %>%
      separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>%
      mutate(entrez=as.numeric(entrez)) %>%
      select(cancer_type, sample, entrez, no_TRUNC_muts, no_NTDam_muts,
             no_GOF_muts, Copy_number, CNVLoss, CNVGain)
  preds = preds %>% left_join(validation_ns, by=c("cancer_type", "sample", "entrez"))

  return(preds)
}

## Exactly the same as bestModelPredict but gives predictions for the whole grid
getGridPred = function(kern, mynu, mygamma, mydegree, dataPath, res_dir, tissue, gains=T, gof=T, install_NCG=TRUE, NCG_path="/Users/fc-8s-imac/ncg-data-update/ncglib/",fname_gene_prop = "Rdata/geneProperties.Rdata"){
    if (!file.exists(res_dir)){
        dir.create(file.path(res_dir))
    }

    ## What features should be excluded
    config = read.table(paste(dataPath, "/config.txt", sep=""), header = F)
    config = config[10:nrow(config), 1]
    s = c("cancer_type", "sample", "entrez", "no_ALL_muts","no_TRUNC_muts", "no_NTDam_muts",
          "no_GOF_muts", "Copy_number", "CNVLoss", "CNVGain")
    s = s[!(s%in%config)]


    dev_mode()
    if(install_NCG) install(NCG_path)
    library(ncglib)
    geneInfo=get_geneinfo(version="NCG5")
    cancerGenes = get_cancer_genes(version="NCG5")
    dev_mode()

    ## Load training and validation sets
    load(paste(dataPath, "/training_set.Rdata", sep=""))
    load(paste(dataPath, "/validation_set.Rdata", sep=""))
    load(paste(dataPath, "/validation_set_noScale.Rdata", sep=""))


    validation_ns = validation_ns %>% tibble::rownames_to_column() %>%
        separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>%
        mutate(entrez=as.numeric(entrez)) %>%
        select(one_of(s))


    ## Fix gene info table from NCG
    geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
    ## Get a cancer gene with all the associated primary sites and cancer sites
    cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
        group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
                                       cancer_site=paste(unique(cancer_site), collapse=",")) %>%
        ungroup
    geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))

    if (kern=="linear"){
        mygamma = 0
    }

    ## Training also contains "N" class? (just in case to be sure between the trials)
    training = training%>%subset(type=="C")

    ## Train for all the grid
    params <- expand.grid(mynu, mygamma, mydegree)
    colnames(params) = c("mynu", "mygamma", "mydegree")

    grid_predictions_table = NULL
    grid_predictions = list()
    novel_drivers_genes_platt_all = NULL
    novel_drivers_genes_platt_inNCG5_all = NULL
    novel_drivers_genes_platt_rst_all = NULL
    new_samples_covered_all = NULL

    for (i in 1:nrow(params)){
        ## Train on the whole set
        cat(i, "\n")
        if (kern=="linear"){
            svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, nu=params$mynu[i])
        }else if (kern=="radial"){
            svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, gamma=params$mygamma[i], nu=params$mynu[i])
        }else if (kern=="polynomial"){
            svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, gamma = params$mygamma[i], nu=params$mynu[i], degree = params$mydegree[i])
        }else if (kern=="sigmoid"){
            svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, gamma=params$mygamma[i], nu=params$mynu[i])
        }

        ## Asses training sensitivity
        trainset = training %>% select(type)
        trainset$fitted = svm.model$fitted
        trainset = trainset %>% select(type, fitted)
        trainset = trainset %>% mutate(f=ifelse(fitted==TRUE, "C", "N"))
        trainset$type = factor(as.character(trainset$type), levels = c("C", "N"))
        trainset$f = factor(as.character(trainset$f), levels = c("C", "N"))
        ## Confusion matrix
        cM <- confusionMatrix(trainset$f,trainset$type, positive = c("C"))
        training_sensitivity = cM$byClass["Sensitivity"]

        ## make prediction
        predictions = predict(svm.model, validation, decision.values = TRUE)
        ## add probability from Platt and Binning
        predictions_platt = plattScaling(svm.model, predictions)
        colnames(predictions_platt) = c("dv", "label", "prob_platt")
        #predictions_binning = binning(svm.model, predictions)
        #colnames(predictions_binning) = c("dv", "label", "prob_binning")
        preds = predictions_platt %>% tibble::rownames_to_column() #%>%
        #left_join(predictions_binning%>%add_rownames(), by=c("rowname","dv", "label"))

        ## Take information from gene info
        preds = preds %>%
            separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>%
            mutate(entrez=as.numeric(entrez)) %>%
            left_join(geneInfo%>%rename(gene_type=cancer_type), by=c("entrez"))

        ## Join also with systems level properties from validaton set
        preds = preds %>% left_join(validation_ns, by=c("cancer_type", "sample", "entrez"))

        ## Save preds
        if (kern=="linear"){
            analysis = paste(kern, ".", params$mynu[i], sep="")
            save(preds, file=paste(res_dir,"/predictions_",analysis, ".Rdata", sep=""))
        }else if (kern=="radial"){
            analysis = paste(kern, ".", params$mynu[i], ".", params$mygamma[i], sep="")
            save(preds, file=paste(res_dir, "/predictions_",analysis, ".Rdata", sep=""))
        }else if (kern=="polynomial"){
            analysis = paste(kern, ".", params$mynu[i], ".", params$mygamma[i], ".", params$mydegree[i], sep="")
            save(preds, file=paste(res_dir, "/predictions_",analysis, ".Rdata", sep=""))
        }else if (kern=="sigmoid"){
            analysis = paste(kern, ".", params$mynu[i], ".", params$mygamma[i], sep="")
            save(preds, file=paste(res_dir, "/predictions_",analysis, ".Rdata", sep=""))
        }


        ## Save the predicted entrez to make the venn diagram
        grid_predictions[[analysis]] = preds%>%subset(label==TRUE & prob_platt >0.95)%>%select(entrez)%>%unique%>%.$entrez

        training_samples = training%>%tibble::rownames_to_column()%>%separate(rowname, into=c("Cancer_type", "sample", "entrez"), sep="\\.")%>%select(sample)%>%unique%>%.$sample
        prediction_samples = preds%>%subset(label==TRUE & prob_platt >0.95)%>%select(sample)%>%unique%>%.$sample

        ## Cutoffs
        cutoffs = c(0.95, 0.85, 0.75)
        new_samples = preds%>%subset(label==TRUE & prob_platt >0.95)%>%select(sample)%>%subset(!(sample%in%training_samples))%>%unique%>%.$sample
        new_samples_no_0.95 = length(new_samples)
        new_samples = preds%>%subset(label==TRUE & prob_platt >0.85)%>%select(sample)%>%subset(!(sample%in%training_samples))%>%unique%>%.$sample
        new_samples_no_0.85 = length(new_samples)
        new_samples = preds%>%subset(label==TRUE & prob_platt >0.75)%>%select(sample)%>%subset(!(sample%in%training_samples))%>%unique%>%.$sample
        new_samples_no_0.75 = length(new_samples)

        ## Store in the vectors to get the unique numbers
#         novel_drivers_genes_platt_all = c(novel_drivers_genes_platt_all, preds%>%subset(label==TRUE & prob_platt >0.95)%>%select(entrez)%>%unique%>%.$entrez)
#         novel_drivers_genes_platt_inNCG5_all = c(novel_drivers_genes_platt_inNCG5_all, preds%>%subset(label==TRUE & prob_platt >0.95 & (gene_type=="cgc" | gene_type=="can"))%>%
#                                                      select(entrez)%>%unique%>%.$entrez)
#         novel_drivers_genes_platt_rst_all = c(novel_drivers_genes_platt_rst_all, preds%>%subset(label==TRUE & prob_platt >0.95 & gene_type=="rst")%>%select(entrez)%>%unique%>%.$entrez)
#         new_samples_covered_all = c(new_samples_covered_all, new_samples)



        df = data.frame(analysis=analysis, nu=params$mynu[i],
                        gamma=params$mygamma[i], degree=params$mydegree[i],
                        training_set_size = nrow(training),
                        training_genes = training%>%tibble::rownames_to_column()%>%separate(rowname, into=c("Cancer_type", "sample", "entrez"), sep="\\.")%>%select(entrez)%>%unique%>%nrow,
                        training_samples = training%>%tibble::rownames_to_column()%>%separate(rowname, into=c("Cancer_type", "sample", "entrez"), sep="\\.")%>%select(sample)%>%unique%>%nrow,
                        training_sensitivity = training_sensitivity,
                        prediction_set_size=nrow(preds),
                        prediction_genes=preds%>%select(entrez)%>%unique%>%nrow,
                        prediction_samples=preds%>%select(sample)%>%unique%>%nrow,
                        novel_drivers_entries=preds%>%subset(label==TRUE)%>%nrow,
                        novel_drivers_genes=preds%>%subset(label==TRUE)%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_genes_inNCG5 = preds%>%subset(label==TRUE & (gene_type=="cgc" | gene_type=="can"))%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_genes_rst = preds%>%subset(label==TRUE & gene_type=="rst")%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_entries_platt_0.95=preds%>%subset(label==TRUE & prob_platt >0.95)%>%nrow,
                        novel_drivers_genes_platt_0.95=preds%>%subset(label==TRUE & prob_platt >0.95)%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_in_samples_mutated_0.95 = preds%>%subset(label==TRUE & prob_platt >0.95)%>%select(sample)%>%unique%>%nrow,
                        novel_drivers_genes_platt_inNCG5_0.95 =preds%>%subset(label==TRUE & prob_platt >0.95 & (gene_type=="cgc" | gene_type=="can"))%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_genes_platt_rst_0.95=preds%>%subset(label==TRUE & prob_platt >0.95 & gene_type=="rst")%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_genes_platt_rst_in_samples_mutated_0.95=preds%>%subset(label==TRUE & prob_platt >0.95 & gene_type=="rst")%>%select(sample)%>%unique%>%nrow,
                        new_samples_covered_0.95 = new_samples_no_0.95,
                        novel_drivers_entries_platt_0.85=preds%>%subset(label==TRUE & prob_platt >0.85)%>%nrow,
                        novel_drivers_genes_platt_0.85=preds%>%subset(label==TRUE & prob_platt >0.85)%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_in_samples_mutated_0.85 = preds%>%subset(label==TRUE & prob_platt >0.85)%>%select(sample)%>%unique%>%nrow,
                        novel_drivers_genes_platt_inNCG5_0.85 =preds%>%subset(label==TRUE & prob_platt >0.85 & (gene_type=="cgc" | gene_type=="can"))%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_genes_platt_rst_0.85=preds%>%subset(label==TRUE & prob_platt >0.85 & gene_type=="rst")%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_genes_platt_rst_in_samples_mutated_0.85=preds%>%subset(label==TRUE & prob_platt >0.85 & gene_type=="rst")%>%select(sample)%>%unique%>%nrow,
                        new_samples_covered_0.85 = new_samples_no_0.85,
                        novel_drivers_entries_platt_0.75=preds%>%subset(label==TRUE & prob_platt >0.75)%>%nrow,
                        novel_drivers_genes_platt_0.75=preds%>%subset(label==TRUE & prob_platt >0.75)%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_in_samples_mutated_0.75 = preds%>%subset(label==TRUE & prob_platt >0.75)%>%select(sample)%>%unique%>%nrow,
                        novel_drivers_genes_platt_inNCG5_0.75 =preds%>%subset(label==TRUE & prob_platt >0.75 & (gene_type=="cgc" | gene_type=="can"))%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_genes_platt_rst_0.75=preds%>%subset(label==TRUE & prob_platt >0.75 & gene_type=="rst")%>%select(entrez)%>%unique%>%nrow,
                        novel_drivers_genes_platt_rst_in_samples_mutated_0.75=preds%>%subset(label==TRUE & prob_platt >0.75 & gene_type=="rst")%>%select(sample)%>%unique%>%nrow,
                        new_samples_covered_0.75 = new_samples_no_0.75
                        )
        grid_predictions_table = rbind(grid_predictions_table, df)

        ## make the plots and save them for each iteration
        for (cutoff in cutoffs){
            cat(cutoff, "\n")
            positive_genes_platt_rst=preds%>%subset(label==TRUE & prob_platt >cutoff & gene_type=="rst")%>%select(entrez)%>%unique%>%nrow
            if (positive_genes_platt_rst==0){
                next
            }else{
                p = plotSysproperties(preds, cutoff, fname_gene_prop, tissue=tissue, install_NCG=install_NCG )
                mylegend<-g_legend(p[[1]])
                for(i in 1:length(p)) p[[i]] = p[[i]] + theme(legend.position="none")
                pdf(file=paste(res_dir, "/predictions_sysproperties_",analysis,"_",cutoff, ".pdf" , sep=""), h= 20, w=24)
                grid.arrange(p[[1]],p[[3]],p[[6]],p[[7]],p[[8]],
                             p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],
                             mylegend, ncol=3)
                dev.off()
                p = plotMolProperties(preds, gains=gains, gof=gof)
                mylegend<-g_legend(p[[1]])
                for(i in 1:length(p)) p[[i]] = p[[i]] + theme(legend.position="none")
                pdf(file=paste(res_dir, "/predictions_molproperties_", analysis,"_",cutoff, ".pdf", sep=""), h= 20, w=24)
                if(gains & gof) {
                  grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]], p[[6]], mylegend, ncol=3)
                }else{
                  grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]], mylegend, ncol=3)
                }
                dev.off()
            }
        }
    }

    ## Read the cv_stats to get test sensitivity
    cv_stats = read.table(file=paste0(dataPath, "/cv_stat_summary.tsv"), header = T, sep="\t")
    cv_stats = cv_stats %>% subset(kernel==kern & type=="Sensitivity" & set=="test") %>% select(analysis, iterations, min, q1, median, mean, q3, max, var)
    grid_predictions_table = grid_predictions_table %>% left_join(cv_stats, by=c("analysis"))


    require(xlsx)
    wb<-createWorkbook(type="xlsx")

    TEXT_STYLE           <- CellStyle(wb) + Font(wb,  heightInPoints=14)
    TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
                              Alignment(rotation=45, horizontal="ALIGN_CENTER") +
                              Border(color="black", position=c("TOP", "LEFT","BOTTOM", "RIGHT"), pen=rep("BORDER_THIN",4))


    write_text_in_cell = function(sheet, rowIndex, text, textStyle){
      rows <-createRow( sheet, rowIndex=rowIndex )
      sheetText <-createCell( rows, colIndex=1 )
      setCellValue( sheetText[[1,1]], text )
      setCellStyle( sheetText[[1,1]], textStyle )
    }

    sheet <- createSheet(wb, sheetName = "predictions")

    text = with( grid_predictions_table, paste0("TRAINING SET. Number_of_genes = ", training_genes[1], " ; Samples = ", training_samples[1], " ; Total_entries = ", training_set_size[1]))
    write_text_in_cell(sheet, 1, text, TEXT_STYLE)

    text = with( grid_predictions_table, paste0("TEST SET.     Number_of_genes = ", prediction_genes[1], " ; Samples = ", prediction_samples[1], " ; Total_entries = ", prediction_set_size[1]))
    write_text_in_cell(sheet, 2, text, TEXT_STYLE)

#     ## Write unique numbers
#     text =  paste0("Unique novel drivers across models: ", length(unique(novel_drivers_genes_platt_all)))
#     write_text_in_cell(sheet, 3, text, TEXT_STYLE)
#     text =  paste0("Unique novel drivers across models (inNCG5): ", length(unique(novel_drivers_genes_platt_inNCG5_all)))
#     write_text_in_cell(sheet, 4, text, TEXT_STYLE)
#     text =  paste0("Unique novel drivers across models (rst): ", length(unique(novel_drivers_genes_platt_rst_all)))
#     write_text_in_cell(sheet, 5, text, TEXT_STYLE)
#     text =  paste0("Unique new samples covered across models: ", length(unique(new_samples_covered_all)))
#     write_text_in_cell(sheet, 6, text, TEXT_STYLE)

    to_remove = c("training_genes","training_samples",'training_set_size','prediction_genes','prediction_samples','prediction_set_size')

    rownames(grid_predictions_table) = NULL

    addDataFrame(grid_predictions_table[, !colnames(grid_predictions_table)%in%to_remove ], sheet, startRow=7, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, row.names = F)

    saveWorkbook(wb,file=paste(res_dir,"/grid_predictions_", kern, ".xlsx", sep=""))

    # write.xlsx(grid_predictions_table, file=paste(res_dir,"/grid_predictions_", kern, ".xlsx", sep=""))

#     ## Print unique numbers
#     cat(paste0("Unique novel drivers across models: ", length(unique(novel_drivers_genes_platt_all))), "\n")
#     cat(paste0("Unique novel drivers across models (inNCG5): ", length(unique(novel_drivers_genes_platt_inNCG5_all))), "\n")
#     cat(paste0("Unique novel drivers across models (rst): ", length(unique(novel_drivers_genes_platt_rst_all))), "\n")
#     cat(paste0("Unique new samples covered across models: ", length(unique(new_samples_covered_all))), "\n")


    ## Make a barplot for the recurrence of genes across models
#     t = table(table(unlist(grid_predictions)))
#     m = melt(t)
#     m[,1] = factor(m[,1])
#     save(grid_predictions, file=paste(res_dir, "/predictions_recurrence_",kern, ".Rdata" , sep=""))
#     pdf(file=paste(res_dir, "/predictions_recurrence_",kern, ".pdf" , sep=""), h= 20, w=24)
#     p = ggplot(m, aes_string(x=colnames(m)[1], y=colnames(m)[2])) +
#         geom_bar(stat="identity",col="black")+
#         geom_text(aes_string(label=colnames(m)[2]), position=position_dodge(width=0.9), vjust=-0.25, size=7)+
#         theme_boss()+theme(axis.line=element_line(colour="black"),
#                            axis.text.x=element_text(size = 12),
#                            axis.text.y=element_text(size = 12)) +
#         ggtitle("Recurrence of genes across models") +
#         xlab("") + ylab("Genes (n)")
#     print(p)
#     dev.off()

}

## Matteo's code for plotting systems-level properties
get_barplot = function(t, main=""){
    t = as.matrix(t)
    s = apply(t,1,sum)
    p = c()
    for(i in 1:ncol(t)) p[i] = fisher.test( matrix(c(t[,i],s),nc=2))$p.value
    names(p) = colnames(t)
    t = t/apply(t,1,sum)
    m = melt(t)
    m$p = p[match(m[,2], names(p))]
    m[,2] = factor(m[,2])
    p = unique(ddply(m, colnames(m)[2], mutate, value = value[which.max(value)], p=round(p[which.max(value)],3), pr = p[which.max(value)]))
    fdr = p.adjust(unique(p[,c(2,5)])$pr,'BH')
    names(fdr) = unique(p[,2])
    p$fdr = fdr[p[,2]]
    p$star = ""
    p$star[p$fdr<=0.1] = "*"
    p$print = paste0(p$p,p$star)
    p[,2] = factor(p[,2])
    r = ggplot(m, aes_string(fill=colnames(m)[1], x=colnames(m)[2], y=colnames(m)[3] )) +
        geom_bar(stat="identity", position=position_dodge(), col="black")+
        geom_text( data=p, aes_string(x = colnames(p)[2], y = colnames(p)[3], label=colnames(p)[8]), nudge_y = 0.1, size=4)+
        theme_boss()+theme(axis.line=element_line(colour="black"),
                           axis.text.x=element_text(size = 12),
                           axis.text.y=element_text(size = 12)) +
        scale_fill_manual(values=c("white","grey"))+ggtitle(main) +
        xlab("")
    r
}
get_boxplot = function(t, main="", log10=F){
    p = wilcox.test(t$value[t[,1]=="pred"], t$value[t[,1]=="rst"] )$p.value
    p <- format(p, scientific = TRUE)
    if (p < 0.05){
        p = paste(p, "*", sep="")
    }
    r = NULL
    if(log10==F){
        r = ggplot(t, aes_string(fill=colnames(t)[1], y=colnames(t)[3], x=colnames(t)[2] )) +geom_boxplot(col="black")+
            theme_boss() + theme(axis.line=element_line(colour="black"),
                                 axis.text.x=element_text(size = 12),
                                 axis.text.y=element_text(size = 12)) +
            scale_fill_manual(values=c("white","grey"))+ggtitle(paste0(main," p=",p))+
            xlab("")
    }else if (log10==T){
        r = ggplot(t, aes_string(fill=colnames(t)[1], y=colnames(t)[3], x=colnames(t)[2] )) +geom_boxplot(col="black")+
            theme_boss() + theme(axis.line=element_line(colour="black"),
                                 axis.text.x=element_text(size = 12),
                                 axis.text.y=element_text(size = 12)) +
            scale_fill_manual(values=c("white","grey"))+ggtitle(paste0(main," p=",p))+scale_y_log10()+
            xlab("")
    }
    r
}
g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

plotSysproperties = function(predictions, cutoff, tissue, gene_class="rst"){

    ## Get positive predictions
    if(gene_class=="all"){
        predictions = predictions %>% subset(label==TRUE & prob_platt>cutoff) %>% select(entrez, label) %>% unique
        print(nrow(predictions))
    }else{
        predictions = predictions %>% subset(label==TRUE & prob_platt>cutoff & gene_type==gene_class) %>% select(entrez, label) %>% unique
        print(nrow(predictions))
    }
    ## Load the file with gene properties
    load("/mnt/lustre/users/k1469280/mourikisa/data/geneProperties_final_mmImputed.Rdata")
    geneProperties=geneProperties_mmImputed
#     ## Get cancer genes to compare
#     dev_mode()
#     if(install_NCG) install(NCG_path)
#     library(ncglib)
#     geneInfo=get_geneinfo(version = "NCG5")
#     cancerGenes = get_cancer_genes(version = "NCG5")
#     ncg_connection = get_ncg_connection()
#     ppi   = ncg_connection %>% tbl("NCG51_interactions") %>% collect(n=Inf)
#     cancer_genes = as.data.frame(subset(get_cancer_genes(), primary_site==tissue))
#     ppi = subset(ppi, entrez%in%cancer_genes$entrez)
#     primary_inter = unique(ppi$entrez2)
#     dev_mode()

    ## Load NCG data from files instead of pulling them directly from the database (to run it in the cluster)
    load("/mnt/lustre/users/k1469280/mourikisa/data/geneInfoNCG5.Rdata")
    load("/mnt/lustre/users/k1469280/mourikisa/data/cancerGenesNCG5.Rdata")
    cancer_genes = as.data.frame(subset(cancerGenes, primary_site==tissue))
    load("/mnt/lustre/users/k1469280/mourikisa/data/ppiNCG51.Rdata")

    ## Fix primary interactors
    ppi = as.data.frame(ppi)
    ppi = ppi %>% group_by(entrez2) %>% summarise(primary_inter=length(unique(entrez))) %>%
        rename(entrez=entrez2)

    ## Fix gene info table from NCG
    geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
    ## Get a cancer gene with all the associated primary sites and cancer sites
    cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
        group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
                                       cancer_site=paste(unique(cancer_site), collapse=",")) %>%
        ungroup
    geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))
    geneProperties = geneProperties %>% rename(entrez=Entrez) %>% left_join(geneInfo, by=c("entrez"))

    ## Bring into the gene properties the predictions
    geneProperties = geneProperties %>% left_join(predictions, by=c("entrez"))

    if(gene_class=="all"){

      geneProperties = geneProperties %>% mutate(type=ifelse(label==TRUE & !is.na(label), "pred", "rst"))

    }else{
      geneProperties = geneProperties %>% subset(cancer_type==gene_class) %>% mutate(type=ifelse(label==TRUE & !is.na(label), "pred", "rst"))
    }

    mypredictions = geneProperties %>% subset(type=="pred") %>% nrow

    myrest = geneProperties %>% subset(type=="rst") %>% nrow

    #geneProperties = geneProperties %>% mutate(type=ifelse(grepl("blood", primary_site), "blood", "rst"))
    ## Bring into the gene properties the primary interactors
    geneProperties = geneProperties %>% left_join(ppi, by=c("entrez"))
    geneProperties$primary_inter[is.na(geneProperties$primary_inter)] = 0


    ## Convert properties to factors
    geneProperties$Genic = factor(geneProperties$Genic, levels=c("0","1"))
    geneProperties$WGD = factor(geneProperties$WGD, levels=c("0","1"))
    geneProperties$hub = factor(geneProperties$hub, levels=c("0","1"))
    geneProperties$central = factor(geneProperties$central, levels=c("0","1"))
    geneProperties$age = factor(geneProperties$age, levels=c("young","old"))
    geneProperties$origin = factor(geneProperties$origin, levels=c("LUCA","Eukaryotes","Opisthokonts","Metazoans","Vertebrates","Mammals","Primates"))
    geneProperties$exp.breadth = factor(geneProperties$exp.breadth, levels=c("Neverexpressed", "OneTissue",  "Selective","Middle", "AlwaysExpressed"  ))
    geneProperties$cancer_type = factor(geneProperties$cancer_type, levels=c("cgc",'can','rst'))
    geneProperties$type = factor(geneProperties$type, levels=c('pred','rst'))
    #geneProperties$type = factor(geneProperties$type, levels=c('blood','rst'))

    ## make the plots

    gimme_table <- function (t) {
      maximum = which(as.numeric(colnames(t))>9)
      if(length(maximum)>1){
        t = cbind(t[,1:(maximum[1]-1) ], "≥10"=apply(t[,maximum],1, sum))
      }
      t
    }

    p = list()
    p[[1]] = get_barplot(with(geneProperties, table( type, Genic) ), paste("Genic", "(", mypredictions, "vs", myrest, ")",sep="") )
    p[[2]] = get_barplot(with(geneProperties, table( type, inCDD) ), "in CDD" )
    p[[3]] = get_barplot( gimme_table(with(geneProperties, table( type, alldomains) )) , "all domains")
    p[[4]] = get_barplot( gimme_table(with(geneProperties, table( type, private_domains) )) , "private domains")
    p[[5]] = get_barplot( gimme_table(with(geneProperties, table( type, commondomains) )) , "common domains")
    p[[6]] = get_barplot(with(geneProperties, table( type, age) ), "age" )
    p[[7]] = get_barplot(with(geneProperties, table( type, origin) ), "origin" )
    p[[8]] = get_barplot(with(geneProperties, table( type, memberofcomplex) ), "Member of complex" )
    p[[9]] = get_barplot( gimme_table(with(geneProperties, table( type, High) )), "High" )
    p[[10]] = get_barplot( gimme_table(with(geneProperties, table( type, Low) )), "Low" )
    p[[11]] = get_barplot( gimme_table(with(geneProperties, table( type, Medium) )), "Medium" )
    p[[12]] = get_barplot( gimme_table(with(geneProperties, table( type, NotExpressed) )), "NotExpressed" )
    p[[13]] = get_barplot(with(geneProperties, table( type, exp.breadth) ), "exp.breadth" )
    p[[14]] = get_boxplot(melt(geneProperties, id.var="type", measure.vars = 'Length.fullrefseq' ),log10=T )
    p[[15]] = get_barplot(with(geneProperties, table( type, WGD) ), "WGD" )
    p[[16]] = get_boxplot(melt(geneProperties, id.var="type", measure.vars = 'degree' ), log10=T )
    p[[17]] = get_boxplot(melt(geneProperties, id.var="type", measure.vars = 'betweenness'),log10=T)
    p[[18]] = get_barplot(with(geneProperties, table( type, hub) ), "Hub" )
    p[[19]] = get_barplot(with(geneProperties, table( type, central) ), "Central" )
    p[[20]] = get_barplot( gimme_table( with(geneProperties, table( type, primary_inter) )) , paste0("Primary interactors of cancer genes of ", tissue))
    return(p)
}

## Plot properties
# pdf(file="predictions_properties.pdf", h= 16, w=20)
# grid.arrange(p[[1]],p[[3]],p[[6]],p[[7]],p[[8]],
#               p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]], ncol=3)
# dev.off()
# write.xlsx(geneProperties,"predictions_properties.xlsx", row.names=F, showNA=F)
## Assess the predictions (back-characterise them in terms of their features etc...)
## I am plotting also the features of the training set (no scaling)
# do.call("grid.arrange", c(plots, ncol=3))
plotMolProperties = function(predictions, gains = T, gof = T, gene_type="rst"){

    ## get numbers for false and true
    falses = predictions %>% subset(label==FALSE) %>% nrow
    trues = predictions %>% subset(label==TRUE) %>% nrow

    pl = function(m, s, main=""){
        m = data.frame(m)

        m$print = round((m$value)*100,digits = 1)
        m[,2] = factor(m[,2])
        p = ggplot(m, aes_string(fill=colnames(m)[1], x=colnames(m)[2], y=colnames(m)[3] )) +
            geom_bar(stat="identity", position=position_dodge(), col="black") +
            geom_text( data=m, aes_string(x = colnames(m)[2], y = colnames(m)[3],
                                          ymax=colnames(m)[3], label=colnames(m)[4]),
                       position = position_dodge(width=1), vjust=-0.7, size=7)+
            theme_boss()+theme(axis.line=element_line(colour="black"),
                               axis.text.x=element_text(size=12),
                               axis.text.y=element_text(size=12))+
            scale_fill_manual(values=c("grey","white"))+ggtitle(main)+
            scale_y_continuous(limits=c(0,1)) +
            xlab("")
    }

    if (gene_type=="rst"){predictions = predictions %>% subset(gene_type=="rst")}

    p = list()
    m = subset(predictions, no_TRUNC_muts!=0 & !((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)))
    if(nrow(m)>0) {
       m = m %>% group_by(label, no_TRUNC_muts) %>% summarise(value=n()) %>% ungroup %>% mutate(value=ifelse(label==FALSE, value/falses, value/nrow(predictions)))
    } else{
      m = data.frame(label=c(FALSE, TRUE), no_TRUNC_muts=c(1,1), value=c(0,0))
    }
    p = c(p, list(pl(m, main=paste0("Truncating mutations \n", predictions%>%select(entrez)%>%unique%>%nrow,
                                    " genes in ", predictions%>%select(sample)%>%unique%>%nrow, " samples = ", nrow(predictions), " entries" ))))

    m = subset(predictions, no_NTDam_muts!=0 & !((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)))
    if(nrow(m)>0) {
      m = m %>% group_by(label, no_NTDam_muts) %>% summarise(value=n()) %>% ungroup %>% mutate(value=ifelse(label==FALSE, value/falses, value/nrow(predictions)))
    } else{
      m = data.frame(label=c(FALSE, TRUE), no_NTDam_muts=c(1,1), value=c(0,0))
    }
    p = c(p, list(pl(m, main="Non-Truncating damaging mutations") ))

    if (gof){
        m = subset(predictions, no_GOF_muts!=0 & !((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)))
        if(nrow(m)>0) {
        m = m %>% group_by(label, no_GOF_muts) %>% summarise(value=n()) %>% ungroup %>% mutate(value=ifelse(label==FALSE, value/falses, value/nrow(predictions)))
        } else{
        m = data.frame(label=c(FALSE, TRUE), no_GOF_muts=c(1,1), value=c(0,0))
        }
        p = c(p, list(pl(m, main="Gain-of-function mutations") ))
    }

    m = subset(predictions, CNVLoss==1 & !((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)))
    if(nrow(m)>0) {
      m = m %>% group_by(label, CNVLoss) %>% summarise(value=n()) %>% ungroup %>% mutate(value=ifelse(label==FALSE, value/falses, value/nrow(predictions)))
    } else{
      m = data.frame(label=c(FALSE, TRUE), CNVLoss=c(1,1), value=c(0,0))
    }
    p = c(p, list(pl(m, main="CNV Losses") ))

    if(gains){
      m = subset(predictions, CNVGain==1 & !((CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)))
      if(nrow(m)>0) {
        m = m %>% group_by(label, CNVGain) %>% summarise(value=n()) %>% ungroup %>% mutate(value=ifelse(label==FALSE, value/falses, value/nrow(predictions)))
      } else{
        m = data.frame(label=c(FALSE, TRUE), CNVGain=c(1,1), value=c(0,0))
      }
      p = c(p, list(pl(m, main="CNV Gains") ))
    }

    m = subset(predictions, (CNVLoss==1 & Copy_number<2) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0) )
    if(nrow(m)>0) {
      m = m %>% group_by(label, CNVLoss) %>% summarise(value=n()) %>% ungroup %>% mutate(value=ifelse(label==FALSE, value/falses, value/nrow(predictions)))
    } else{
      m = data.frame(label=c(FALSE, TRUE), CNVLoss=c(1,1), value=c(0,0))
    }
    p = c(p, list(pl(m, main="CNV Losses & (TRUNC OR NTDam muts)")))

    return(p)
}
# pdf(file="predictions_molecular_properties.pdf", h= 16, w=20)
# grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]], ncol=3)
# dev.off()


################################################
# Feature Ranking with SVM-RFE
################################################
svmrfeFeatureRanking = function(training, mynu, mygamma, mydegree, kern){
  x = training%>%select(-type)%>%data.matrix()
  y = training%>%.$type
  n = ncol(x)

  survivingFeaturesIndexes = seq(1:n)
  featureRankedList = vector(length=n)
  rankedFeatureIndex = n

  while(length(survivingFeaturesIndexes)>0){
    #train the support vector machine
    if(kern=="linear"){
        svmModel = svm(x[, survivingFeaturesIndexes], y, nu = mynu, cachesize=500,  scale=F, type="one-classification", kernel=kern )
    }

    if (kern=="linear"){
        svm.model <- svm(x[, survivingFeaturesIndexes], y, kernel=kern, type="one-classification", scale=FALSE, nu=mynu)
    }else if (kern=="radial"){
        svm.model <- svm(x[, survivingFeaturesIndexes], y, kernel=kern, type="one-classification", scale=FALSE, gamma = mygamma, nu=mynu)
    }else if (kern=="polynomial"){
        ## probably changing degree
        svm.model <- svm(x[, survivingFeaturesIndexes], y, kernel=kern, type="one-classification", scale=FALSE, gamma = mygamma, nu=mynu, degree = mydegree)
    }else if (kern=="sigmoid"){
        svm.model <- svm(x[, survivingFeaturesIndexes], y, kernel=kern, type="one-classification", scale=FALSE, gamma = mygamma, nu=mynu)
    }

    #compute the weight vector
    w = t(svmModel$coefs)%*%svmModel$SV

    #compute ranking criteria
    rankingCriteria = w * w

    #rank the features
    ranking = sort(rankingCriteria, index.return = TRUE)$ix

    #update feature ranked list
    featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
    rankedFeatureIndex = rankedFeatureIndex - 1

    #eliminate the feature with smallest ranking criterion
    survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]]

  }

  return (featureRankedList)
}

##################################################
# Feature Ranking with Average Multiclass SVM-RFE
##################################################
svmrfeFeatureRankingForMulticlass = function(x,y){
  n = ncol(x)

  survivingFeaturesIndexes = seq(1:n)
  featureRankedList = vector(length=n)
  rankedFeatureIndex = n

  while(length(survivingFeaturesIndexes)>0){
    #train the support vector machine
    svmModel = svm(x[, survivingFeaturesIndexes], y, cost = 10, cachesize=500,  scale=F, type="C-classification", kernel="linear" )

    #compute the weight vector
    multiclassWeights = svm.weights(svmModel)

    #compute ranking criteria
    multiclassWeights = multiclassWeights * multiclassWeights
    rankingCriteria = 0
    for(i in 1:ncol(multiclassWeights))rankingCriteria[i] = mean(multiclassWeights[,i])

    #rank the features
    (ranking = sort(rankingCriteria, index.return = TRUE)$ix)

    #update feature ranked list
    (featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]])
    rankedFeatureIndex = rankedFeatureIndex - 1

    #eliminate the feature with smallest ranking criterion
    (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
    cat(length(survivingFeaturesIndexes),"\n")
  }

}

################################################
# This function gives the weights of the hyperplane
################################################
svm.weights<-function(model){
  w=0
  if(model$nclasses==2){
    w=t(model$coefs)%*%model$SV
  }else{    #when we deal with OVO svm classification
    ## compute start-index
    start <- c(1, cumsum(model$nSV)+1)
    start <- start[-length(start)]

    calcw <- function (i,j) {
      ## ranges for class i and j:
      ri <- start[i] : (start[i] + model$nSV[i] - 1)
      rj <- start[j] : (start[j] + model$nSV[j] - 1)

      ## coefs for (i,j):
      coef1 <- model$coefs[ri, j-1]
      coef2 <- model$coefs[rj, i]
      ## return w values:
      w=t(coef1)%*%model$SV[ri,]+t(coef2)%*%model$SV[rj,]
      return(w)
    }

    W=NULL
    for (i in 1 : (model$nclasses - 1)){
      for (j in (i + 1) : model$nclasses){
        wi=calcw(i,j)
        W=rbind(W,wi)
      }
    }
    w=W
  }
  return(w)
}

#################################################
# Feature weights for a range of hyperparameters
#################################################
getFeatureWeights = function(trainingPath, kern, nuRange){
  load(trainingPath)
  training = training %>% subset(type=="C")
  df = NULL
  count = 1
  for (mynu in nuRange){
    frl = svmrfeFeatureRanking(training, mynu, kern)
    features = names(training)[frl]
    dff = data.frame(features=features)
    dff$rank = 1:nrow(dff)
    colnames(dff) = c("features", paste("rank_", mynu, sep=""))
    if (count==1){
      df = dff
    }else{
      df = df %>% left_join(dff, by=c("features"))
    }
    count = count+1
  }
  return(df)
}

## I usually plot weights by doing:
# df %>% gather(model, rank, -features) %>%
#   ggplot(aes(x=model, y=factor(features))) +
#   geom_tile(aes(fill=rank)) + coord_equal()+
#   scale_fill_gradient(low = "red",high = "yellow") +
#   theme(
#     axis.text.x=element_text(angle=90, colour = "black"),axis.text.y=element_text(colour = "black"), plot.title = element_text(vjust=2)
#   ) +
#   xlab("") + ylab("") + ggtitle("Feature ranking for \n linear nu=0.5 (refined features)")


## Prepare tables for each cancer type and put them in the MySQL database
# for (ct in cancer_types){
#     cat(ct, "\n")
#     df = prepareTableML(df=total_table, cancer_type=ct)
#     write.table(df, file=paste("/home/mourikisa/novel_driver_prediction/Rdata/total_table_TCGA_19014_prepared_", ct, ".tsv", sep=""),
#                 row.names = F, col.names = F, sep = "\t", quote = F)
# }
prepareTableML <- function(df, cancer_type, geneProperties_dir="Rdata/geneProperties.Rdata",...){
    df = df %>% subset(Cancer_type==ct)
    ## Fix TPM, where there is NA put -1 cause in the database is numeric
    df$TPM[is.na(df$TPM)] = -1

    cat("Fixing copy number variations according to new cut-offs...", "\n")
    df <- df %>% mutate(CNV_type=ifelse(Copy_number>=4, "Gain", ifelse(Copy_number<2, "Loss", NA)))

    ## Replace numbers with names here
    df[,c("no_NSI_muts", "no_TRUNC_muts",
          "no_NTDam_muts", "no_NTDamFunction_muts",
          "no_NTDamCons_muts", "no_NTDamSC_muts",
          "no_GOF_muts")][is.na(df[,c("no_NSI_muts", "no_TRUNC_muts",
                                      "no_NTDam_muts", "no_NTDamFunction_muts",
                                      "no_NTDamCons_muts", "no_NTDamSC_muts",
                                      "no_GOF_muts")])] <- 0

    ## Join with gene Properties
    cat("Joining table with systems-level properties...", "\n")
    load(geneProperties_dir)
    df <- df %>% left_join(geneProperties, by=c("Entrez"))

    ## Convert categorical features to multiple factors
    cat("Performing cleaning of categorical variables...", "\n")
    ## CNV type
    df <- df %>%
        mutate(CNVGain=ifelse(is.na(CNV_type), 0, ifelse(CNV_type=="Gain",1, 0)),
               CNVLoss=ifelse(is.na(CNV_type), 0, ifelse(CNV_type=="Loss",1, 0))) %>%
        select(-CNV_type)

    ## Copy number (where copy number is NA, put copy number equal to 2)
    df$Copy_number[which(is.na(df$Copy_number))] <- 2

    ## Convert the length to bp instead of aa
    df <- df %>% mutate(Length.fullrefseq=Length.fullrefseq*3)

    ## expr_class
    df <- df %>%
        mutate(ExpT_ME=ifelse(is.na(expr_class), 0, ifelse(expr_class=="ME",1, 0)),
               ExpT_HE=ifelse(is.na(expr_class), 0, ifelse(expr_class=="HE",1, 0)),
               ExpT_LE=ifelse(is.na(expr_class), 0, ifelse(expr_class=="LE",1, 0)),
               ExpT_NE=ifelse(is.na(expr_class), 0, ifelse(expr_class=="notExpr",1, 0)),
               ExpT_NET=ifelse(is.na(expr_class), 0, ifelse(expr_class=="notExpr(T)",1, 0))) %>%
        select(-expr_class)

    ## age
    df <- df %>%
        mutate(old=ifelse(is.na(age), 0, ifelse(age=="old",1, 0)),
               young=ifelse(is.na(age), 0, ifelse(age=="young",1, 0))) %>%
        select(-age)

    ## origin
    df <- df %>%
        mutate(luca=ifelse(is.na(origin), 0, ifelse(origin=="LUCA",1, 0)),
               eukaryotes=ifelse(is.na(origin), 0, ifelse(origin=="Eukaryotes",1, 0)),
               metazoans=ifelse(is.na(origin), 0, ifelse(origin=="Metazoans",1, 0)),
               vertebrates=ifelse(is.na(origin), 0, ifelse(origin=="Vertebrates",1, 0)),
               opisthokonts=ifelse(is.na(origin), 0, ifelse(origin=="Opisthokonts",1, 0)),
               mammals=ifelse(is.na(origin), 0, ifelse(origin=="Mammals",1, 0)),
               primates=ifelse(is.na(origin), 0, ifelse(origin=="Primates", 1, 0))) %>% ## Initially omitted from the models - added here for reference
        select(-origin)

    ## exp.breadth.class
    df <- df %>%
        mutate(selective=ifelse(is.na(exp.breadth), 0, ifelse(exp.breadth=="Selective",1, 0)),
               always.expressed=ifelse(is.na(exp.breadth), 0, ifelse(exp.breadth=="AlwaysExpressed",1, 0)),
               middle=ifelse(is.na(exp.breadth), 0, ifelse(exp.breadth=="Middle",1, 0)),
               one.tissue=ifelse(is.na(exp.breadth), 0, ifelse(exp.breadth=="OneTissue",1, 0)),
               never.expressed=ifelse(is.na(exp.breadth), 0, ifelse(exp.breadth=="Neverexpressed",1, 0))) %>%
        select(-exp.breadth)

    ## Degree, Betweenness, Hub, Central have NAs
    df <- df %>% mutate(degree=ifelse(is.na(degree), 0, degree), betweenness=ifelse(is.na(betweenness), 0, betweenness),
                        hub=ifelse(is.na(hub), 0, hub), central=ifelse(is.na(central), 0, central))

    ## Before you add (chnage NAs in these columns to 0)
    df[,c("High", "Low", "Medium", "NotExpressed")][is.na(df[,c("High", "Low", "Medium", "NotExpressed")])] <- 0
    df <- df %>% ungroup %>% mutate(tot.tissues=High+Low+Medium)
    df <- data.frame(df)

    cat("Converting features to factors...", "\n")
    fcols <- c("Genic", "memberofcomplex",
               "WGD", "hub", "central", "CNVGain", "CNVLoss",
               "ExpT_ME", "ExpT_HE", "ExpT_LE", "ExpT_NE",
               "ExpT_NET", "old", "young", "luca", "eukaryotes",
               "metazoans", "vertebrates", "opisthokonts",
               "mammals", "selective", "always.expressed",
               "middle", "one.tissue", "never.expressed")
    cols <- which(colnames(df) %in% fcols)
    for(i in cols){
        df[,i] = factor(df[,i], levels = c(0,1))
    }

#     df <- df %>% mutate(key=paste(Cancer_type, Sample, Entrez, sep=".")) %>%
#         select(-Cancer_type, -Sample, -Entrez)
#     rnames <- df$key
#     df <- df %>% select(-key)
#     df <- data.frame(df)
#     row.names(df) <- rnames

    return(df)
}

## Function to extract training set (Vog genes and NCG FP)
## Note that dplyr::filter is used - subset is very slow
## You need to run this function before the training because after scaling you cannot filter based on mutations etc
## Pay extra attention when you subset for training (no_GOF_muts is NAs or 0s)
get_VogNCGFP_training <- function(df, fp_dir="/mnt/lustre/users/k1469280/mourikisa/data/NCG_false_positives.txt", ...){

    ## Make sure you load the unscaled table because only there you can subset
    ## based on 0 mutations...etc, the scaled one has no 0s
    ## Load the table and get the rownames, then subset the scaled ttable based on the rownames
    #load("")

    ## Get false positive genes
    false_positive_genes <- read.table(fp_dir, header = T, sep = "\t")
    false_positive_genes <- false_positive_genes %>% select(entrez, symbol)

    ## Training set extraction
    cat("Extracting training dataset...", "\n")
    training <- df %>% dplyr::filter((vogel=="Vog.Oncogene" | vogel=="Vog.TS") & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | CNVGain==1 | CNVLoss==1 | ExpT_NET==1))
    training <- training %>% mutate(type=ifelse(vogel=="Vog.Oncogene", "O", "T"))
    cat("Adding NCG False Positives...", "\n")
    trainingFP <- df %>% dplyr::filter(Entrez %in% false_positive_genes$entrez)
    #trainingFP <- trainingFP %>% mutate(rowname=paste(Cancer_type, Sample, Entrez, sep=".")) %>% select(-Cancer_type, -Sample, -Entrez)
    trainingFP <- trainingFP %>% mutate(type="N")

    ## Final training set
    training <- rbind(training, trainingFP)

    ## General housekeeping of the training set
    training$type <- as.factor(training$type)

    ## Exclude here any features that you dont want to include in the
    ## training/prediction
    #rnames <- training$rowname
    #training <- training %>% select(-TPM, -inCDD, -vogel, -rowname)
    #training = data.frame(training)
    #row.names(training) <- rnames
    #training = training %>% select(rowname,type)
    return(training)
}

## Load the total table first to avoid re-loading everytime
## Usually I load total_table_muts_cnv
mutProfile = function(df, entrez, symbol){

    mut_profile=list()
    lof = df %>% subset(Entrez==entrez & (!is.na(no_TRUNC_muts) | !is.na(no_NTDam_muts) | (CNV_type=="Loss" & Copy_number==0))) %>% select(Cancer_type, Sample) %>%
        unique %>% group_by(Cancer_type) %>% summarise(n=n()) %>% ungroup %>% data.frame()

    lof2samples = df %>% select(Cancer_type, Sample) %>% unique %>% group_by(Cancer_type) %>% summarise(total_pats=n()) %>%
        ungroup %>% left_join(lof, by=c("Cancer_type"))

    ##Make the NAs, 0s
    lof2samples$n[is.na(lof2samples$n)] <- 0
    lof2samples = lof2samples %>% mutate(label=paste(Cancer_type, "(", total_pats, ")", sep=""), freq=(n/total_pats)*100)

    ## sort based on the freq
    lof2samples$label = factor(as.character(lof2samples$label), level = lof2samples$label[order(lof2samples$freq,decreasing = F)],ordered = T)

    #pdf(file="/Users/fc-8s-imac/Desktop/smarca4_LOF_TCGA.pdf",height = unit(8.27, "inches"), width = unit(11.69, "inches"), useDingbats = F)
    p = ggplot(lof2samples%>%arrange(freq), aes(x=label, y=freq)) + geom_bar(stat="identity") + theme_boss() +
        theme(axis.text.x=element_text(angle=90),
              axis.title.y=element_text(size = 12)) +
        ylab(paste("% samples with LOF mutation in ", symbol, sep="")) + xlab("") +
        ggtitle(paste("Percentage of samples with LOF mutations in ",symbol, " across \n 31 cancer types", sep=""))

    mut_profile[["plot"]] = p
    mut_profile[["table"]] = lof2samples
    #dev.off()
    #write.table(lof2samples, file="/Users/fc-8s-imac/Desktop/smarca4_LOF_TCGA.tsv", quote=F, row.names = F, sep="\t")

    return(mut_profile)
}

excludeTraining <- function(training="", df){
    load(training)
    df <- df[!(row.names(df) %in% row.names(df)),]
    return(df)
}

get_TCGA_cancer_types <- function(){
    cancer_types <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM",
                      "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD",
                      "LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM",
                      "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
    return(cancer_types)
}

get_TCGA_cancer_types_exprWithNormal <- function(){
    ctWithNormals <- c("BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA",
                       "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD",
                       "LUSC", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM",
                       "STAD", "THCA", "THYM", "UCEC")
    return(ctWithNormals)
}

################################################################
# check sensitivity for final predictions
# when CGCs altered in samples are considered True positives
# I will do it when consider different cuttofs in proabilities
# 95%, 96%, 97%, 98%, 99%, 100%
################################################################
checkPredictions = function(model, validation_set, geneInfo, plattTraining){

    ## Fix gene info table from NCG
    geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
    colnames(geneInfo) = c("entrez", "geneType")
    ## Make predictions and get the associated probability
    predictions = predict(svm.model, validation_set%>%select(-type), decision.values = TRUE)
    preds_platt = plattScaling(plattTraining, predictions)
    colnames(preds_platt) = c("dv", "label", "prob_platt")
    preds_binning = binning(svm.model, predictions)
    colnames(preds_binning) = c("dv", "label", "prob_binning")
    preds = preds_platt %>% tibble::rownames_to_column() %>%
        left_join(preds_binning%>%tibble::rownames_to_column(), by=c("rowname","dv", "label"))
    preds = preds %>%
        separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>%
        mutate(entrez=as.numeric(entrez)) %>%
        left_join(geneInfo, by=c("entrez"))

    cgcs_obs = preds %>% subset(geneType=="cgc") %>% nrow
    cgcs_unq = preds %>% subset(geneType=="cgc") %>% select(entrez) %>%
        unique %>% nrow
    samples = preds %>% select(sample) %>% unique %>% nrow
    cat("Samples: ", samples, "\n")
    cat("CGC genes: ", cgcs_obs, "\n")
    cat("CGC unique genes: ", cgcs_unq, "\n")

    platt_cuttoff = 0.95
    platt_sens_list_obs = NULL
    while (platt_cuttoff<=1){
        ## Sensitivity overall (observations)
        cat("Cut-off: ", platt_cuttoff, "\n")
        cat(paste("CGC observations: ", cgcs_obs, sep=""), "\n")
        genes = preds %>% subset(prob_platt > platt_cuttoff & geneType=="cgc") %>% nrow
        cat(paste("Observations with Platt's p > ", platt_cuttoff,": ", genes, sep=""), "\n")
        sens = genes/cgcs_obs
        platt_sens_list_obs = c(platt_sens_list_obs, sens)
        cat("Sensitivity (observations): ", sens, "\n")

        ## Sensitivity unique genes
        cat(paste("CGC unique genes: ", cgcs_unq, sep=""), "\n")
        genes = preds %>% subset(prob_platt > platt_cuttoff) %>%
            subset(geneType=="cgc") %>% select(entrez) %>% unique %>% nrow
        cat(paste("Unique genes with Platt's p > ", platt_cuttoff,": ", genes, sep=""), "\n")
        sens = genes/cgcs_unq
        cat("Sensitivity (unique CGC genes): ", sens, "\n")

        cat("\n")
        platt_cuttoff = platt_cuttoff + 0.01
    }
    return(platt_sens_list_obs)

}



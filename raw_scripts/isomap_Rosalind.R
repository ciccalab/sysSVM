rm(list=ls())

## Install packages
pkgsR  =  c('plyr',
            'tidyr',
            'xlsx',
            'dplyr',
            'Rtsne')
for (pkgR in pkgsR){
    if (!pkgR %in% rownames(installed.packages())) { 
        install.packages(pkgR, lib="~/R/x86_64-pc-linux-gnu-library/3.3", repos='https://cran.ma.imperial.ac.uk/', dependencies=TRUE, INSTALL_opts = c('--no-lock'))
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

#update.packages(repos=biocinstallRepos(), ask=FALSE)

rm(pkg,pkgR,pkgs,pkgsR)
## Load dplyr always last
detach("package:dplyr", unload=TRUE)
library(dplyr)

## Define command line arguments
args = commandArgs(trailingOnly = TRUE)

## Parameters
CANCER_TYPE = args[1]
TISSUE_NCG = args[2]
TISSUE_GTEX = args[3] ## Code commented out so it is not used at the moment
GOF=T
GAINS=T
if (args[4]=="F" | args[4]=="FALSE"){
    GAINS=F
}
CANCER_DIR = paste0(args[5], "/", CANCER_TYPE)
CONFIG_PATH = args[6]

ISOMAP=T
if (args[7]=="F"  | args[7]=="FALSE"){
    ISOMAP=F
}

#source(CONFIG_PATH)
 
load(paste0(CANCER_DIR, "/scores.Rdata"))


if (ISOMAP){
     
    ## Get features and scores (for dimensionality reduction)
    load(paste0(CANCER_DIR, "/training_set.Rdata"))
    load(paste0(CANCER_DIR, "/validation_set.Rdata"))

    cohort = rbind(training%>%tibble::rownames_to_column()%>%separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.") %>% data.frame(), 
                   validation%>%tibble::rownames_to_column()%>%separate(rowname, into=c("cancer_type", "sample", "entrez"), sep="\\.")%>%mutate(type="P")%>%data.frame())
    
    feature_names = names(validation)
    cohort[,feature_names] = sapply(cohort[,feature_names], function(x) as.numeric(as.character(x))) ## as.character is used to maintain factor's levels 
    
    cohort = cohort %>% mutate(entrez=as.numeric(entrez)) %>% left_join(scores[["scores"]])
    cohort$gene_type[is.na(cohort$gene_type)] = "cgc" ## because cgc is not in the scores file
    
    ## ------
    ## tsne
    ## ------
    tsne_data = Rtsne(as.matrix(cohort[,feature_names]), check_duplicates = FALSE)
    
    tsne_coords = tsne_data$Y %>% as.data.frame()
    colnames(tsne_coords) = c("tSNE1", "tSNE2")
    cohort = cbind(cohort, tsne_coords)
    
    ggplot(cohort, aes(tSNE1, tSNE2)) + 
        geom_jitter(data=cohort%>%subset(type=="P"), aes(shape=type, colour=score), width = 0.01, height = 0) +
        geom_jitter(data=cohort%>%subset(type=="C"), aes(tSNE1, tSNE2, colour=type), width = 0.01, height = 0) +
        geom_density2d(data=cohort%>%subset(type=="C"),contour = T, aes(colour = "#2b8cbe")) + ## density based on training
        scale_color_manual(values = c("blue", "black", "grey50")) +
        #geom_density2d(data=cohort%>%subset(type=="P"),contour = T, aes(colour = "#636363")) +
        theme_boss() + 
        theme(axis.line.x = element_line(color="black"),
              axis.line.y = element_line(color="black"))
    
    require(flowCore)
    require(spade)
    require(igraph)
    output_dir = "~/athena/data/OAC/129_OAC/dimensionality_reduction/spade_out"
    write.FCS(flowFrame(as.matrix(cohort[,c(feature_names, "tSNE1", "tSNE2")])), filename = "~/athena/data/OAC/129_OAC/dimensionality_reduction/OAC_tSNE.fcs")
    SPADE.driver("~/athena/data/OAC/129_OAC/dimensionality_reduction/OAC_tSNE.fcs", out_dir = output_dir, k=100, cluster_cols = c("tSNE1", "tSNE2"))
    mst_graph = igraph:::read.graph(paste(output_dir, "mst.gml", sep=.Platform$file.sep), format="gml")
    SPADE.plot.trees(mst_graph, output_dir, out_dir = paste(output_dir, "pdf", sep=.Platform$file.sep), layout = igraph:::layout.kamada.kawai(mst_graph))
    
    ## Plot
    plot(mst_graph, edge.arrow.size=.4,vertex.label=NA, vertex.size=5)
    
    ## Now I basically need to read SPADE's output and draw the heatmap to identify regions of the map with certain properties
    
    
    ## ------
    ## Isomap
    ## ------
    cat("Running Isomap...", "\n")
    iso_dim2 = Isomap(data=as.matrix(cohort[,feature_names]), dims=2)
    iso_dim2 = iso_dim2 %>% data.frame()
    colnames(iso_dim2) = c("X1", "X2")
    
    ## Add the coordinates
    export = cbind(cohort, iso_dim2)
    
    write.xlsx(export%>%data.frame(), file=paste0(CANCER_DIR, "/scores_plus_ISOMAP.xlsx"), row.names = F)
    message("Isomap run finished")
}else{
    message("Finished without Isomap")
}




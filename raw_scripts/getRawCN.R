## Function to get the raw copy number from ascat's output
library(dplyr)
library(tidyr)


getCNRaw = function(cnvs){
    ## There are three locations for the data
    loc1 = "/home/mourikisa/data/OAC/129_OAC/ascat/"
    loc2 = "/home/mourikisa/data/OAC/71_OAC/ascat/"
    loc3 = "/home/mourikisa/data/OAC/87_OAC/66_ICGC/ascat/"
    
    res = NULL
    count = 1
    for(sample_id in unique(cnvs$sample)){
        
        cat(paste0(count, ". ", sample_id), "\n")
        loc = NULL
        if(sample_id%in%list.files(loc1)){loc=loc1}
        if(sample_id%in%list.files(loc2)){loc=loc2}
        if(sample_id%in%list.files(loc3)){loc=loc3}
        
        if(is.null(loc)){stop("No sample found.")}
        
        ## Get the file with the raw copy number
        pipe = list.dirs(paste0(loc, sample_id), recursive = F) ## version of pipeline not always the same
        pipe = pipe[!grepl("parsing_and_annotation", pipe)]
        if(length(pipe[grepl("g1k", pipe)])>0){pipe = pipe[grepl("g1k", pipe)]}
        if(length(pipe)>1){pipe = pipe[grepl("g1k", pipe) & grepl("1.5", pipe)]} ## Get the latest version of the pipeline
        if(length(pipe)==0){stop("No pipes found")}
        pipe = unlist(strsplit(pipe, "/"))
        pipe = pipe[length(pipe)]

        temp_cnv_fn = paste0(loc, sample_id, "/", pipe, "/results/")
        raw_fn = grep("copynumber.txt", list.files(temp_cnv_fn, full.names = T), value = T)
        
        raw = read.table(raw_fn, sep = "\t", header = T)
        
        ## Finally get the regions we care about
        export = cnvs %>% subset(sample==sample_id)
        d = cnvs %>% subset(sample==sample_id) %>% select(Chromosome, Start, End) %>% unique %>% data.frame()
        
        cat(paste0("Regions: ", nrow(d) ), "\n")
        
        d$Raw_CN = apply(d, 1, function(x){
            chr = gsub("chr", "", x[1])
            start = as.numeric(x[2])
            end = as.numeric(x[3])
            
            cn_raw = raw %>% subset(Position>=start & Position<=end & Chromosome==chr) %>% .$Raw.copy.number %>% mean()
            return(cn_raw)
        })
        
        cat(nrow(export), "\n")
        export = export %>% left_join(d)
        cat(nrow(export), "\n")
        res = rbind(res, export)
        
        count = count + 1
    }
    
    return(res)
    
}
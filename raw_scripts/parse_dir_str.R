#!/usr/bin/Rscript

## utility script to create the commands to for the mutation annotation part of the pipeline

## Change the main directory of samples accordingly
mainDir = "~/rosalind_lustre/mourikisa/data/OAC/87_OAC/21_literature/"

if(!dir.exists(paste0(mainDir, "/Logs"))){
    dir.create(paste0(mainDir, "/Logs"))
}

sample_dirs = strsplit(list.dirs(paste0(mainDir, "/strelka"), recursive = F), "/")
sample_dirs = unlist(lapply(sample_dirs, function(x) x[length(x)]))
mut_fn = "/strelka/filtered_results/TBA.snp.pass.vcf"
indel_fn = "/strelka/results/passed.somatic.indels.vcf"
cnv_fn = "/results/TBA.copynumber.caveman.csv"
stats_fn = "/results/TBA.samplestatistics.csv"
qsub_cmd = "qsub -N mut_annot -pe smp 1 -l h_vmem=12G "
cmds = NULL
for(s in sample_dirs){
    ##muts
    cat(s, "\n")
    pipe = list.dirs(paste0(mainDir, "/strelka/", s), recursive = F) ## version of pipeline not always the same
    pipe = pipe[!grepl("parsing_and_annotation", pipe)]
    ## The following 2 lines is for ICGC (not applicable to all other cohorts) - comment them in for the main cohort
    #if(length(pipe[grepl("g1k", pipe)])>0){pipe = pipe[grepl("g1k", pipe)]}
    #if(length(pipe)>1){pipe = pipe[grepl("g1k", pipe) & grepl("1.5", pipe)]} ## Get the latest version of the pipeline
    if(length(pipe)==0 | length(pipe)>1){stop("No or multiple pipes found in strelka")}
    pipe = unlist(strsplit(pipe, "/"))
    pipe = pipe[length(pipe)]
    temp_m_fn = paste0(pipe, gsub("TBA", s, mut_fn))
    temp_i_fn = paste0(pipe, indel_fn)
    
    ## cnvs
    pipe = list.dirs(paste0(mainDir, "/ascat/", s), recursive = F) ## version of pipeline not always the same
    pipe = pipe[!grepl("parsing_and_annotation", pipe)]
    #if(length(pipe[grepl("g1k", pipe)])>0){pipe = pipe[grepl("g1k", pipe)]}
    #if(length(pipe)>1){pipe = pipe[grepl("g1k", pipe) & grepl("1.5", pipe)]} ## Get the latest version of the pipeline
    if(length(pipe)==0 | length(pipe)>1){stop("No or multiple pipes found in ascat")}
    pipe = unlist(strsplit(pipe, "/"))
    pipe = pipe[length(pipe)]
    ## Not all batches with the same pattern of file names
    #temp_cnv_fn = paste0(pipe, gsub("TBA", paste0(s, "__", pipe, "__"), cnv_fn))
    #temp_stats_fn = paste0(pipe, gsub("TBA", paste0(s, "__", pipe, "__"), stats_fn))
    temp_cnv_fn = paste0(pipe, gsub("TBA", s, cnv_fn))
    temp_stats_fn = paste0(pipe, gsub("TBA", s, stats_fn))
    
    
    ##svs
    # pipe = list.dirs(paste0(mainDir, "/manta/", s), recursive = F) ## version of pipeline not always the same
    # pipe = pipe[!grepl("parsing_and_annotation", pipe)]
    # if(length(pipe[grepl("g1k", pipe)])>0){pipe = pipe[grepl("g1k", pipe)]}
    # if(length(pipe)>1){pipe = pipe[grepl("g1k", pipe) & grepl("1.5", pipe)]} ## Get the latest version of the pipeline
    # if(length(pipe)==0){stop("No pipes found")}
    # 
    # ## grab the file with manta calls because not all files have similar names
    # temp_sv_fn = list.files(paste0(pipe, "/results/filtered"), recursive = F, pattern="manta.filtered")
    # if(length(temp_sv_fn)>1){stop("Multiple or no manta calls found")}
    # if(length(temp_sv_fn)<1){
    #     print(paste0("No files found: ", s))
    #     next    
    # }
    # 
    # pipe = unlist(strsplit(pipe, "/"))
    # pipe = pipe[length(pipe)]
    
    
    ## Fix mainDir for the cluster
    mD = gsub("~/rosalind_lustre/", "/mnt/lustre/users/k1469280/", mainDir)
    
    temp_cmd = paste0(
        qsub_cmd,
        "-o ",
        mD, "/Logs/", s, ".out -e ",
        mD, "/Logs/", s, ".err ./../../data_parsing.sh",
        " -s ", mD, "/strelka/", s,
        " -m ", mD, "/strelka/", s, "/",temp_m_fn, 
        " -i ", mD, "/strelka/", s, "/", temp_i_fn, 
        " -c ", mD, "/ascat/", s, "/", temp_cnv_fn,
        " -t ", mD, "/ascat/", s, "/", temp_stats_fn,
        #" -v ", mD, "/manta/", s, "/", pipe, "/results/filtered/", temp_sv_fn,
        " -g hg19"
    )
    cmds = rbind(cmds, temp_cmd)
}

write.table(cmds, file=paste0(mainDir, "/data_parsing_cmds.txt"), quote = F, row.names = F, col.names = F)



## check sample names across data types
sample_dirs_mut = strsplit(list.dirs(paste0(mainDir, "/strelka"), recursive = F), "/")
sample_dirs_mut = unlist(lapply(sample_dirs_mut, function(x) x[length(x)]))
sample_dirs_cnv = strsplit(list.dirs(paste0(mainDir, "/ascat"), recursive = F), "/")
sample_dirs_cnv = unlist(lapply(sample_dirs_cnv, function(x) x[length(x)]))
sample_dirs_sv = strsplit(list.dirs(paste0(mainDir, "/manta"), recursive = F), "/")
sample_dirs_sv = unlist(lapply(sample_dirs_sv, function(x) x[length(x)]))
library(gplots)
venn(list(mut=sample_dirs_mut, cnv=sample_dirs_cnv, sv=sample_dirs_sv))




## and some checking of the run
sample_dirs = list.dirs(mainDir, recursive = F)
count = 0
failed = 0
samples_failed = NULL
for(s in sample_dirs){
    cat(s, "\n")
    if(file.exists(paste0(s, "/mutation_annotation.err"))){
        con <- file(paste0(s, "/mutation_annotation.err"))
        lines = readLines(con)
        if(length(lines)>10){
            cat(paste(lines[(length(lines)-10):length(lines)], collapse = "\n"), "\n\n")
        }else{
            cat(paste(lines, collapse = "\n"), "\n\n")
        }
        
        if(sum("Execution halted"%in%lines)>0){
            samples_failed = c(samples_failed, s)
            failed = failed +1
        }
        close(con)
        count = count +1
    }
}
cat(paste0(count, " files have started"), "\n")
cat(paste0(failed, " files have failed"))
print(samples_failed)




###############################################################################
## Delete left-over files during trials (Just in case I need a reset)
mainDir = "~/athena/data/OAC/71_OAC/strelka"
sample_dirs = list.dirs(mainDir, recursive = F)
oldwd = getwd()
for(s in sample_dirs){
    setwd(paste0(s, "/parsing_and_annotation"))
    unlink("annovarnsyn.onco")
    unlink("annovarsyn.onco")
    unlink("annovaroncodriveclust-results.tsv")
}
setwd(oldwd)
###############################################################################






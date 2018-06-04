rm(list=ls())

## Define command line arguments
args = commandArgs(trailingOnly = TRUE)

## Parameters
CANCER_TYPE = args[1]
MYSQL=F
if (args[2]=="T" | args[2]=="TRUE"){
    MYSQL=T
}
INPUT_FN = args[3]

EXCLUDED_FEATURES = unlist(strsplit(args[4], split = ","))
## Get gains from command line
GAINS=T
if (args[5]=="F" | args[5]=="FALSE"){
    GAINS=F
}

SOURCE_DIR=args[6]
## Result dir - if not exist create it otherwise error DO NOT OVERWRITE
RES_DIR = paste0(args[7], "/", CANCER_TYPE)
if (dir.exists(RES_DIR)){
    stop("Results directory already exists!")
}else{
    dir.create(RES_DIR)
}

CONFIG_PATH = args[8]

## this is for scaling purposes (to scale with the mean and sd of the main cohort - we assume that they are representative)
MEANS_SDS = args[9]
if (MEANS_SDS=="F" | MEANS_SDS=="FALSE"){
    MEANS_SDS=NULL
}else {
    load(MEANS_SDS)
    MEANS_SDS = mean_sd
}

NU = seq(0.05, 0.9, 0.05)
GAMMA = 2^seq(-7,4)
DEGREE = c(3, 4, 9)
NCORE = 3

## Analysis pipeline
source(CONFIG_PATH)

write.config(EXCLUDED_FEATURES, RES_DIR)

createDescribeTrainingCGC(cancer_type=CANCER_TYPE, res_dir=RES_DIR, gains=GAINS, mysql=MYSQL, input_fn=INPUT_FN, means_sds = MEANS_SDS)

runCommands(CANCER_TYPE, "linear", NU, 0, 0,
            save_dir=RES_DIR,
            fname="jobs_linear.txt", ncore=NCORE,
            res_dir=RES_DIR,
            source_dir=SOURCE_DIR)

runCommands(CANCER_TYPE, "radial", NU, GAMMA, 0,
            save_dir=RES_DIR,
            fname="jobs_radial.txt", ncore=NCORE,
            res_dir=RES_DIR,
            source_dir=SOURCE_DIR)

runCommands(CANCER_TYPE, "polynomial", NU, GAMMA, DEGREE,
            save_dir=RES_DIR,
            fname="jobs_polynomial.txt", ncore=NCORE,
            res_dir=RES_DIR,
            source_dir=SOURCE_DIR)

runCommands(CANCER_TYPE, "sigmoid", NU, GAMMA, 0,
            save_dir=RES_DIR,
            fname="jobs_sigmoid.txt", ncore=NCORE,
            res_dir=RES_DIR,
            source_dir=SOURCE_DIR)
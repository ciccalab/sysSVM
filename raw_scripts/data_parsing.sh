#!/bin/sh
## Set the working Directorty to be the current one, i.e. where you are submitting the script
##$ -cwd
##$ -j y

## Set the SHELL type to be sh##
#$ -S /bin/sh

## Set the Parallel Environment (in this case mpi)
##$ -pe smp 1


while getopts ':s:m:i:c:t:v:g:' flag; do
  case "${flag}" in
    s) sample_dir="${OPTARG}" ;;
    m) snvs="${OPTARG}" ;;
    i) indels="${OPTARG}" ;;
    c) cnvs="${OPTARG}" ;;
    t) stats="${OPTARG}" ;;
    ##v) svs="${OPTARG}" ;;
    g) genome="${OPTARG}" ;;
  esac
done

echo "Sample directory: ${sample_dir}"


## CONFIGURATION
## --------------

module load general/R/3.2.1
module load bioinformatics/bedtools2/2.25.0

Rscript='/mnt/lustre/users/k1469280/mourikisa/data/OAC/data_parsing.R'

## Run script
#Rscript $Rscript ${sample_dir} ${snvs} ${indels} ${cnvs} ${stats} ${svs} ${genome}
Rscript $Rscript ${sample_dir} ${snvs} ${indels} ${cnvs} ${stats} ${genome}

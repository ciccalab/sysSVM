#!/bin/sh
## Set the working Directorty to be the current one, i.e. where you are submitting the script
#$ -cwd
#$ -j y

## Set the SHELL type to be sh##
#$ -S /bin/sh

## Set the Parallel Environment (in this case mpi)
##$ -pe mpi 1

#$ -e /users/k1469280/Logs/
#$ -o /users/k1469280/Logs/

## CONFIGURATION
## --------------

module load bioinformatics/R/3.3.3

## Prepare parameters
ct=$1
tissue_ncg=$2
tissue_gtex=$3
gains=$4
res_dir=$5
config_path=$6
prediction_set=$7
prediction_set_ns=$8
stats=$9
columns=$10

Rscript='predict_new_samples.R'

# ------------
# RUN PIPELINE
# ------------

## STEP 01: Prepare result file, training/validation set and training commands
Rscript $Rscript $ct $tissue_ncg $tissue_gtex $gains $res_dir $config_path $prediction_set $prediction_set_ns $stats $columns



#!/bin/sh
## Set the working Directorty to be the current one, i.e. where you are submitting the script
#$ -cwd
#$ -j y

## Set the SHELL type to be sh##
#$ -S /bin/sh

## Set the Parallel Environment (in this case mpi)
##$ -pe mpi 1

#$ -e /home/mourikisa/Logs/
#$ -o /home/mourikisa/Logs/


##$ -M matteo.cereda@kcl.ac.uk
##$ -m e

## CONFIGURATION
## --------------

module load R

## Prepare parameters
n=$1
times=$2
save_dir=$3
chunk=$4

Rscript='randomise_GSEA.R'

# ----
# RUN
# ----

Rscript $Rscript $n $times $save_dir $chunk


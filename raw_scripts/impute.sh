#!/bin/sh
## Set the working Directorty to be the current one, i.e. where you are submitting the script
#$ -cwd
#$ -j y

## Set the SHELL type to be sh##
#$ -S /bin/sh

## Set the Parallel Environment (in this case mpi)
##$ -pe mpi 1

#$ -e /home/mourikisa/novel_driver_prediction/Logs/
#$ -o /home/mourikisa/novel_driver_prediction/Logs/

##$ -e /home/ceredam/novel_driver_prediction/Logs/
##$ -o /home/ceredam/novel_driver_prediction/Logs/


##$ -M matteo.cereda@kcl.ac.uk
##$ -m e

## CONFIGURATION
## --------------

module load R

## Prepare parameters
iteration=$1

Rscript='impute.R'

## Run imputation
Rscript $Rscript $iteration




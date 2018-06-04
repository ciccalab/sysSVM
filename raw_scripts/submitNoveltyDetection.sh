#!/bin/sh
## Set the working Directorty to be the current one, i.e. where you are submitting the script
#$ -cwd
#$ -j y

## Set the SHELL type to be sh##
#$ -S /bin/sh

## Set the Parallel Environment (in this case mpi)
##$ -pe mpi 1

##$ -e /home/mourikisa/novel_driver_prediction/Logs/
##$ -o /home/mourikisa/novel_driver_prediction/Logs/

##$ -e /home/ceredam/novel_driver_prediction/Logs/
##$ -o /home/ceredam/novel_driver_prediction/Logs/


##$ -M matteo.cereda@kcl.ac.uk
##$ -m e

## CONFIGURATION
## --------------

module load bioinformatics/R/3.3.3

ct=$1
kern=$2
mynu=$3
mygamma=$4
mydegree=$5
resDir=$6
script_path=$7
Rscript='/modelTrainingCV_NoveltyDetection.R'

# ---
# RUN
# ---
Rscript $script_path$Rscript $ct $kern $mynu $mygamma $mydegree $resDir $path

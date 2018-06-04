#!/bin/sh
## Set the working Directorty to be the current one, i.e. where you are submitting the script
#$ -cwd
#$ -j y

## Set the SHELL type to be sh##
#$ -S /bin/sh

## Set the Parallel Environment (in this case mpi)
##$ -pe mpi 1

##$ -e /users/k1469280/Logs/
##$ -o /users/k1469280/Logs/



##$ -M matteo.cereda@kcl.ac.uk
##$ -m e

## CONFIGURATION
## --------------

module load bioinformatics/R/3.3.3

## Prepare parameters
ct=$1
mysql=$2
input_fn=$3
excluded_features=$4
gains=$5
source_dir=$6
res_dir=$7
config_path=$8
mean_sd=$9

Rscript='prepare.R'

# ------------
# RUN PIPELINE
# ------------

## STEP 01: Prepare result file, training/validation set and training commands
Rscript $Rscript $ct $mysql $input_fn $excluded_features $gains $source_dir $res_dir $config_path $mean_sd

## STEP 02: Run training commands
cd $res_dir/$ct
chmod +x jobs_linear.txt
chmod +x jobs_radial.txt
chmod +x jobs_polynomial.txt
chmod +x jobs_sigmoid.txt



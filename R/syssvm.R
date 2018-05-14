#' Run sysSVM
#'
#' \code{syssvm} uses e1071 to run an one-class SVM
#'
#'
#' @param input.file path to input file
#' @param output.dir path to output directory where sysSVM results will be written
#' @param exclude.features features to be excluded
#' @param models check
#' @param scaling.factors object with sysSVM scaling factors. Use this to make predictions using already-made sysSVM models
#' @param cv k-fold cross-validation [default is 3]
#' @param iters check
#' @param kernels check
#' @param cores check
#' @param check check
#' @param reproduce set to TRUE to reproduce previously published results
#'
#'
#'
#'
#'
#'
#'

syssvm <- function(input.file=NULL, output.dir=NULL, exclude.features=NULL, models=NULL,
                   scaling.factors=NULL, cv=3, iters=100,
                   kernels=c("linear", "polynomial", "radial", "sigmoid"),
                   cores=2, ncg.tissue.name=NULL, reproduce=FALSE) {

  if(is.null(input.file) | is.null(output.dir)){
    stop("sysSVM ERROR: Input/Output file unknown, please provide input file and output directory")
  }

  trainingMode = TRUE
  if(!is.null(models) & !is.null(scaling.factors)){
    cat("Training mode disabled, using given models to predict cancer genes...", "\n")
    trainingMode = FALSE
  }

  source("R/dependencies.R")
  source("R/functions.R")

  ## Create training and prediction sets
  # createDescribeTrainingCGC(input.file=input.file, output.dir=output.dir, exclude.features=exclude.features,
  #                           trainingMode=trainingMode, models=models, scaling.factors=scaling.factors)
  #
  #
  # runNoveltyDetection(output.dir=output.dir, cv=cv, iters=iters, kernels=kernels, cores=cores)

  scoreGenes(ncg.tissue.name=ncg.tissue.name, output.dir=output.dir,
             exclude.features=exclude.features, reproduce=reproduce)

}

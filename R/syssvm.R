#' Run sysSVM
#'
#' \code{syssvm} uses e1071 to run an one-class SVM
#'
#'
#' @param object an object, whose class determines
#'
#'
#'
#'

syssvm <- function(input.file=NULL, output.dir=NULL, exclude.features=NULL, models=NULL,
                   scaling.factors=NULL, cv=3, iters=100,
                   kernels=c("linear", "polynomial", "radial", "sigmoid"),
                   cores=2) {

  if(is.null(input.file) | is.null(output.dir)){
    stop("sysSVM ERROR: Input/Output file unknown, please provide input file and output directory")
  }

  trainingMode = TRUE
  if(!is.null(models) & !is.null(scaling.factors)){
    cat("Training mode disabled, using given models to predict cancer genes...", "\n")
    trainingMode = FALSE
  }

  #source("dependecies.R")
  source("R/functions.R")
  source("R/dependencies.R")

  ## Create training and prediction sets
  createDescribeTrainingCGC(input.file=input.file, output.dir=output.dir, exclude.features=exclude.features,
                            trainingMode=trainingMode, models=models, scaling.factors=scaling.factors)


  runNoveltyDetection(output.dir=output.dir, cv=cv, iters=iters, kernels=kernels, cores=cores)

}

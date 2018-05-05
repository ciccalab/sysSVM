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

syssvm <- function(input.file=NULL, output.dir=NULL, exclude.features=NULL, models=NULL, scaling.factors=NULL) {

  if(is.null(input.file) | is.null(output.dir)){
    stop("sysSVM ERROR: Input/Output file unknown, please provide input file and output directory")
  }

  trainingMode = TRUE
  if(!is.null(models) & !is.null(scaling.factors)){
    cat("Training mode disabled, using given models to predict cancer genes...", "\n")
    trainingMode = FALSE
  }

  #source("dependecies.R")
  source("R/prepare.R")
  source("R/dependencies.R")

  createDescribeTrainingCGC(input.file=input.file, output.dir=output.dir, exclude.features=exclude.features,
                            trainingMode=trainingMode, models=models, scaling.factors=scaling.factors)


}

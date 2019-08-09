#' Run sysSVM
#'
#' \code{syssvm} uses e1071 to run an one-class SVM
#'
#'
#' @param input.file Path to input file [default is "example_data/oac_ML_input_formatted.tsv"]
#' @param output.dir Path to output directory where sysSVM results will be written [default is "OAC"]
#' @param exclude.features Features to be excluded from training/prediction [default is c("young", "no_ALL_muts")]
#' @param models Functionality to use already trained models for new samples
#' @param scaling.factors Object with sysSVM scaling factors. Use this to make predictions using already-trained sysSVM models
#' @param cv k-fold cross-validation [default is 3]
#' @param iters Iterations to be performed during ross-validation [default is 1000]
#' @param step Step to summarise results (i.e. every step iterations) [default is 100]
#' @param top.rank Sys-candidates to be considered per patient [default is 10]
#' @param kernels Kernels to be used for training/prediction [default is c("linear", "polynomial", "radial", "sigmoid")]
#' @param cores Multi-core functionality [default is 2]
#' @param ncg.tissue.name NCG tissue for refinement of sys-candidates [default is "esophagus"]
#' @param refine set to TRUE to refine scores based on gtex expression and previously published results
#'
#' @examples syssvm()
#'
#'
#'
#'
#'

syssvm <- function(output.dir="OAC", exclude.features=c("young", "no_NSI_muts"),
                   models=NULL,
                   scaling.factors=NULL,
                   cv=3, iters=1000, step=100, top.rank=10,
                   kernels=c("linear", "polynomial", "radial", "sigmoid"),
                   cores=2, ncg.tissue.name="esophagus", refine=TRUE) {

  trainingMode = TRUE
  if(!is.null(models) & !is.null(scaling.factors)){
    cat("Training mode disabled, using given models to predict cancer genes...", "\n")
    trainingMode = FALSE
  }

  ## Create training and prediction sets
  createDescribeTrainingCGC(output.dir=output.dir, exclude.features=exclude.features,
                            trainingMode=trainingMode, models=models, scaling.factors=scaling.factors)


  runNoveltyDetection(output.dir=output.dir, cv=cv, iters=iters, kernels=kernels, cores=cores)

  scoreGenes(ncg.tissue.name=ncg.tissue.name, output.dir=output.dir,
             iters=iters, step=step, top.rank=top.rank,
             exclude.features=exclude.features, refine=refine)

}

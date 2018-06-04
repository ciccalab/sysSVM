#!/usr/bin/env Rscript

rm(list=ls())

usage = " 
- Script to train SVM one-class classification given a set of hyperparameters
(c, gamma and nu)
- Training is refined using functions in config.R and the functions within
- Result directory is the place where the results will be stored
- Becareful when runnign multiple cancer types to give result directory 
unique fr each cancer type
- Usually the script sun using qsub commands created from run commands function 
in config.R, make the file with the commands executable and submit it to the cluster
- The script also saves the best model based on accuracy trained on the whole training set
- In the optimization function of one-class svm there is only 
"

pkgs = c('plyr','tidyr',
         'caret',
         'ggplot2',
         'grid',
         'gridExtra',
         'dplyr',
         'doMC',
         'xlsx',
         'scales',
         'doParallel',
         'foreach',
         'RMySQL',
         'pROC',
         'e1071',
         'scales',
         'RColorBrewer',
         'sampling')


for (pkgR in pkgs) 
    if (!pkgR %in% rownames(installed.packages())) { 
        cat("Installing: ", pkgR, "\n")
        install.packages( pkgR, repos="http://cran.us.r-project.org")
        library(pkgR, character.only = TRUE, quietly = TRUE)
    } else { 
        library(pkgR, character.only = TRUE, quietly = TRUE)
    }

args = commandArgs(trailingOnly = TRUE)

for(i in 1:length(args)) cat(i, "\t", args[i], '\n')

## This is for when we try different kernels in parallel
## Give the kernel in the name of the job
cancer_type = args[1]
kern = args[2]
#mycost = args[3]
mynu = args[3] 
mygamma = args[4]
mydegree = args[5]
res_dir <- args[6]

## Constants
CV=3
LO=1 ## Leave out i.e here I take CV-1 folds for the train and 1 for the test
ITERS=100

registerDoMC(cores = 3)


#mycost = as.numeric(mycost)
mynu = as.numeric(mynu)
mygamma = as.numeric(mygamma)

create.output.folder = function (res_dir, no) {
    
    if(length(grep("Results",list.dirs(res_dir, recursive = F, full.names=F)))==0){
        dir.create(paste(res_dir, "/Results", sep=""))
    }
    timer    = format(Sys.time(), "%Y-%m-%d.%H_%M_%S")
    if(kern=="linear"){
        job_dir = paste(res_dir, "/Results/", kern,".",mynu, sep="" )
        if(!file.exists(job_dir)){
            dir.create(job_dir)
        }
    }else if (kern=="radial"){
        job_dir = paste(res_dir, "/Results/", kern,".", mynu, ".", mygamma, sep="")
        if(!file.exists(job_dir)){
            dir.create(job_dir)
        }
    }else if (kern=="polynomial"){
        job_dir = paste(res_dir, "/Results/", kern,".", mynu, ".", mygamma, ".", mydegree, sep="")
        if(!file.exists(job_dir)){
            dir.create(job_dir)
        }
    }else if (kern=="sigmoid"){
        job_dir = paste(res_dir, "/Results/", kern,".", mynu, ".", mygamma, sep="")
        if(!file.exists(job_dir)){
            dir.create(job_dir)
        }
    }
    iteration_dir = paste(job_dir, "/iteration",".", no, sep="")
    if(!file.exists(iteration_dir)){
        dir.create(iteration_dir)
    }
    return(iteration_dir)
}

## get the data from path (training)
load(paste(res_dir, "/training_set.Rdata", sep = ""))

## get the negative observations (negatives)
#load(paste(res_dir, "/negative_set.Rdata", sep = ""))
## Then probably we need to restrict the number of negatives according to some criteria
## Because this set is usually very big
## For now just take approximately LO/CV of the training set
#negatives = negatives %>% slice(1:ceiling((nrow(training)*LO)/CV))

#cat(paste(cv, "-fold cross-validation...", sep=""),"\n")
pgenes <- training %>% add_rownames %>% separate(rowname, into=c("Cancer_type", "Sample", "Entrez"), sep="\\.") %>% subset(type=="C") %>% select(Entrez) %>% unique %>% .$Entrez
#fpgenes <- training %>% add_rownames %>% separate(rowname, into=c("Cancer_type", "Sample", "Entrez"), sep="\\.") %>% subset(type=="N") %>% select(Entrez) %>% unique %>% .$Entrez

## create the set of genes (this is for when I do 3x3)
#cvTrainSets = split(pgenes, ceiling(seq_along(pgenes)/(length(pgenes)/cv)))


#for (i in 1:length(cvTrainSets)){
for (i in 1:ITERS){  
  iteration_dir <- create.output.folder(res_dir, i)
    
  ## For when I do 3x3 as I did in the past
  #traingenes <- as.vector(unlist(cvTrainSets[-i]))
  #testgenes <- as.vector(unlist(cvTrainSets[i]))

  traingenes = sample(pgenes,ceiling(((CV-LO)/CV)*length(pgenes)))
  testgenes = pgenes[!(pgenes%in%traingenes)]
  
  ## Extract trainset and testset
  trainset <- training %>% tibble::rownames_to_column() %>% separate(rowname, into=c("Cancer_type", "Sample", "Entrez"), sep="\\.") %>% subset(Entrez %in% traingenes)
  testset <-training %>% tibble::rownames_to_column() %>% separate(rowname, into=c("Cancer_type", "Sample", "Entrez"), sep="\\.") %>% subset(Entrez %in% testgenes)

  ## Make Cancer_type, Sample and Entrez row names
  trainset <- trainset %>% mutate(key=paste(Cancer_type, Sample, Entrez, sep=".")) %>% 
      select(-Cancer_type, -Sample, -Entrez)
  rnames <- trainset$key
  trainset <- trainset %>% select(-key) %>% data.frame()
  row.names(trainset) <- rnames

  testset <- testset %>% mutate(key=paste(Cancer_type, Sample, Entrez, sep=".")) %>% 
      select(-Cancer_type, -Sample, -Entrez)
  rnames <- testset$key
  testset <- testset %>% select(-key) %>% data.frame()
  row.names(testset) <- rnames
    
  ## Exclude columns that are NaN - this should be changed to make the column all 0s
  trainset=Filter(function(x)!all(is.nan(x)), trainset)
  testset=Filter(function(x)!all(is.nan(x)), testset)
  
  ## Store the train and testset
  save(trainset, file=paste(iteration_dir, "/trainset.Rdata", sep=""))
  save(testset, file=paste(iteration_dir, "/testset.Rdata", sep=""))
    
  ## Report to the log
  #cat(paste(no, length(traingenes), paste(traingenes, collapse = ","), length(testgenes), paste(traingenes, collapse = ","), sep = "\t"), "\n")
    
  ## -----------------------------------------------------
  ##                  SVM
  ## -----------------------------------------------------
  cat("Model training...", "\n")
    
  ## e1071 implementation
  if (kern=="linear"){
    svm.model <- svm(type ~ ., data = trainset, kernel=kern, type="one-classification", scale=FALSE, nu=mynu)
  }else if (kern=="radial"){
    svm.model <- svm(type ~ ., data = trainset, kernel=kern, type="one-classification", scale=FALSE, gamma = mygamma, nu=mynu)
  }else if (kern=="polynomial"){
    ## probably changing degree
    svm.model <- svm(type ~ ., data = trainset, kernel=kern, type="one-classification", scale=FALSE, gamma = mygamma, nu=mynu, degree = mydegree)
  }else if (kern=="sigmoid"){
    svm.model <- svm(type ~ ., data = trainset, kernel=kern, type="one-classification", scale=FALSE, gamma = mygamma, nu=mynu)  
  }
  save(svm.model, file=paste(iteration_dir, "/svmModel.Rdata", sep=""))

  cat("Prediction using the trained model...", "\n")
  ## e1071
  prediction <- predict(svm.model, testset%>%select(-type), decision.values = TRUE, probability = TRUE)
  ## save prediction in a file
  save(prediction, file=paste(iteration_dir, "/prediction.Rdata", sep=""))
  rm(svm.model, prediction, trainset, testset) 
}
#sink()

## Assess the instances of the model and train the best model on the whole dataset
## Get to the jobs directory
cat("Estimating performance metrics...", "\n")
modelDir = paste(unlist(strsplit(iteration_dir, "\\/"))[-length(unlist(strsplit(iteration_dir, "\\/")))], collapse="/")
setwd(modelDir)
iter_dirs <- list.dirs(recursive = F)
cv_stats = NULL
for (iter_dir in iter_dirs){
  setwd(iter_dir)
  cv_analysis <- tail(unlist(strsplit(modelDir, "\\/")), n=1)
  iter <- as.numeric(tail(unlist(strsplit(iter_dir, split = "\\.")), n=1))
  load("trainset.Rdata")
  load("testset.Rdata")
  load("svmModel.Rdata")
  load("prediction.Rdata")
  
  noTrain = trainset %>% nrow
  ## ---------------------- Test set -------------------------------------
  ## Concatenate the predictions to the true values
  testset$prediction = prediction
  testset = testset %>% select(type, prediction)
  testset = testset %>% mutate(p=ifelse(prediction==TRUE, "C", "N"))
  testset$type = factor(as.character(testset$type), levels = c("C", "N"))
  testset$p = factor(as.character(testset$p), levels = c("C", "N"))
  
  ## Performance metrics
  ## Confusion matrix
  cM <- confusionMatrix(testset$p,testset$type, positive = c("C"))
  cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, kernel=kern,
                                         nu=mynu, gamma=mygamma, degree=mydegree,
                                         iteration=iter,
                                         type="accuracy", trainSize=noTrain,
                                         class="overall", set = "test",
                                         value=cM$overall[["Accuracy"]]))
  cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, kernel=kern,
                                         nu=mynu, gamma=mygamma, degree=mydegree,
                                         iteration=iter,
                                         type="Sensitivity", trainSize=noTrain,
                                         class="overall", set = "test", 
                                         value=cM$byClass["Sensitivity"]))
  cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, kernel=kern,
                                         nu=mynu, gamma=mygamma, degree=mydegree,
                                         iteration=iter,
                                         type="Specificity", trainSize=noTrain,
                                         class="overall", set = "test", 
                                         value=cM$byClass["Specificity"]))

  ## ---------------------- Train set ------------------------------------- 
  trainset$fitted = svm.model$fitted
  trainset = trainset %>% select(type, fitted)
  trainset = trainset %>% mutate(f=ifelse(fitted==TRUE, "C", "N"))
  trainset$type = factor(as.character(trainset$type), levels = c("C", "N"))
  trainset$f = factor(as.character(trainset$f), levels = c("C", "N"))
    
  ## Confusion matrix
  cM <- confusionMatrix(trainset$f,trainset$type, positive = c("C"))
  cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, kernel=kern,
                                         nu=mynu, gamma=mygamma, degree=mydegree,
                                         iteration=iter,
                                         type="accuracy", trainSize=noTrain,
                                         class="overall", set = "train",
                                         value=cM$overall[["Accuracy"]]))
  cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, kernel=kern,
                                         nu=mynu, gamma=mygamma, degree=mydegree,
                                         iteration=iter,
                                         type="Sensitivity", trainSize=noTrain,
                                         class="overall", set = "train", 
                                         value=cM$byClass["Sensitivity"]))
  cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, kernel=kern,
                                         nu=mynu, gamma=mygamma, degree=mydegree,
                                         iteration=iter,
                                         type="Specificity", trainSize=noTrain,
                                         class="overall", set = "train", 
                                         value=cM$byClass["Specificity"]))
  ## To return to the previous directory
  setwd(modelDir)
  ## And clean
  rm(svm.model, prediction, trainset, testset)
}

## write metrics
write.table(cv_stats, file=paste(modelDir,"/cv_stats.tsv",sep=""), quote=F, row.names = F, sep="\t")









## Prepare training and prediction set
#############
## Functions
#############

## Marking of training/prediction set
markTrainCGC <- function(df){

  geneInfo=sysSVM:::geneInfo
  cancerGenes=sysSVM:::cancerGenes

  ## Get information for cgc/cans and primary sites
  ## Fix gene info table from NCG
  geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
  ## Get a cancer gene with all the associated primary sites and cancer sites
  cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
    group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
                                   cancer_site=paste(unique(cancer_site), collapse=",")) %>%
    ungroup
  geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez")) %>%
    rename(gene_type=cancer_type)

  ## Get false positive genes
  false_positive_genes <- sysSVM:::false_positive_genes
  false_positive_genes <- false_positive_genes %>% select(entrez, symbol)

  ## Get primary site information for the df
  df = df %>% left_join(geneInfo, by=c("entrez"="entrez"))

  ## Training set extraction
  cat("Marking genes for training...", "\n")
  #df = df %>% mutate(TP=ifelse((vogel=="Vog.Oncogene" | vogel=="Vog.TS") & (no_TRUNC_muts!=0 | no_NTDam_muts!=0 | no_GOF_muts!=0 | CNVGain==1 | CNVLoss==1 | ExpT_NET==1), T, F))
  #if( df%>%grepl(tissue, primary_site)%>%nrow==0 ) stop("Tissue not present")


  df = df %>% mutate(TP=ifelse((gene_type=="cgc") & (no_TRUNC_muts!=0 |
                                                       no_NTDam_muts!=0 |
                                                       no_GOF_muts!=0 |
                                                       (CNVGain==1) |
                                                       (Copy_number==0) |
                                                       ((Copy_number==1) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) |
                                                       BND!=0 |
                                                       INS!=0 |
                                                       INV!=0) & !(entrez%in%false_positive_genes$entrez), T, F))


  cat("Adding type annotation for training...", "\n")

  df = df %>% mutate(type=ifelse(TP==T, "C", NA))

  df = df %>% select(-TP)
  ## For cancer-specific training set
  df = df %>% select(-gene_type, -primary_site, -cancer_site)
  df$type = as.factor(df$type)

  ## Exclude genes without any alterations from training & prediction set  (I now added the exclusion of the NCG FPs as an option by default they are not excluded)
  df = df %>% subset((no_TRUNC_muts!=0 |
                        no_NTDam_muts!=0 |
                        no_GOF_muts!=0 |
                        (CNVGain==1) |
                        (Copy_number==0) |
                        ((Copy_number==1) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)) |
                        BND!=0 |
                        INS!=0 |
                        INV!=0) )

  return(df)
}

## Scaling
getScaledTable <- function(df, type="cs", models=NULL, scaling.factors=NULL){ #

  ## is.nan for data frames
  is.nan.data.frame <- function(x){do.call(cbind, lapply(x, is.nan))}

  cat("Feature scaling: All factors are excluded from scaling...\n")

  factor_cols <- names(sapply(df, function(x) class(x))[sapply(df, function(x) class(x))=="factor"])

  df_f <- df[,names(df)%in%factor_cols]

  df <-df[,!(names(df)%in%factor_cols)]

  if(is.null(scaling.factors)){
    ## Get mean and sd for scaling new observations
    scaling.factors = sapply(df, function(cl) list(means=mean(cl,na.rm=TRUE), sds=sd(cl,na.rm=TRUE)))
    cat("Scaling...", "\n")
    df <- data.frame(lapply(df, function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))
    df <- cbind(df, df_f)
    #df <- data.frame(lapply(df, function(x) rescale(x, to=c(-1,1)))) ## only cs scaling for now
  }else{
    load(scaling.factors)
    cat("Scaling...", "\n")
    if (type=="cs"){
      df <- data.frame(lapply(seq_along(df), function(y,n,i) {
        (y[i]-unlist(mean_sd["means",n[i]]))/unlist(mean_sd["sds",n[i]])
      }, y=df, n=names(df)))
      df <- cbind(df, df_f)
    }
  }

  return(list(df_scaled=df, scaling.factors=scaling.factors))
}


## Create training and prediction sets
createDescribeTrainingCGC = function(output.dir=NULL, exclude.features=NULL, trainingMode=TRUE, models=NULL, scaling.factors=NULL){

  require(dplyr)
  require(tidyr)

  ## Get data
  d = sysSVM:::oac_data

  FEATURES = colnames(d)[!colnames(d)%in%exclude.features]

  d = d %>% select(one_of(c("sample", "entrez", FEATURES))) %>% data.frame()
  ## We need to make it factors
  cat("Converting binary features to factors", "\n")

  ## Get which features have only two levels and onvert them to factors
  to.be.factors = lapply(d, function(x){
    length(unique(x))==2
  })

  fcols=colnames(d)[unlist(to.be.factors)]
  fcols = fcols[!fcols%in%c("sample", "entrez")]

  cols <- which(colnames(d) %in% fcols)
  for(i in cols){
    d[,i] = factor(d[,i], levels = c(0,1))
  }


  ## Mark training observations before you scale
  if(trainingMode){
    cat("Creating training and prediction set", "\n")

    if (!file.exists(output.dir)){
      dir.create(file.path(output.dir))
    }else{
      stop("sysSVM ERROR: Output directory exists. Can't overwrite...")
    }

    d = markTrainCGC(df=d)

    rnames = paste0(d$sample, ".", d$entrez)
    d=d %>% select(-sample, -entrez) %>% data.frame()
    rownames(d) = rnames
    ## Save training and prediction set without scaling
    training_ns = d %>% subset(!is.na(type))
    save(training_ns, file=paste(output.dir,"/training_set_noScale.Rdata", sep=""))

    prediction_ns = d %>% subset(is.na(type)) %>% select(-type)
    save(prediction_ns, file=paste(output.dir,"/prediction_set_noScale.Rdata", sep=""))


    ## Scale dataset
    d = getScaledTable(df = d, scaling.factors = scaling.factors)
    scaling.factors = d[["scaling.factors"]]
    ## Save means and sds
    if(!is.null(scaling.factors)){
      save(scaling.factors, file=paste(output.dir,"/scaling.factors.Rdata", sep=""))
    }

    d = d[["df_scaled"]]
    training = d %>% subset(!is.na(type))
    save(training, file=paste(output.dir,"/training_set.Rdata", sep=""))

    prediction = d %>% subset(is.na(type)) %>% select(-type)
    save(prediction, file=paste(output.dir,"/prediction_set.Rdata", sep=""))

  }else{
    ## Scale dataset
    if(is.null(scaling.factors)){
      stop("sysSVM ERROR: Please provide scaling values or train sysSVM", "\n")
    }
    cat("Scaling data and preparing table for prediction...", "\n")

    if (!file.exists(output.dir)){
      dir.create(file.path(output.dir))
    }else{
      stop("sysSVM ERROR: Output directory exists. Can't overwrite...")
    }

    save(d, file=paste(output.dir,"/prediction_set_noScale.Rdata", sep=""))
    d = getScaledTable(df = d, scaling.factors = scaling.factors)

    d = d[["df_scaled"]]
    save(d, file=paste(output.dir,"/prediction_set.Rdata", sep=""))
  }

}


## Run novelty detection
runNoveltyDetection = function(output.dir=NULL, cv=3, iters=100,
                               kernels=c("linear", "polynomial", "radial", "sigmoid"),
                               cores=2){

  require(doSNOW)
  require(parallel)

  ## Constants
  NU = seq(0.05, 0.9, 0.05)
  GAMMA = 2^seq(-7,4)
  DEGREE = c(3, 4, 9)

  LO=1 ## Leave out i.e here I take CV-1 folds for the train and 1 for the test
  CV=cv

  ## Start with linear
  cat("Grid search for:", "\n")

  create.output.folder = function (output.dir, kern, mynu, mygamma, mydegree, no) {

    if(length(grep("CV",list.dirs(output.dir, recursive = F, full.names=F)))==0){
      dir.create(paste(output.dir, "/CV", sep=""))
    }
    timer    = format(Sys.time(), "%Y-%m-%d.%H_%M_%S")
    if(kern=="linear"){
      job_dir = paste(output.dir, "/CV/", kern,".",mynu, sep="" )
      if(!file.exists(job_dir)){
        dir.create(job_dir)
      }
    }else if (kern=="radial"){
      job_dir = paste(output.dir, "/CV/", kern,".", mynu, ".", mygamma, sep="")
      if(!file.exists(job_dir)){
        dir.create(job_dir)
      }
    }else if (kern=="polynomial"){
      job_dir = paste(output.dir, "/CV/", kern,".", mynu, ".", mygamma, ".", mydegree, sep="")
      if(!file.exists(job_dir)){
        dir.create(job_dir)
      }
    }else if (kern=="sigmoid"){
      job_dir = paste(output.dir, "/CV/", kern,".", mynu, ".", mygamma, sep="")
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
  load(paste(output.dir, "/training_set.Rdata", sep = ""))
  pgenes <- training %>% add_rownames %>% separate(rowname, into=c("sample", "entrez"), sep="\\.") %>% subset(type=="C") %>% select(entrez) %>% unique %>% .$entrez

  cv_stats_full = NULL
  for(k in kernels){
    cat(k, "\n")
    if(k=="linear"){
      total = length(NU)

      cl <- makeCluster(cores)
      registerDoSNOW(cl)
      # create progress bar
      pb <- txtProgressBar(max = total, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      result <- foreach(i = 1:total, .combine = rbind, .packages=c('caret', 'sampling', 'e1071', 'tidyr', 'dplyr'),
                        .options.snow = opts) %dopar%
                        {

                            Sys.sleep(0.1)
                            mynu=NU[i]
                            mygamma=0
                            mydegree=0

                            cv_stats_grid = NULL
                            for(y in 1:iters){
                              iteration_dir <- create.output.folder(output.dir, k, mynu, mygamma, mydegree, y)
                              print(iteration_dir)
                              traingenes = sample(pgenes,ceiling(((CV-LO)/CV)*length(pgenes)))
                              testgenes = pgenes[!(pgenes%in%traingenes)]

                              ## Extract trainset and testset
                              trainset <- training %>% tibble::rownames_to_column() %>% separate(rowname, into=c("sample", "entrez"), sep="\\.") %>% subset(entrez %in% traingenes)
                              testset <-training %>% tibble::rownames_to_column() %>% separate(rowname, into=c("sample", "entrez"), sep="\\.") %>% subset(entrez %in% testgenes)

                              ## Make Cancer_type, Sample and Entrez row names
                              trainset <- trainset %>% mutate(key=paste(sample, entrez, sep=".")) %>%
                                select(-sample, -entrez)
                              rnames <- trainset$key
                              trainset <- trainset %>% select(-key) %>% data.frame()
                              row.names(trainset) <- rnames

                              testset <- testset %>% mutate(key=paste(sample, entrez, sep=".")) %>%
                                select(-sample, -entrez)
                              rnames <- testset$key
                              testset <- testset %>% select(-key) %>% data.frame()
                              row.names(testset) <- rnames

                              ## Exclude columns that are NaN - this may be changed in later versions to make the column all 0s
                              trainset=Filter(function(x)!all(is.nan(x)), trainset)
                              testset=Filter(function(x)!all(is.nan(x)), testset)

                              ## Store the train and testset
                              save(trainset, file=paste(iteration_dir, "/trainset.Rdata", sep=""))
                              save(testset, file=paste(iteration_dir, "/testset.Rdata", sep=""))

                              ## -----------------------------------------------------
                              ##                  SVM
                              ## -----------------------------------------------------

                              svm.model <- svm(type ~ ., data = trainset, kernel=k, type="one-classification", scale=FALSE, nu=mynu)
                              save(svm.model, file=paste(iteration_dir, "/svmModel.Rdata", sep=""))

                              prediction <- predict(svm.model, testset%>%select(-type), decision.values = TRUE, probability = TRUE)
                              ## save prediction in a file
                              save(prediction, file=paste(iteration_dir, "/prediction.Rdata", sep=""))

                              ## Get the performance
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
                              cv_stats = NULL
                              cM <- confusionMatrix(testset$p,testset$type, positive = c("C"))
                              cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                     nu=mynu, gamma=mygamma, degree=mydegree,
                                                                     iteration=y,
                                                                     type="accuracy", trainSize=noTrain,
                                                                     class="overall", set = "test",
                                                                     value=cM$overall[["Accuracy"]]))
                              cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                     nu=mynu, gamma=mygamma, degree=mydegree,
                                                                     iteration=y,
                                                                     type="Sensitivity", trainSize=noTrain,
                                                                     class="overall", set = "test",
                                                                     value=cM$byClass["Sensitivity"]))
                              cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                     nu=mynu, gamma=mygamma, degree=mydegree,
                                                                     iteration=y,
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
                              cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                     nu=mynu, gamma=mygamma, degree=mydegree,
                                                                     iteration=y,
                                                                     type="accuracy", trainSize=noTrain,
                                                                     class="overall", set = "train",
                                                                     value=cM$overall[["Accuracy"]]))
                              cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                     nu=mynu, gamma=mygamma, degree=mydegree,
                                                                     iteration=y,
                                                                     type="Sensitivity", trainSize=noTrain,
                                                                     class="overall", set = "train",
                                                                     value=cM$byClass["Sensitivity"]))
                              cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                     nu=mynu, gamma=mygamma, degree=mydegree,
                                                                     iteration=y,
                                                                     type="Specificity", trainSize=noTrain,
                                                                     class="overall", set = "train",
                                                                     value=cM$byClass["Specificity"]))
                              cv_stats_grid = rbind(cv_stats_grid, cv_stats)
                            }
                          return(cv_stats_grid)
                        }
      close(pb)
      stopCluster(cl)
      cv_stats_full = rbind(cv_stats_full, result)
    }

    if(k=="polynomial"){
      param_grid=expand.grid(NU, GAMMA, DEGREE)
      total = nrow(param_grid)

      cl <- makeCluster(cores)
      registerDoSNOW(cl)
      # create progress bar
      pb <- txtProgressBar(max = total, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      result <- foreach(i = 1:total, .combine = rbind, .packages=c('caret', 'sampling', 'e1071', 'tidyr', 'dplyr'),
                        .options.snow = opts) %dopar%
                        {

                          Sys.sleep(0.1)
                          mynu=param_grid[i, 1]
                          mygamma=param_grid[i, 2]
                          mydegree=param_grid[i, 3]

                          cv_stats_grid = NULL
                          for(y in 1:iters){
                            iteration_dir <- create.output.folder(output.dir, k, mynu, mygamma, mydegree, y)

                            traingenes = sample(pgenes,ceiling(((CV-LO)/CV)*length(pgenes)))
                            testgenes = pgenes[!(pgenes%in%traingenes)]

                            ## Extract trainset and testset
                            trainset <- training %>% tibble::rownames_to_column() %>% separate(rowname, into=c("sample", "entrez"), sep="\\.") %>% subset(entrez %in% traingenes)
                            testset <-training %>% tibble::rownames_to_column() %>% separate(rowname, into=c("sample", "entrez"), sep="\\.") %>% subset(entrez %in% testgenes)

                            ## Make Cancer_type, Sample and Entrez row names
                            trainset <- trainset %>% mutate(key=paste(sample, entrez, sep=".")) %>%
                              select(-sample, -entrez)
                            rnames <- trainset$key
                            trainset <- trainset %>% select(-key) %>% data.frame()
                            row.names(trainset) <- rnames

                            testset <- testset %>% mutate(key=paste(sample, entrez, sep=".")) %>%
                              select(-sample, -entrez)
                            rnames <- testset$key
                            testset <- testset %>% select(-key) %>% data.frame()
                            row.names(testset) <- rnames

                            ## Exclude columns that are NaN - this may be changed in later versions to make the column all 0s
                            trainset=Filter(function(x)!all(is.nan(x)), trainset)
                            testset=Filter(function(x)!all(is.nan(x)), testset)

                            ## Store the train and testset
                            save(trainset, file=paste(iteration_dir, "/trainset.Rdata", sep=""))
                            save(testset, file=paste(iteration_dir, "/testset.Rdata", sep=""))

                            ## -----------------------------------------------------
                            ##                  SVM
                            ## -----------------------------------------------------

                            svm.model <- svm(type ~ ., data = trainset, kernel=k, type="one-classification", scale=FALSE, gamma = mygamma, nu=mynu, degree = mydegree)
                            save(svm.model, file=paste(iteration_dir, "/svmModel.Rdata", sep=""))

                            prediction <- predict(svm.model, testset%>%select(-type), decision.values = TRUE, probability = TRUE)
                            ## save prediction in a file
                            save(prediction, file=paste(iteration_dir, "/prediction.Rdata", sep=""))

                            ## Get the performance
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
                            cv_stats = NULL
                            cM <- confusionMatrix(testset$p,testset$type, positive = c("C"))
                            cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                   nu=mynu, gamma=mygamma, degree=mydegree,
                                                                   iteration=y,
                                                                   type="accuracy", trainSize=noTrain,
                                                                   class="overall", set = "test",
                                                                   value=cM$overall[["Accuracy"]]))
                            cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                   nu=mynu, gamma=mygamma, degree=mydegree,
                                                                   iteration=y,
                                                                   type="Sensitivity", trainSize=noTrain,
                                                                   class="overall", set = "test",
                                                                   value=cM$byClass["Sensitivity"]))
                            cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                   nu=mynu, gamma=mygamma, degree=mydegree,
                                                                   iteration=y,
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
                            cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                   nu=mynu, gamma=mygamma, degree=mydegree,
                                                                   iteration=y,
                                                                   type="accuracy", trainSize=noTrain,
                                                                   class="overall", set = "train",
                                                                   value=cM$overall[["Accuracy"]]))
                            cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                   nu=mynu, gamma=mygamma, degree=mydegree,
                                                                   iteration=y,
                                                                   type="Sensitivity", trainSize=noTrain,
                                                                   class="overall", set = "train",
                                                                   value=cM$byClass["Sensitivity"]))
                            cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                   nu=mynu, gamma=mygamma, degree=mydegree,
                                                                   iteration=y,
                                                                   type="Specificity", trainSize=noTrain,
                                                                   class="overall", set = "train",
                                                                   value=cM$byClass["Specificity"]))
                            cv_stats_grid = rbind(cv_stats_grid, cv_stats)
                          }
                          return(cv_stats_grid)
                        }
      close(pb)
      stopCluster(cl)
      cv_stats_full = rbind(cv_stats_full, result)
    }

    if(k=="radial" | k=="sigmoid"){
      param_grid=expand.grid(NU, GAMMA)
      total = nrow(param_grid)

      cl <- makeCluster(cores)
      registerDoSNOW(cl)
      # create progress bar
      pb <- txtProgressBar(max = total, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      result <- foreach(i = 1:total, .combine = rbind, .packages=c('caret', 'sampling', 'e1071', 'tidyr', 'dplyr'),
                        .options.snow = opts) %dopar%
                        {

                          Sys.sleep(0.1)
                          mynu=param_grid[i, 1]
                          mygamma=param_grid[i, 2]
                          mydegree=0
                          cv_stats_grid = NULL
                          for(y in 1:iters){
                            iteration_dir <- create.output.folder(output.dir, k, mynu, mygamma, mydegree, y)

                            traingenes = sample(pgenes,ceiling(((CV-LO)/CV)*length(pgenes)))
                            testgenes = pgenes[!(pgenes%in%traingenes)]

                            ## Extract trainset and testset
                            trainset <- training %>% tibble::rownames_to_column() %>% separate(rowname, into=c("sample", "entrez"), sep="\\.") %>% subset(entrez %in% traingenes)
                            testset <-training %>% tibble::rownames_to_column() %>% separate(rowname, into=c("sample", "entrez"), sep="\\.") %>% subset(entrez %in% testgenes)

                            ## Make Cancer_type, Sample and Entrez row names
                            trainset <- trainset %>% mutate(key=paste(sample, entrez, sep=".")) %>%
                              select(-sample, -entrez)
                            rnames <- trainset$key
                            trainset <- trainset %>% select(-key) %>% data.frame()
                            row.names(trainset) <- rnames

                            testset <- testset %>% mutate(key=paste(sample, entrez, sep=".")) %>%
                              select(-sample, -entrez)
                            rnames <- testset$key
                            testset <- testset %>% select(-key) %>% data.frame()
                            row.names(testset) <- rnames

                            ## Exclude columns that are NaN - this may be changed in later versions to make the column all 0s
                            trainset=Filter(function(x)!all(is.nan(x)), trainset)
                            testset=Filter(function(x)!all(is.nan(x)), testset)

                            ## Store the train and testset
                            save(trainset, file=paste(iteration_dir, "/trainset.Rdata", sep=""))
                            save(testset, file=paste(iteration_dir, "/testset.Rdata", sep=""))

                            ## -----------------------------------------------------
                            ##                  SVM
                            ## -----------------------------------------------------

                            svm.model <- svm(type ~ ., data = trainset, kernel=k, type="one-classification", scale=FALSE, gamma = mygamma, nu=mynu)
                            save(svm.model, file=paste(iteration_dir, "/svmModel.Rdata", sep=""))

                            prediction <- predict(svm.model, testset%>%select(-type), decision.values = TRUE, probability = TRUE)
                            ## save prediction in a file
                            save(prediction, file=paste(iteration_dir, "/prediction.Rdata", sep=""))

                            ## Get the performance
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
                            cv_stats = NULL
                            cM <- confusionMatrix(testset$p,testset$type, positive = c("C"))
                            cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                   nu=mynu, gamma=mygamma, degree=mydegree,
                                                                   iteration=y,
                                                                   type="accuracy", trainSize=noTrain,
                                                                   class="overall", set = "test",
                                                                   value=cM$overall[["Accuracy"]]))
                            cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                   nu=mynu, gamma=mygamma, degree=mydegree,
                                                                   iteration=y,
                                                                   type="Sensitivity", trainSize=noTrain,
                                                                   class="overall", set = "test",
                                                                   value=cM$byClass["Sensitivity"]))
                            cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                   nu=mynu, gamma=mygamma, degree=mydegree,
                                                                   iteration=y,
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
                            cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                   nu=mynu, gamma=mygamma, degree=mydegree,
                                                                   iteration=y,
                                                                   type="accuracy", trainSize=noTrain,
                                                                   class="overall", set = "train",
                                                                   value=cM$overall[["Accuracy"]]))
                            cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                   nu=mynu, gamma=mygamma, degree=mydegree,
                                                                   iteration=y,
                                                                   type="Sensitivity", trainSize=noTrain,
                                                                   class="overall", set = "train",
                                                                   value=cM$byClass["Sensitivity"]))
                            cv_stats <- rbind(cv_stats, data.frame(kernel=k,
                                                                   nu=mynu, gamma=mygamma, degree=mydegree,
                                                                   iteration=y,
                                                                   type="Specificity", trainSize=noTrain,
                                                                   class="overall", set = "train",
                                                                   value=cM$byClass["Specificity"]))
                            cv_stats_grid = rbind(cv_stats_grid, cv_stats)
                          }
                          return(cv_stats_grid)
                        }
      close(pb)
      stopCluster(cl)
      cv_stats_full = rbind(cv_stats_full, result)
    }

  }
  ## write metrics
  write.table(cv_stats_full, file=paste(output.dir,"/cv_stats.tsv",sep=""), quote=F, row.names = F, sep="\t")
}


## Calculate the score of each gene
## This function will run every 100 iters and calculate the most frequent set of genes
getScore = function(path=NULL,
                    linear_preds=NULL, linearCVS=NULL, linearVar=NULL,
                    radial_preds=NULL, radialCVS=NULL, radialVar=NULL,
                    polynomial_preds=NULL, polynomialCVS=NULL, polynomialVar=NULL,
                    sigmoid_preds=NULL, sigmoidCVS=NULL, sigmoidVar=NULL){


  ## Get the number of predictions before intersection
  all_preds = NULL
  all_scores = NULL

  if(!is.null(linear_preds) & !is.null(linearCVS)){
    load(paste0(path,linear_preds))
    linear = preds
    rm(preds)
    ## Push it to the all_preds
    all_preds = rbind(all_preds, linear %>% subset(label) %>% mutate(kernel="linear", training="cgc") %>% data.frame())

    ## Kernel-specific score
    linear = linear %>% arrange(sample, dv) %>% group_by(sample) %>% mutate(rank=row_number(-dv)) %>% ungroup
    linear = linear %>% group_by(sample) %>% mutate(N=max(rank), Nlog10=log10(N)) %>% ungroup
    linear = linear %>% mutate(cvs=linearCVS, var=linearVar) %>% mutate(ratio=cvs/var)
    linear = linear %>% mutate(score=-log10(rank/N)*cvs)
    linear = linear %>% mutate(kernel="linear") %>% data.frame()
    all_scores = rbind(all_scores, linear)
  }

  if(!is.null(radial_preds) & !is.null(radialCVS)){
    load(paste0(path,radial_preds))
    radial = preds
    rm(preds)
    ## Push it to the all_preds
    all_preds = rbind(all_preds, radial %>% subset(label) %>% mutate(kernel="radial", training="cgc") %>% data.frame())

    ## Kernel-specific score
    radial = radial %>% arrange(sample, dv) %>% group_by(sample) %>% mutate(rank=row_number(-dv)) %>% ungroup
    radial = radial %>% group_by(sample) %>% mutate(N=max(rank), Nlog10=log10(N)) %>% ungroup
    radial = radial %>% mutate(cvs=radialCVS, var=radialVar) %>% mutate(ratio=cvs/var)
    radial = radial %>% mutate(score=-log10(rank/N)*cvs)
    radial = radial %>% mutate(kernel="radial") %>% data.frame()
    all_scores = rbind(all_scores, radial)
  }

  if(!is.null(polynomial_preds) & !is.null(polynomialCVS)){
    load(paste0(path,polynomial_preds))
    polynomial = preds
    rm(preds)
    ## Push it to the all_preds
    all_preds = rbind(all_preds, polynomial %>% subset(label) %>% mutate(kernel="polynomial", training="cgc") %>% data.frame())

    ## Kernel-specific score
    polynomial = polynomial %>% arrange(sample, dv) %>% group_by(sample) %>% mutate(rank=row_number(-dv)) %>% ungroup
    polynomial = polynomial %>% group_by(sample) %>% mutate(N=max(rank), Nlog10=log10(N)) %>% ungroup
    polynomial = polynomial %>% mutate(cvs=polynomialCVS, var=polynomialVar) %>% mutate(ratio=cvs/var)
    polynomial = polynomial %>% mutate(score=-log10(rank/N)*cvs)
    polynomial = polynomial %>% mutate(kernel="polynomial") %>% data.frame()
    all_scores = rbind(all_scores, polynomial)
  }

  if(!is.null(sigmoid_preds) & !is.null(sigmoidCVS)){
    load(paste0(path,sigmoid_preds))
    sigmoid = preds
    rm(preds)
    ## Push it to the all_preds
    all_preds = rbind(all_preds, sigmoid %>% subset(label) %>% mutate(kernel="sigmoid", training="cgc") %>% data.frame())

    ## Kernel-specific score
    sigmoid = sigmoid %>% arrange(sample, dv) %>% group_by(sample) %>% mutate(rank=row_number(-dv)) %>% ungroup
    sigmoid = sigmoid %>% group_by(sample) %>% mutate(N=max(rank), Nlog10=log10(N)) %>% ungroup
    sigmoid = sigmoid %>% mutate(cvs=sigmoidCVS, var=sigmoidVar) %>% mutate(ratio=cvs/var)
    sigmoid = sigmoid %>% mutate(score=-log10(rank/N)*cvs)
    sigmoid = sigmoid %>% mutate(kernel="sigmoid") %>% data.frame()
    all_scores = rbind(all_scores, sigmoid)
  }

  ## Final score for each gene across kernels
  scores = all_scores %>% group_by(sample, entrez, gene_type, N, Nlog10) %>%
    summarise(sum_of_scores=sum(score), kernels_used=length(unique(kernel))) %>%
    mutate(score=(1/(kernels_used*Nlog10))*sum_of_scores) %>% ungroup %>% data.frame()


  ## Add from how many kernels they have predicted as TRUE
  preds2kernels = all_preds %>% subset(label) %>% group_by(sample, entrez) %>% summarise(kernels_predicted=paste(unique(kernel), collapse=","), kernels_predicted_no=length(unique(kernel))) %>% ungroup
  scores = scores %>% left_join(preds2kernels)
  scores$kernels_predicted_no[is.na(scores$kernels_predicted_no)] = 0

  return(list(all_scores=all_scores, scores=scores))

}

## Platt scaling
plattScaling = function(model, predictions){
  ## Impementation was taken from pseudocode in:
  ## Probabilistic Outputs for Support Vector Machines and Comparisons to Regularized Likelihood Methods, Platt, 1999
  ## Construct the training set for Platt's sigmoid

  decision_values = model$decision.values %>% as.vector()
  labels = model$fitted %>% as.vector()
  df = data.frame(dv=decision_values, lbs = labels)

  ## the prediction at the moment is a predict object coming from svm predict
  dv = attr(predictions, "decision.values") %>% as.vector()
  preds = data.frame(dv = dv)
  preds$label = as.vector(predictions)
  rownames(preds) = names(predictions)

  ## Input parameters
  out = df$dv
  target = df$lbs
  prior1 = df%>%subset(lbs==TRUE)%>%nrow
  prior0 = df%>%subset(lbs==FALSE)%>%nrow
  ## Output
  ## A, B = parameters of sigmoid

  A = 0
  B = log((prior0+1)/(prior1+1))
  hiTarget = (prior1+1)/(prior1+2)
  loTarget = 1/(prior0+2)
  lambda = 1e-3
  olderr = 1e300


  pp = rep((prior1+1)/(prior0+prior1+2), nrow(df))
  count = 0

  for (it in 1:100){
    a=0
    b=0
    c=0
    d=0
    e=0

    ## First compute the Hessian & gradient of error function
    ## with respect to A & B
    for (i in 1:length(pp)){
      if (target[i]==TRUE){
        t = hiTarget
      }else{
        t = loTarget
      }

      d1 = pp[i]-t
      d2 = pp[i]*(1-pp[i])
      a = a + out[i]*out[i]*d2
      b = b + d2
      c = c + out[i]*d2
      d = d + out[i]*d1
      e = e + d1
    }
    ## If gradient is really tiny, then stop
    if ((abs(d) < 1e-9) & (abs(e) < 1e-9)){break}

    oldA = A
    oldB = B
    err = 0

    ## Loop until goodness of fit increases
    while(TRUE){
      det = (a+lambda)*(b+lambda)-c*c
      if (det==0){ ## If determinant if Hessian is 0
        ## Increase stabilizer
        lambda = lambda*10
        next
      }
      A = oldA + ((b+lambda)*d-c*e)/det
      B = oldB + ((a+lambda)*e-c*d)/det

      ## Now compute the goodness of fit
      err = 0
      for (i in 1:length(pp)){
        p = 1/(1+exp(out[i]*A+B))
        pp[i] = p
        ## At this step, make sure log(0) returns -200
        if (p<=1.383897e-87){
          err = err - t*(-200)+(1-t)*log(1-p)
        }else if (p==1){
          err = err - t*log(p)+(1-t)*(-200)
        }else{
          err = err - t*log(p)+(1-t)*log(1-p)
        }
        if(err==-Inf) browser()
      }
      if (err < olderr*(1+1e-7)){
        lambda = lambda*0.1
        break
      }
      ## Error did not decrease: increase stabilizer by factor of 10
      ## and try again
      lambda = lambda* 10
      if (lambda >= 1e6){ ## something is broken. Give up
        break
      }
    }
    diff = err-olderr
    scale = 0.5*(err+olderr+1)
    if (diff > -1e-3*scale & diff < 1e-7*scale){
      count = count + 1
    }else{
      count = 0
    }
    olderr = err
    if (count==3){
      break
    }


  }

  applyPlatt = function(x){
    p = 1/(1+exp(x*A+B))
  }

  preds$prob = apply(preds, 1, function(x) applyPlatt(x[1]))
  return(preds)
}

## Reproduce OAC results
getSysCansOAC = function(output.dir=NULL, working.dir = NULL, gtex.tissue=NULL,
                         scores_fn="scores.Rdata",
                         exclude_expr=c("Not Expressed")){

  if(is.null(gtex.tissue)){
    stop("Please provide tissue")
  }

  ## Get the data we need for the definition of the sys cans
  ## Load geneInfo to get the symbols as well
  geneInfo=sysSVM:::geneInfo
  geneInfoNCG5 = geneInfo %>% select(entrez, symbol) %>% unique

  ## Load the CGC that we will consider drivers
  cat("Getting CGC annotation", "\n")
  cgcs = sysSVM:::cgcs
  ## refine to those that will be retained
  retained_CGC = cgcs %>% subset(RET.DIS=="RET")
  retained_CGC_trans = retained_CGC %>% subset(grepl("Trans", Drivers.to.be.retained))
  ## Get training and validation set
  load(paste0(output.dir, "/training_set_noScale.Rdata"))
  training_ns = training_ns%>%tibble::rownames_to_column()%>%
    separate(rowname, into=c("sample", "entrez"), sep="\\.") %>%
    mutate(entrez=as.numeric(entrez)) %>% data.frame()
  ## Refine the training set to include only the ones we decided to coonsider drivers
  cat("Refining CGC drivers", "\n")
  cgc_drivers = NULL
  for (i in 1:nrow(retained_CGC)){
    en = retained_CGC$entrez[i]
    gene_type= retained_CGC$GENETICS[i]
    symbol = retained_CGC$symbol[i]
    retained_drivers = retained_CGC$Drivers.to.be.retained[i]
    #message(paste(symbol, gene_type, sep = "\t"))
    if (gene_type=="O"){
      if(grepl("Trans", retained_drivers)){
        dd = training_ns %>% subset(entrez==en) %>% subset(CNVGain==1 | no_NTDam_muts!=0 | no_GOF_muts!=0 | BND==1 | INV==1 | INS==1)
      }else{
        dd = training_ns %>% subset(entrez==en) %>% subset(CNVGain==1 | no_NTDam_muts!=0 | no_GOF_muts!=0)
      }
      cgc_drivers = rbind(cgc_drivers, dd)
    }else if (gene_type=="TS"){
      if(grepl("Trans", retained_drivers)){
        dd = training_ns %>% subset(entrez==en) %>% subset(Copy_number==0 | no_TRUNC_muts!=0 | no_NTDam_muts!=0 | (Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)) | BND==1 | INV==1 | INS==1)
      }else{
        dd = training_ns %>% subset(entrez==en) %>% subset(Copy_number==0 | no_TRUNC_muts!=0 | no_NTDam_muts!=0 | (Copy_number==1 & (no_TRUNC_muts>0 | no_NTDam_muts>0)))
      }
      cgc_drivers = rbind(cgc_drivers, dd)
    }else if (gene_type=="O/TS"){
      if(grepl("Trans", retained_drivers)){
        ## All drivers are considered
        dd = training_ns %>% subset(entrez==en)
      }else{
        ## Here we need to exlude all the genes with only SVs
        dd = training_ns %>% subset(entrez==en) %>% subset(!(no_TRUNC_muts==0 & no_NTDam_muts==0 & no_GOF_muts==0 & Copy_number!=0 & CNVGain==0 & (BND!=0 | INS!=0 | INV!=0)))
      }
      cgc_drivers = rbind(cgc_drivers, dd)
    }else{
      stop("Unknown gene type")
    }
  }

  ## Now add in the CGC drivers the CGCs identified from previous studies
  previous = sysSVM:::previous
  ## Get only those that were not considered already
  previous = previous %>% subset(to.be.considered & !already.considered)
  ## Get the entrez for these genes
  previous = previous %>% left_join(geneInfo%>%select(entrez, symbol), by=c("Gene.symbol"="symbol"))

  for(i in 1:nrow(previous)){
    en = previous$entrez[i]
    gene_type= previous$Type[i]
    symbol = previous$Gene.symbol[i]

    if(gene_type=="SVs"){
      dd = training_ns %>% subset(entrez==en) %>% subset(BND==1 | INV==1 | INS==1 | CNVLoss==1 | CNVGain==1)
    }else if(gene_type=="CNVs"){
      dd = training_ns %>% subset(entrez==en) %>% subset(CNVLoss==1 | CNVGain==1)
    }
    cgc_drivers = rbind(cgc_drivers, dd)
  }


  load(paste0(output.dir, "/prediction_set_noScale.Rdata"))
  prediction_ns = prediction_ns%>%tibble::rownames_to_column()%>%
    separate(rowname, into=c("sample", "entrez"), sep="\\.")%>%
    mutate(type="P")%>%mutate(entrez=as.numeric(entrez))%>%data.frame()
  cohort = rbind(cgc_drivers, prediction_ns) ## Instead of the whole training, we use only the subset of the refined CGC
  ## Get scores
  load(paste0(working.dir, "/",scores_fn))
  scores = scores[["scores"]]

  cohort = cohort %>% left_join(scores)
  cohort$gene_type[is.na(cohort$gene_type)] = "cgc"
  cohort = cohort %>% left_join(geneInfoNCG5)

  samples = unique(cohort$sample)
  sys_cans = NULL
  if (!is.null(gtex.tissue)){
    cat("Expression is considered for the selection of sys cans", "\n")
    ## Load expression data (GTeX)
    cat("Loading GTeX expression data", "\n")
    gtex=sysSVM:::gtex
    gtex = gtex %>% select(entrez, tissue, exp_level) %>% spread(tissue, exp_level)
    ## Checkpoint
    if(!(gtex.tissue%in%colnames(gtex))){
      stop("Wrong tissue name. Please provide a valid GTeX tissue name.")
    }
    ## Select tissue of the cancer type
    gtex = gtex[,c("entrez", gtex.tissue)]
    colnames(gtex) = c("entrez", "gtex_tissue_expression")
    ## bring in the cohort table the expression
    cohort = cohort %>% left_join(gtex)


    total = length(samples)
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    for (i in 1:length(samples)){
      Sys.sleep(0.1)
      ## When you get the candidates here you need to take into account expression in normal tissues and in the cancer sample

      ## Here I get the ranks for the non-CGC genes. The CGCs will be added as patient.rank==NA
      ## I also add 2 ranks, one with all the genes and one without amplifications
      d = cohort %>% subset(sample==samples[i] & gene_type!="cgc" & !gtex_tissue_expression%in%exclude_expr & !is.na(gtex_tissue_expression)) %>%
        arrange(desc(score)) %>% mutate(patient.rank=row_number(-score))
      d_withoutAmp = d %>% subset(sample==samples[i] & gene_type!="cgc" & CNVGain==0 & !gtex_tissue_expression%in%exclude_expr & !is.na(gtex_tissue_expression)) %>%
        arrange(desc(score)) %>% mutate(patient.rank.no.amp=row_number(-score))
      ## Bring in the rank without amplification
      d = suppressMessages(d %>% left_join(d_withoutAmp))

      ## Add the CGCs
      d_cgcs = cohort %>% subset(sample==samples[i] & gene_type=="cgc" & !gtex_tissue_expression%in%exclude_expr & !is.na(gtex_tissue_expression))
      d = plyr::rbind.fill(d, d_cgcs)

      ## Add them to the big table
      sys_cans = rbind(sys_cans, d)
      setTxtProgressBar(pb, i)
    }
  }

  ## Add false positive annotation
  false_positive_genes <- sysSVM:::false_positive_genes
  false_positive_genes <- false_positive_genes %>% select(entrez) %>% .$entrez
  sys_cans = sys_cans %>% mutate(NCG_FP=ifelse(entrez%in%false_positive_genes, TRUE, FALSE))

  return(sys_cans)

}


## Score genes
scoreGenes = function(ncg.tissue.name=NULL,
                      output.dir=NULL, exclude.features=NULL,
                      iters=NULL,
                      step=NULL,
                      top.rank=10,
                      refine=TRUE,
                      gtex.tissue="Esophagus"){
  require(e1071)
  require(caret)
  cat("Summarising cross-validation", "\n")
  cv_stats_fn = paste0(output.dir, "/cv_stats.tsv")
  cv_stats = read.table(cv_stats_fn, header = T, sep = "\t")

  outDir = paste0(output.dir, "/Results/")
  dir.create(outDir) ## Write results in a separate directory

  geneInfo=sysSVM:::geneInfo
  cancerGenes=sysSVM:::cancerGenes

  ## Check point for NCG tissues
  ncg_tissues = cancerGenes %>% select(primary_site) %>% unique %>% .$primary_site
  if (!(ncg.tissue.name%in%ncg_tissues)){
    stop("Please provide a valid tissue for NCG5!")
  }

  ## Load training and prediction sets
  load(paste(output.dir, "/training_set.Rdata", sep=""))
  load(paste(output.dir, "/prediction_set.Rdata", sep=""))
  load(paste(output.dir, "/prediction_set_noScale.Rdata", sep=""))

  ## What features should be excluded
  s = c("sample", "entrez", "no_ALL_muts","no_TRUNC_muts", "no_NTDam_muts",
        "no_GOF_muts", "Copy_number", "CNVLoss", "CNVGain")
  if(!is.null(exclude.features)){
    s = s[!(s%in%exclude.features)]
  }

  prediction_ns = prediction_ns %>% tibble::rownames_to_column() %>%
    separate(rowname, into=c("sample", "entrez"), sep="\\.") %>%
    mutate(entrez=as.numeric(entrez)) %>%
    select(one_of(s))

  ## Fix gene info table from NCG
  geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
  ## Get a cancer gene with all the associated primary sites and cancer sites
  cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
    group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
                                   cancer_site=paste(unique(cancer_site), collapse=",")) %>%
    ungroup
  geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))

  ## Training
  training = training%>%subset(type=="C")

  ## Filter out columns that are NaN - this should be changed to replace the column with all 0s
  training=Filter(function(x)!all(is.nan(x)), training)
  prediction=Filter(function(x)!all(is.nan(x)), prediction)

  ## Check the steps to summarise the results
  steps = seq(0,iters,step)
  steps = steps[2:length(steps)]
  if(steps[length(steps)]!=iters){steps[length(steps)]=iters}

  geneLists = list()
  recurrence_table = NULL

  for(st in steps){
    cat(paste0("Current step: Scoring at ", st, " iterations..."))
    outDir2 = paste0(output.dir, "/Results/", st, "/")
    dir.create(outDir2) ## Write results in a separate directory
    cv_stats_summary <- cv_stats %>% subset(!is.na(value) & type=="Sensitivity" & set=="test" & iteration<=st) %>%
      group_by(kernel, nu, gamma, degree) %>% summarize(iterations=n(),
                                                        min=min(value),
                                                        q1=quantile(value)[2],
                                                        median=median(value),
                                                        mean=mean(value),
                                                        q3=quantile(value)[4],
                                                        max=max(value),
                                                        var=stats::var(value)) %>% ungroup
    cv_stats_summary$analysis = paste(cv_stats_summary$kernel, cv_stats_summary$nu, cv_stats_summary$gamma, cv_stats_summary$degree, sep=".")

    write.table(cv_stats_summary, file=paste0(outDir2,"/cv_stats_summary.tsv"), quote=F, row.names=F, sep="\t")

    ## -------------------
    ## Get the best models
    ## -------------------
    best_models = NULL
    for (k in unique(cv_stats_summary$kernel)){
      d = cv_stats_summary %>% subset(kernel==k) %>% data.frame()
      bm = d %>% arrange(desc(mean)) %>% slice(1:5) %>% mutate(var_rank=row_number(var)) %>% subset(var_rank==1)
      best_models = rbind(best_models, bm)
    }
    write.table(best_models, file=paste0(outDir2,"/best_models.tsv"), quote=F, row.names=F, sep="\t")

    for (i in 1:nrow(best_models)){
      kern = as.character(best_models$kernel[i])
      analysis = best_models$analysis[i]
      mynu=best_models$nu[i]
      mygamma=best_models$gamma[i]
      mydegree=best_models$degree[i]

      if (kern=="linear"){
        svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, nu=mynu)
      }else if (kern=="radial"){
        svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, gamma=mygamma, nu=mynu)
      }else if (kern=="polynomial"){
        svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, gamma=mygamma, nu=mynu, degree=mydegree)
      }else if (kern=="sigmoid"){
        svm.model <- svm(type ~ ., data = training, kernel=kern, type="one-classification", scale=FALSE, gamma=mygamma, nu=mynu)
      }

      ## Asses training sensitivity
      trainset = training %>% select(type)
      trainset$fitted = svm.model$fitted
      trainset = trainset %>% select(type, fitted)
      trainset = trainset %>% mutate(f=ifelse(fitted==TRUE, "C", "N"))
      trainset$type = factor(as.character(trainset$type), levels = c("C", "N"))
      trainset$f = factor(as.character(trainset$f), levels = c("C", "N"))
      ## Confusion matrix
      cM <- confusionMatrix(trainset$f,trainset$type, positive = c("C"))
      training_sensitivity = cM$byClass["Sensitivity"]

      ## make prediction
      predictions = predict(svm.model, prediction, decision.values = TRUE)
      ## add probability from Platt and Binning
      predictions_platt = plattScaling(svm.model, predictions)
      colnames(predictions_platt) = c("dv", "label", "prob_platt")
      preds = predictions_platt %>% tibble::rownames_to_column()

      ## Take information from gene info
      preds = preds %>%
        separate(rowname, into=c("sample", "entrez"), sep="\\.") %>%
        mutate(entrez=as.numeric(entrez)) %>%
        left_join(geneInfo%>%rename(gene_type=cancer_type), by=c("entrez"))

      ## Join also with systems level properties from validaton set
      preds = preds %>% left_join(prediction_ns, by=c("sample", "entrez"))


      save(svm.model, file=paste(outDir2,"/model_",analysis, ".Rdata", sep=""))
      save(preds, file=paste(outDir2,"/predictions_",analysis, ".Rdata", sep=""))
    }

    ## Get the scores (based on the best models above)
    LINEAR_PREDS = best_models$analysis[best_models$kernel=="linear"]
    LINEAR_CVS = round(best_models$mean[best_models$kernel=="linear"], digits = 2)
    LINEAR_VAR = round(best_models$var[best_models$kernel=="linear"], digits = 2)

    POLYNOMIAL_PREDS = best_models$analysis[best_models$kernel=="polynomial"]
    POLYNOMIAL_CVS = round(best_models$mean[best_models$kernel=="polynomial"], digits = 2)
    POLYNOMIAL_VAR = round(best_models$var[best_models$kernel=="polynomial"], digits = 2)

    RADIAL_PREDS = best_models$analysis[best_models$kernel=="radial"]
    RADIAL_CVS = round(best_models$mean[best_models$kernel=="radial"], digits = 2)
    RADIAL_VAR = round(best_models$var[best_models$kernel=="radial"], digits = 2)

    SIGMOID_PREDS = best_models$analysis[best_models$kernel=="sigmoid"]
    SIGMOID_CVS = round(best_models$mean[best_models$kernel=="sigmoid"], digits = 2)
    SIGMOID_VAR = round(best_models$var[best_models$kernel=="sigmoid"], digits = 2)

    #save.image(file = paste0(output.dir, "/namespace.Rdata"))

    scores = getScore(path=outDir2,
                      linear_preds = paste0("/predictions_", LINEAR_PREDS, ".Rdata"), linearCVS = LINEAR_CVS, linearVar = LINEAR_VAR,
                      radial_preds = paste0("/predictions_", RADIAL_PREDS, ".Rdata"), radialCVS = RADIAL_CVS, radialVar = RADIAL_VAR,
                      polynomial_preds = paste0("/predictions_", POLYNOMIAL_PREDS, ".Rdata"), polynomialCVS = POLYNOMIAL_CVS, polynomialVar = POLYNOMIAL_VAR,
                      sigmoid_preds = paste0("/predictions_", SIGMOID_PREDS, ".Rdata"), sigmoidCVS = SIGMOID_CVS, sigmoidVar = SIGMOID_VAR)

    save(scores, file=paste0(outDir2, "/scores.Rdata"))

    if(refine){
      ## Add refined CGCs, exclude "Not expressed" genes and define patient ranks
      ## Then I use the getSysCans (source it from config.R) to get a data frame with the ranks. Command as follows:
      syscan = getSysCansOAC(output.dir = output.dir, working.dir = outDir2, gtex.tissue = gtex.tissue)
      ## And then save it
      save(syscan, file=paste0(outDir2, "/syscan.Rdata"))
      cat("\n")
    }## proper function to calculate ranks will be added here (instead of getSysCansOAC)

    ## Check genes and decide whether to stop or not
    top_syscans = syscan %>% subset(patient.rank<=top.rank)
    geneLists[[as.character(st)]] = unique(top_syscans$entrez)
  }

  unique_predictions = NULL
  count = 1
  for(ul in unique(geneLists)){
    l = lapply(geneLists, function(x) sum(x==ul)==length(ul))
    oc = sum(unlist(l))
    iterations_predicted = paste0(names(which(l ==TRUE)), collapse=",")
    d = data.frame(list = count, occurences = oc, iterations=iterations_predicted)
    unique_predictions = rbind(unique_predictions, d)
    count = count + 1
  }

  unique_predictions = unique_predictions %>% dplyr::arrange(desc(occurences)) %>% dplyr::mutate(steps=length(steps), frequency=occurences/steps)
  write.table(unique_predictions, file=paste0(output.dir, "/gene_prediction_frequency.tsv"), row.names = F, sep="\t", quote = F)
  lapply(unique_predictions, class)
  ## And finally export the most frequent list (from the bigest number of iterations to approximate the true values of scores/sensitivity of the models)
  check_max = which(unique_predictions$frequency==max(unique_predictions$frequency))
  find_iter = as.numeric(unlist(strsplit(as.character(unique_predictions$iterations[check_max]), ",")))
  find_iter = max(find_iter)
  file.copy(paste0(output.dir, "/Results/", find_iter, "/syscan.Rdata"), paste0(output.dir, "/syscan.Rdata"))
  file.copy(paste0(output.dir, "/Results/", find_iter, "/scores.Rdata"), paste0(output.dir, "/scores.Rdata"))


}


## Prepare training and prediction set
#############
## Functions
#############

## Marking of training/prediction set
markTrainCGC <- function(df, fp_dir="example_data/NCG_false_positives.txt",
                         geneInfo_fn="example_data/geneInfoNCG5_2.Rdata",
                         cancerGenes_fn="example_data/cancerGenesNCG5_2.Rdata"){

  load(geneInfo_fn)
  load(cancerGenes_fn)

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
  false_positive_genes <- read.table(fp_dir, header = T, sep = "\t")
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
createDescribeTrainingCGC = function( input.file=NULL, output.dir=NULL, exclude.features=NULL, trainingMode=TRUE, models=NULL, scaling.factors=NULL){

  ## Get data
  if(!is.null(input.file)){
    d = read.table(input.file, sep="\t", header = T)

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
    ## Save training and validation set without scaling
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

  ## Constants
  NU = seq(0.05, 0.9, 0.05)
  GAMMA = 2^seq(-7,4)
  DEGREE = c(3, 4, 9)

  LO=1 ## Leave out i.e here I take CV-1 folds for the train and 1 for the test
  CV=cv

  ## Start with linear
  cat("Grid search for:", "\n")

  create.output.folder = function (output.dir, kern, mynu, mygamma, mydegree, no) {

    if(length(grep("Results",list.dirs(output.dir, recursive = F, full.names=F)))==0){
      dir.create(paste(output.dir, "/Results", sep=""))
    }
    timer    = format(Sys.time(), "%Y-%m-%d.%H_%M_%S")
    if(kern=="linear"){
      job_dir = paste(output.dir, "/Results/", kern,".",mynu, sep="" )
      if(!file.exists(job_dir)){
        dir.create(job_dir)
      }
    }else if (kern=="radial"){
      job_dir = paste(output.dir, "/Results/", kern,".", mynu, ".", mygamma, sep="")
      if(!file.exists(job_dir)){
        dir.create(job_dir)
      }
    }else if (kern=="polynomial"){
      job_dir = paste(output.dir, "/Results/", kern,".", mynu, ".", mygamma, ".", mydegree, sep="")
      if(!file.exists(job_dir)){
        dir.create(job_dir)
      }
    }else if (kern=="sigmoid"){
      job_dir = paste(output.dir, "/Results/", kern,".", mynu, ".", mygamma, sep="")
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

                              svm.model <- svm(type ~ ., data = trainset, kernel=k, type="one-classification", scale=FALSE, nu=mynu)
                              save(svm.model, file=paste(iteration_dir, "/svmModel.Rdata", sep=""))

                              prediction <- predict(svm.model, testset%>%select(-type), decision.values = TRUE, probability = TRUE)
                              ## save prediction in a file
                              save(prediction, file=paste(iteration_dir, "/prediction.Rdata", sep=""))

                              ## Get the performance
                              noTrain = trainset %>% nrow
                              cv_analysis = tail(unlist(strsplit(iteration_dir, "\\/")), n=1)
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
                              cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, kernel=k,
                                                                     nu=mynu, gamma=mygamma, degree=mydegree,
                                                                     iteration=y,
                                                                     type="accuracy", trainSize=noTrain,
                                                                     class="overall", set = "test",
                                                                     value=cM$overall[["Accuracy"]]))
                              cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, kernel=k,
                                                                     nu=mynu, gamma=mygamma, degree=mydegree,
                                                                     iteration=y,
                                                                     type="Sensitivity", trainSize=noTrain,
                                                                     class="overall", set = "test",
                                                                     value=cM$byClass["Sensitivity"]))
                              cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, kernel=k,
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
                              cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, kernel=k,
                                                                     nu=mynu, gamma=mygamma, degree=mydegree,
                                                                     iteration=y,
                                                                     type="accuracy", trainSize=noTrain,
                                                                     class="overall", set = "train",
                                                                     value=cM$overall[["Accuracy"]]))
                              cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, kernel=k,
                                                                     nu=mynu, gamma=mygamma, degree=mydegree,
                                                                     iteration=y,
                                                                     type="Sensitivity", trainSize=noTrain,
                                                                     class="overall", set = "train",
                                                                     value=cM$byClass["Sensitivity"]))
                              cv_stats <- rbind(cv_stats, data.frame(analysis=cv_analysis, kernel=k,
                                                                     nu=mynu, gamma=mygamma, degree=mydegree,
                                                                     iteration=y,
                                                                     type="Specificity", trainSize=noTrain,
                                                                     class="overall", set = "train",
                                                                     value=cM$byClass["Specificity"]))
                              return(cv_stats)
                            }

                        }
      close(pb)
      stopCluster(cl)
      cv_stats_full = rbind(cv_stats_full, result)
    }

  }
  ## write metrics
  write.table(cv_stats_full, file=paste(output.dir,"/cv_stats.tsv",sep=""), quote=F, row.names = F, sep="\t")
}





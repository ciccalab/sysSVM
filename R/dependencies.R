pkgsR  =  c('plyr',
            'tidyr',
            'caret',
            'e1071',
            'sampling',
            'doSNOW',
            'dplyr',
            'readxl')
for (pkgR in pkgsR){
  if (!pkgR %in% rownames(installed.packages())) {
    install.packages(pkgR, dependencies=TRUE, INSTALL_opts = c('--no-lock'), repos = "https://cran.ma.imperial.ac.uk/")
    library(pkgR, character.only = TRUE, quietly = TRUE)
  } else {
    library(pkgR, character.only = TRUE, quietly = TRUE)
  }
}

# source("http://bioconductor.org/biocLite.R")
#
# pkgs  =  c("RDRToolbox")
# for (pkg in pkgs){
#   if (!pkg %in% rownames(installed.packages())) {
#     biocLite(pkg)
#     library(pkg, character.only = TRUE, quietly = TRUE)
#   } else {
#     library(pkg, character.only = TRUE, quietly = TRUE)
#   }
# }

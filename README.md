# sysSVM
Patient-specific cancer driver prediction using support vector machines and systems biology

### Example data were imported using the following command
```r
> install.packages("devtools")
> library(devtools)
> devtools::use_data(cgcs, cancerGenes, previous, gtex, geneInfo, geneProperties_mmImputed, false_positive_genes, oac_data, internal = TRUE)
```

### Install 

```r
> install.packages("devtools")
> library(devtools)
> install_github("thmourikis/sysSVM")
```

#### Run test data ####
```r
> library(sysSVM)
> options(stringsAsFactors = F)
> syssvm()
```

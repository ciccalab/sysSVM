# sysSVM
Patient-specific cancer driver prediction using support vector machines and systems biology

### Install ### (when it goes public)

```r
> install.packages("devtools")
> library(devtools)
> install_github("thmourikis/sysSVM")
```

#### Run test data ####
```r
> syssvm(input.file = "example_data/oac_ML_input_formatted.tsv", output.dir = "OAC", exclude.features = c("young", "no_ALL_muts"), cores=2, iters = 10, kernels=c("linear"), ncg.tissue.name = "esophagus", reproduce=T)
```

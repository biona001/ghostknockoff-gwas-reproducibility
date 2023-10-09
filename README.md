# GhostKnockoff GWAS reproducibility

This repo contains source codes used in the paper **In silico identification of putative causal genetic variants** by He et al. 

Note that code here is *purely for review purposes*. We are preparing a more robust and computationally efficient pipeline in an alternate repo which will be available within the next 1-2 months. The final repo will contain a software pipeline that is more general, easy to install, user-friendly, and scalable than the source code here. 

## Example: Running the Alzheimers Diseases analyses

```shell
$ Rscript --vanilla GKL_RunAnalysis_All.R arg1 arg2 arg3 arg4 arg5
```
where 
+ `arg1`: study number (integer between 1 to 20 representing different Alzhimers disease studies)
+ `arg2`: path to Z score file which includes Z scores as well as the chr/pos and ref/alt alleles
+ `arg3`: path to summary statistics file for the different AD studies (in particular, including sample size information)
+ `arg4`: path to the cS2G file, which maps each SNP to the closet gene
+ `arg5`: Output directory

For example, on the Sherlock cluster, one can execute:

```shell
$ Rscript --vanilla GKL_RunAnalysis_All.R 1 /oak/stanford/groups/zihuai/ESGWAS_lasso/Z_scores/AD_meta/AD_Zscores_Meta.txt /oak/stanford/groups/zihuai/ESGWAS_lasso/AD_Analysis/SummaryStatInfo.txt /oak/stanford/groups/zihuai/XinranQi/cS2G_UKBB/topcS2G_allVariants/topcS2GGene_allVariants.csv /scratch/users/bbchu/AD_meta/Results/
```

## Dependnecies (concise)

The pipeline was tested on Stanford's Sherlock cluster which loads the following modules 
+ `R/4.0.2`
+ `cmake/3.24.2`
+ `harfbuzz/1.4.8`
+ `fribidi/1.0.12`
+ `libgit2/1.1.0`
+ `openssl/3.0.7`

The main GhostBasil pipeline is implemented in the `GKL_RunAnalysis_All.R` file, which depends on the following `R` packages

+ `data.table` v1.14.8
+ `ghostbasil` v0.1.3 (This is currently not publicaly available)
+ `Matrix` v1.6-0
+ `susieR` 0.12.35
+ `rhdf5` v2.34.0

## Full `sessionInfo()` output:

```R
> sessionInfo()
R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /share/software/user/open/openblas/0.2.19/lib/libopenblasp-r0.2.19.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] susieR_0.12.35    rhdf5_2.34.0      ghostbasil_0.1.3  Matrix_1.6-0     
[5] data.table_1.14.8

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.10        magrittr_2.0.3     tidyselect_1.2.0   munsell_0.5.0     
 [5] colorspace_2.1-0   lattice_0.20-41    R6_2.5.1           rlang_1.1.0       
 [9] fansi_1.0.4        plyr_1.8.9         dplyr_1.1.1        tools_4.0.2       
[13] grid_4.0.2         gtable_0.3.3       utf8_1.2.3         irlba_2.3.5.1     
[17] cli_3.6.1          matrixStats_0.63.0 tibble_3.2.1       lifecycle_1.0.3   
[21] crayon_1.5.2       mixsqp_0.3-48      Rhdf5lib_1.12.1    ggplot2_3.4.2     
[25] vctrs_0.6.1        rhdf5filters_1.2.1 codetools_0.2-16   glue_1.6.2        
[29] pillar_1.9.0       compiler_4.0.2     generics_0.1.3     scales_1.2.1      
[33] reshape_0.8.9      pkgconfig_2.0.3
```

## Data dependencies

+ The pipeline for the meta-analysis of Alzheimers Diseases requires 
    - [AD_Zscores_Meta.txt]() (1.8GB)
    - [topcS2GGene_allVariants.csv]() (143MB), and 
    - [pre-computed knockoff statistics]() (15GB) (please see **Knockoff generation** section below for details)
as inputs. Please download the files and put them in the folder `ghostknockoff-gwas-reproducibility/data`.
+ Note the primary pipeline implemented in `GKL_RunAnalysis_All.R` additionally requires as input a list of typed SNPs present on the UK Biobank genotyping chip. At the moment, this input file is *hard coded* into the source code, because we are not sure if we can distribute the list of typed variants freely. 

## Knockoff generation

+ Processing of LD panels (including downloading and importing the data matrices) is carried out by an independent software [EasyLD.jl](https://github.com/biona001/EasyLD.jl).
+ Knockoff optimization problem was carried out by [Knockoffs.jl](https://github.com/biona001/Knockoffs.jl). 
+ Note that in our final software release, pre-computed knockoff statistics will be made available for download, so users will not have to manually install `EasyLD.jl` nor `Knockoffs.jl` to carry out this step.

## Contact

For questions regarding this repo, please reach out to Zihuai He (`zihuai@stanford.edu`) or Benjamin Chu (`bbchu@stanford.edu`). 

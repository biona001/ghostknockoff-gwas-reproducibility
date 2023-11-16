# GhostKnockoff GWAS reproducibility

This repo contains source codes used in the paper **In silico identification of putative causal genetic variants** by He et al. 

Note that code here is *purely for review purposes*. We are preparing a more robust and computationally efficient pipeline in an alternate repo which will be available within the next 1-2 months. The final repo will contain a software pipeline that is more general, easy to install, user-friendly, and scalable than the source code here. 

## Installation

To run the provided examples, one need to 
1. Install all required `R` packages
2. Download pre-computed knockoff statistics
They steps are described below. 

### Installing `R` dependencies

To run the provided scripts, we need to first install the [ghostbasil](https://github.com/JamesYang007/ghostbasil) `R` package. This involves:

1. Clone the repository 
```
git clone https://github.com/JamesYang007/ghostbasil
```
2. (optional) on cluster environments, we found it necessary to load the following modules
```
module load mpfr/4.1.0 gcc/12.1.0 libgit2/1.1.0 system harfbuzz fribidi
```
3. (optional) for Mac users, 
```
brew install libomp
```
4. Within `R`, install the `ghostbasil` package
```
$ cd ghostbasil/R
$ R
> library(devtools)
> install()
```

Please also install the following `R` packages

+ `data.table`
+ `Matrix`
+ `susieR`
+ `rhdf5`
+ `plyr`
+ `dplyr`
+ `CMplot`

### Download Required Data

Please [download this data](https://drive.google.com/file/d/1_ajlxFWE2MCSgBXDgDbeZh9Lq721WANA/view?usp=drive_link) and unzip it:
```
unzip data.zip
```
After unzipping, you will find the following files:
- `EUR` directory (9.3G) contains pre-computed knockoff statistics (for EUR ancestry) stored in `.h5` format as well as summaries for each block (please see **Knockoff generation** section below for details)
- `AD_Zscores_Meta.txt` (1.8GB) contains study-specific Z scores for each SNP as well as basic allelic information (CHR/POS/REF/ALT/...etc)
- `topcS2GGene_allVariants.csv` (143MB) contains the nearest gene for each SNP
- `SummaryStatInfo.txt` (4KB) contains summaries for the 10 Alzheimer Disease studies (sample size, human genome build...etc)

**Notes:**
+ Data is stored on google drive and does NOT support `wget`. Please download it manually on the browser. Alternatively, you can try installing [gdown](https://stackoverflow.com/questions/25010369/wget-curl-large-file-from-google-drive).
+ Because the files are large, sometimes unzipping throws a warning `error: invalid zip file with overlapped components (possible zip bomb)`. Please try with 
```
UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE unzip data.zip
```

## Example: Running the Alzheimers Diseases analyses

```shell
$ Rscript --vanilla GKL_RunAnalysis_All.R arg1 arg2 arg3 arg4 arg5 arg6
```
where 
+ `arg1`: Integer between 1 to 11 representing different Alzhimers disease studies. The first 10 studies are summarized in `data/SummaryStatInfo.txt` while study 11 is a meta-analysis of the 10 studies.
+ `arg2`: path to pre-computed knockoff statistics (see item 1 under [Download Required Data](https://github.com/biona001/ghostknockoff-gwas-reproducibility#download-required-data))
+ `arg3`: path to Z score file which includes Z scores as well as the chr/pos and ref/alt alleles (see item 2 under [Download Required Data](https://github.com/biona001/ghostknockoff-gwas-reproducibility#download-required-data))
+ `arg4`: path to summary statistics file for the different AD studies (see item 3 under [Download Required Data](https://github.com/biona001/ghostknockoff-gwas-reproducibility#download-required-data))
+ `arg5`: path to the cS2G file, which maps each SNP to the closet gene (see item 4 under [Download Required Data](https://github.com/biona001/ghostknockoff-gwas-reproducibility#data-dependencies))
+ `arg6`: Output directory

For example, to run the meta-analysis GWAS result for Alzheimers Disease:

```shell
$ Rscript --vanilla GKL_RunAnalysis_All.R 11 data/EUR data/AD_Zscores_Meta.txt data/SummaryStatInfo.txt data/topcS2GGene_allVariants.csv Results
```

## Dependencies

The pipeline was tested on Stanford's Sherlock cluster which loads the following modules 
+ `R/4.0.2`
+ `cmake/3.24.2`
+ `harfbuzz/1.4.8`
+ `fribidi/1.0.12`
+ `libgit2/1.1.0`
+ `openssl/3.0.7`

The main GhostBasil pipeline is implemented in the `GKL_RunAnalysis_All.R` file, which depends on the following `R` packages

+ `data.table` v1.14.8
+ `ghostbasil` v0.1.3
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

## Knockoff generation

+ Processing of LD panels (including downloading and importing the data matrices) is carried out by an independent software [EasyLD.jl](https://github.com/biona001/EasyLD.jl).
+ Knockoff optimization problem was carried out by [Knockoffs.jl](https://github.com/biona001/Knockoffs.jl). 
+ Note that in our final software release, pre-computed knockoff statistics will be made available for download, so users will not have to manually install `EasyLD.jl` nor `Knockoffs.jl` to carry out this step.

## Contact

For questions regarding this repo, please reach out to Zihuai He (`zihuai@stanford.edu`) or Benjamin Chu (`bbchu@stanford.edu`). 

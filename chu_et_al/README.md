# Group Knockoff Reproducibility

This repo contains source codes to reproduce the simulated and real data experiments used in the paper [Second-order group knockoffs with applications to GWAS](https://arxiv.org/abs/2310.15069). All scripts are provided as individual jupyter notebooks, which can be viewed directly in the web browser. To execute the scripts, one can do so directly within jupyter notebooks, or copy the code into the Julia REPL. 

Note that these source codes are *purely for review and reproducibility purposes*. Users interested in running the overall pipeline should proceed to the main [GhostKnockoffGWAS page](https://github.com/biona001/GhostKnockoffGWAS). 

## Software versions

We developed a number of softwares to enable group knockoff simulation and analysis. Installing these Julia packages will automatically install the relevant Julia dependencies.

+ `Knockoffs.jl` (v1.1.1): https://github.com/biona001/Knockoffs.jl
+ `EasyLD.jl` (pre-release): https://github.com/biona001/EasyLD.jl
+ `GhostKnockoffGWAS.jl` (pre-release): https://github.com/biona001/GhostKnockoffGWAS

For the real-data analysis, we also used 

+ `ghostbasil` (v0.1.3): https://github.com/JamesYang007/ghostbasil
+ `liftOver` (v1.14.0): https://bioconductor.org/packages/release/workflows/html/liftOver.html

We used `Julia` v1.8.4 and `R` v4.0.2 for all analyses. 

## Simulated experiments

All scripts (including scripts to make figures) are available in jupyter notebooks starting with `compare`.

+ Figure 1 is a combination of 5 simulations, where the simulations are provided in standalone notebooks (e.g. `compare_group_knockoff_AR1.ipynb`), while the figure is generated in `compare_group_knockoff_figure1.ipynb`. 
+ Figure 2 code is provided in `compare_group_knockoff_panUKB3.ipynb`
+ Figure 3 code is provided in `compare_exchangeability.ipynb`
+ Figure 4 and Table 1 code is provided in `compare_group_knockoff_timings.ipynb`
+ Figure 5 code is detailed below
+ Figure S1 code is provided in `compare_group_knockoff_group_vs_ungroup.ipynb`
+ Figure S2 code is provided in `compare_group_knockoff_panUKB_block_sizes.ipynb`
+ Figure S3 code is provided in `compare_group_knockoff_panUKB.ipynb`
+ Figure S4, similar to Figure 1, is a combination of 5 simulations. Individual simulation code is provided in standalone notesbooks ending in `_marginal.ipynb`
+ Table S1: each entry was searched manually on the [GWAS catalog](https://www.ebi.ac.uk/gwas/)
+ Figure S5 code is provided in `compare_group_knockoff_trace_and_mineval.ipynb`

## Albuminuria GWAS

The analysis shown in Figure 5 is conducted in several steps

+ `ghostknockoff-part0.ipynb`: Downloads the Pan-UKB LD matrices 
+ `ghostknockoff-part1.ipynb`: Reads the downloaded LD matrices and solves the knockoff optimization problem across 1703 genomic regions 
+ `ghostknockoff-part2.ipynb`: Reads Albuminuria Z-scores and intermediate knockoff results (e.g. S matrix) from part 1. Then it generates knockoff Z-scores and performs knockoff filter
+ `ghostknockoff-manhattan.ipynb`: Creates the Manhattan plot as featured in our paper

## Questions?

If anything is unclear, or you cannot reproduce our results, please file an issue on GitHub and I will take a look at it asap. You can also reach me at `bbchu@stanford.edu`.

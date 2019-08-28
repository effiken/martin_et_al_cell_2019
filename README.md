## Making the figures
### Requirements

Tested on macOS

1. R
2. R packages: 
   - viridis
   - RColorBrewer
   - gplots
   - Matrix
   - Matrix.utils
   - seriation
   - [scDissector](https://github.com/effiken/scDissector)
3. Downloaded and unzipped version of this repository  on a local path.


### Running the scripts in R

1. Assuming martin_et_al_cell_2019 is the local path of the repository we need to load the script files:

`source("martin_et_al_cell_2019/scripts/figures_main.R")`

2. And then make figures 1-5 and supplumentary figures:

`make_martin_et_al_figures("martin_et_al_cell_2019/",download_data = F)`


### Output

Figure will be generated in:
  - martin_et_al_cell_2019/output/main_figures/
  - martin_et_al_cell_2019/output/supp_figures/
  
Tables will be generated in:
  - martin_et_al_cell_2019/output/tables/

## Clustering

### Requirements

Tested on linux LSF HPC. Due to lack of support of some of the depdendencies, the script cannot run on macOS.

1. R
2. R packages:
   - Matrix
   - Matrix.utils
   - gplots
   - seriation
   - [tglkmeans](https://bitbucket.org/tanaylab/tglkmeans)
   - [scDissector](https://github.com/effiken/scDissector)
3. Downloaded and unzipped version of this repository  on a local path.

### Running the scripts in R

Assuming martin_et_al_cell_2019 is the local path of the repository, the following script will run the clustering distributedly on LSF:

`source("martin_et_al_cell_2019/scripts/clustering/run_clustering_ileum.r")`

Alternatively, clustering can be run locally:

`source("martin_et_al_cell_2019/scripts/clustering/run_clustering_ileum_local.r")`

Note: Each run of the clustering might produce slightly different results due to different random seeds.

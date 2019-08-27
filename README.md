# Making the figures
## Requirements

1. R
2. R packages: 
   - viridis
   - RColorBrewer
   - RColorBrewer
   - gplots
   - Matrix
   - Matrix.utils
   - seriation
   - [scDissector](https://github.com/effiken/scDissector)
3. Downloaded and unzipped version of this repository  on a local path.


## Running the scripts in R

1. Assuming GitHub/martin_et_al_cell_2019 is the local path of the repository we need to load the script files:

`source("GitHub/martin_et_al_cell_2019/scripts/figures_main.R")`

2. And then make figures 1-5 and supplumentary figures:

`make_martin_et_al_figures("~/Downloads/martin_et_al_cell_2019/",download_data = F)`




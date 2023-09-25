## R-codes
for preprint egusphere-2023-1267, https://doi.org/10.31223/X54S90

The codes are written by Takahito Mitsui. Some parts are written following earlywarnings package in R libraries by Vasilis Dakos: https://github.com/earlywarningtoolbox/earlywarnings-R 

### ews_data.R  for Figs 1, 2 & 3
1. This code uses Rasmussen_et_al_2014_QSR_Table_2.xlsx and seierstad.xlsx for input data
2. make a directry/folder named 'ews_data' at the same directry level of ews_data.R
3. install R-libraries listed in the top of the R-script if needed   
4. change the number 'k' (k=1, 2, ..., 12) to specify which proxy and interstedials are analyzed.
5. choose parameters for calculating CSD indicators (variance, lag.1 autocorrelation, lambda, etc): ns (the number of surrogate data), smoothing (gaussian or loess), bandwidth (0-100%), winsize (0-100%)     
6. run ews_data.R. Output Figures are created in 'ews_data'

### ews_data_summary.R  for Fig 4
- The basic elemnts in this code are the same as ews_data.R
- make a directry/folder named 'ews_data_summary' at the same directry level of ews_data_summary.R
- change the number 'k' (k=1, 2, ..., 6) to specify which proxy is analyzed.
- run ews_data_summary.R. Output Figures are created in 'ews_data_summary'
  
### R-Codes  for Fig. 5 
are included in several folders:
- schematic_stommel2 (Stommel model Fig. 5a)
- schematic_fhn2  (Stochastic slow-fast oscillations in Figs 5b & 5b)
- schematic_hopf2 (Hopf bifurcation mechanism in Figs 5d & 5e)
- schematic_mmo2  (Mixed-mode oscillation mechanism in Figs 5f & 5g)

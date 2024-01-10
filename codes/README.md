## R-codes
for preprint egusphere-2023-1267, https://doi.org/10.31223/X54S90

The codes are written by Takahito Mitsui. Some parts are written following earlywarnings package in R libraries by Vasilis Dakos: https://github.com/earlywarningtoolbox/earlywarnings-R 

### ews_data.R  for Figs 1, 2 and 3 (and Figs S1-S10)
1. This code uses Rasmussen_et_al_2014_QSR_Table_2.xlsx and seierstad.xlsx for input data
2. make a directry/folder named 'ews_data' at the same directry level of ews_data.R
3. install R-libraries listed in the top of the R-script if needed   
4. change the number 'k' (k=1, 2, ..., 12) to specify which proxy and interstedials are analyzed. For example, k=1 corresponds to d18O data from GI-25 to GI-20 (Fig. 2), and k=2 corresponds to d18O data from GI-19.2 to GI-1
5. choose parameters for calculating CSD indicators (variance, lag.1 autocorrelation, lambda, etc): ns (the number of surrogate data), smoothing (gaussian or loess), bandwidth (0-100%), winsize (0-100%)     
6. run ews_data.R. Output Figures are created in 'ews_data'

### ews_data_summary.R  for Fig 4 (and Figs S11-S22)
- The basic elemnts in this code are the same as ews_data.R
- make a directry/folder named 'ews_data_summary' at the same directry level of ews_data_summary.R
- change the number 'k' (k=1, 2, ..., 6) to specify which proxy is analyzed.
- run ews_data_summary.R. Output Figures 'result*.dat' are created in 'ews_data_summary'
- Then generate Fig 4d (result_matrix.eps) by result.R in ews_data_summary
  
### R-Codes for Fig. 5 
are included in several folders as follows:
- schematic_stommel2/schematic_stommel.R (Stommel model Fig. 5a)
- schematic_fhn2/schematic_fhn.R  (Stochastic slow-fast oscillations in Figs 5b & 5b)
- schematic_hopf2/schematic_hopf.R (Hopf bifurcation mechanism in Figs 5d & 5e)
- schematic_mmo2/schematic_mmo.R  (Mixed-mode oscillation mechanism in Figs 5f & 5g)

### R-Code for Fig. 6 
- ews_stommel_rate.R in folder schematic_stommel2 

### ews_data_summary_further_test.R for Appendix A 
- Program to calculate the probability of observing robust precursor signals for 5000 phase-randomized surrogtes in Appendix A
- The outputs of this R-code (like proxy1_DO**_case*.dat) are generated in the directry ews_data_summary_further_test 
- Step 1. In the code, set index 'i' to specify which DO is examined: i=2 for DO-25 or i=20 for DO-12 (those two are considered in the paper)  
- Step 2. In the code, set index 'case' from 1 to 10 (this just splits 5000 experiments into ten 500 experiments allowing 'parallel' computations). Then run this R-code for each case.
- Step 3. Finally go to directry ews_data_summary_further_test. Then run result.R to obtain the probability of observing robust precursor signals for variance: the shown numbers are probabilites of getting robust SPS for the variance (>=15 times over 30 parameter sets), lag-1AC (>=15), variance (>15) and lag-1 AC (>15), respectively.


### Table S1
- is generated from ews_data_summary/duration.R

### ews_data_summary_rebound.R for Fig. S23
- make a directry/folder named 'ews_data_summary_rebound' at the same directry level of ews_data_summary_rebound.R
- change the number 'k' (k=1, 2, ..., 6) to specify which proxy is analyzed.
- run ews_data_summary_rebound.R. Output Figures 'result*.dat' are created in 'ews_data_summary_rebound'
- Then generate result_matrix_rebound.eps by result.R in ews_data_summary_rebound

### ews_data_5yr.R for Figs S24-S28
- ews_data_5yr.R generates Figs S25-S28 in the directory ews_data_5yr
- change k to select a proxy and interstadials from NGRIP d18O and log Ca2+

### R-codes for Figs S29-S31
are included in several folders as follows:
- schematic_fhn2/ews_fhn.R  (Fig. S29)
- schematic_hopf2/ews_hopf.R (Fig. S30)
- schematic_mmo2/ews_mmo.R  (Fig. S31)

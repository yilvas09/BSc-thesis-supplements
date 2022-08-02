# Supplements for my BSc thesis

This repository contains supplementing `MATLAB` codes and figures of my BSc thesis (Effect of Stochastic Demand on Exit and Entry - Applying Real Option Theory), where 2 different stochasic dynamics of market demand (**geometric Brownian motion** and **mean-reverting process**) are studied.

## Geometric Brownian Market Demand

The `gBm_baseline` contains a equilibrium market model with geometric Brownian demand. Codes and figures in this folder serve as a benchmark for the following mean-reverting case and are replications of models proposed by *Pindyck (2009)* and *Bustamante and Donangelo (2017)*.

To compute entry and exit threshold in a geometric Brownian market:
* `plot_threhold.m`

To compute the value of incumbent firms:
* `plot_val_func_GBM.m`

To plot risk-loading functions w.r.t. omega under different degrees of concentration (i.e. different numbers of incumbents) in *Bustamante and Donangelo (2017)*:
* `plot_betas_GBM.m`

`mean_reverting/figures` contains output figures of the scripts.

## Mean-reverting Market Demand

The `mean_reverting` directory contains my main contribution of studying mean-reverting market demands.

To perform comparative static analysis in mean-reverting markets for each parameter:
*	`Batch_E_MRV.m` : comparative statics for exit sunk cost, $E$
*	`Batch_eta_MRV.m` : comparative statics for reverting speed, $\eta$
*	`Batch_I_MRV.m` : comparative statics for entry sunk cost, $I$
*	`Batch_n_MRV.m` : comparative statics for number of active firms, $n$
*	`Batch_ph_MRV.m` : comparative statics for fixed operating cost, $\phi_i$
*	`Batch_sgm_MRV.m` : comparative statics for volatility rate of market size, $\sigma$
*	`Batch_tht_st_MRV.m` : comparative statics for mean market size, $\theta^*$

To perform additional numerical experiments discussed in appendix:
*	`Experiment_1.m` : replicate oligopoly effect
*	`Experiment_2.m` : replicate dilution effect

To replicate figure that plots entry and exit thresholds under different numbers of active firms:
*	`Plot_Threshold_MRV.m`

To replicate figure that plots value functions under different numbers of active firms:
*	`Plot_Value_functions_MRV.m`

To replicate figure that plots systematic risk loadings of incumbents under different numbers of active firms:
*	`Plot_Betas_MRV.m`

`mean_reverting/auxilaries` contains subroutines for the above `MATLAB` scripts.

`mean_reverting/figures` contains output figures of the scripts.

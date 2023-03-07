# Methods to link population data to climatic drivers 

Code linked to the preprint [Predictive ability of climate on vital rates and population dynamics depends on vital rate and model complexity](https://doi.org/10.1101/2022.03.11.484031)

Key points to use this repository:

* This repository contains all of the work done on this project since fall 2017. To establish a stronger link to the preprint, start from file _wrapper.R_. 
* We did not include the data to reproduce the results in this repository, because it occupies _XXXXX MB_. To obtain these model results, please contact the corresponding author, or wait for the published article.
* The STAN code is contained in directory _/stan_
* The files used on the supercomputer are in directory _/supercomp_
* To establish a connection with the manuscript figures:
  + Figure 1 in main manuscript: _R/mod_sel_continuous_summ.R_ produces _weight_by_mod_and_climvar_box_dot.tiff_
  + Figure 2 in main manuscript: _R/mod_sel_continuous_summ.R_ produces _weight_nullmodel_by_response_dot.tiff_
  + Figure S1 in appendix: _R/Figure1.R_, produces file _fig1_vertical.png_
  + Figure S2 in appendix: _R/prior_pc/prior_pc_log_lambda.R_, produces file _normal_y_sim.tiff_
  + Figure S3 in appendix: _R/prior_pc/prior_pc_log_beta.R_, produces file _beta_y_sim.tiff.tiff_
  + Figure S4 in appendix: _R/prior_pc/prior_pc_log_gamma.R_, produces file _gamma_y_sim.tiff.tiff_
  + Figure S5 in appendix: _R/mod_diagnostics.R_, produces file _prop_issue_mod_vr.tiff_
  + Figure S6 in appendix: _R/mod_diagnostics.R_, produces file _prop_issue_rep_yr.tiff_
  + Figure S7 in appendix: _R/mod_diagnostics.R_, produces file _neff_rhat.tiff_
  + Figure S8 in appendix: _R/check_bayes.R_, produces file _p_val_llam_mod.tiff_
  + Figure S9 in appendix: _R/y_vs_ante_case_study.R_, produces file _Astragalus_thygensis_NEW_year2.png_
  + Figure S10 in appendix: _R/mod_sel_continuous_summ.R_, produces file _lppd_diff_precip_2022.2_order.png_
  + Figure S11 in appendix: _R/mod_sel_continuous_summ.R_, produces file _lppd_diff_airt_2022.2_order.png_
  + Figure S12 in appendix: _R/mod_sel_continuous_summ.R_, produces file _weight_by_mod_resp_box_dot_prec-temp.tiff_
  


# file to (for now) mock the final wrapper file
# 1. Get/Format original data
# 2. Fit climate models
# 3. Check models: diagnostics, Bayes-checks, y~x
# 4. Model selection results
library(tidyverse)

# 1. Get/Format original data --------------------

# scrape CHELSA data
# source()
# scrape PRISM data (post 2013)
source( 'scrape_prism.R' )

# format the raw matrices for Cryptantha flava
source( 'format_raw_cryptantha_flava.R' )


# 2. Prior predictive checks ---------------------------------------------

source( 'prior_pc_beta.R' )
source( 'prior_pc_gamma.R' )
source( 'prior_pc_log_lambda.R' )


# 3. Fit climate models --------------------------------------------------

# NOTE: Example of how to fit the actual models 
  # models were fit on a supercomputer with the code in /simulation directory
source( 'mov_win_vr.R' )


# 4. Check models: diagnostics, Bayes-checks, y~x ----------------------

# Possibly import results in the local repository
source( 'import_results.R' )

# model diagnostics and 
source( 'mod_diagnostics.R' )

# Perform Bayesian checks
source( 'check_bayes.R' )

# Store plots of model predictors versus actual data
source( 'y_vs_ante.R' )

# Produce the figure for the case study
source( 'y_vs_ante_case_study.R' )


# 5. Model selection results -------------------------------------------

# Main model selection results, and associated figures 
source( 'mod_sel_continuous_summ.R' )


# 6. Supplementary material on model plots -------------------------------------

# Write code for model plots
source( 'R/pdfs/plots/write_rmd_linear_figures.R' )

# "Knit" the Rmarkdown files for temperature figures
rmarkdown::render('R/pdfs/plots/ppt_plots.Rmd',
                  output_dir  = 'C:/CODE/climate_drivers_methods/results/pdfs/',#getwd(),
                  output_file = 'ppt_model_plots.pdf')


# "Knit" the Rmarkdown files for temperature figures
rmarkdown::render('R/pdfs/plots/tmp_plots.Rmd',
                  output_dir  = 'C:/CODE/climate_drivers_methods/results/pdfs/',#getwd(),
                  output_file = 'tmp_model_plots.pdf')


# 7. Supplementary material on convergence plots -------------------------------

# Write code for convergence plots
source( 'R/pdfs/converg/write_rmd_coverg_figures.R' )

# "Knit" the Rmarkdown files for temperature figures
rmarkdown::render('R/pdfs/converg/ppt_converg_plots.Rmd',
                  output_dir  = 'C:/CODE/climate_drivers_methods/results/pdfs/',#getwd(),
                  output_file = 'ppt_converg_plots.pdf')

# "Knit" the Rmarkdown files for temperature figures
rmarkdown::render('R/pdfs/converg/tmp_converg_plots.Rmd',
                  output_dir  = 'C:/CODE/climate_drivers_methods/results/pdfs/',#getwd(),
                  output_file = 'tmp_converg_plots.pdf')

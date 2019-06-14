# file to (for now) mock the final wrapper file
# 1. Get/Format original data
# 2. Fit climate models
# 3. Check models: diagnostics, Bayes-checks, y~x
# 4. Model selection results


# 1. Get/Format original data --------------------

# scrape CHELSA data
# source()
# scrape PRISM data (post 2013)
source('scrape_prism.R')

# format the raw matrices for Cryptantha flava
source('format_raw_cryptantha_flava.R')


# 2. Fit climate models --------------------------------------------------

# fit the actual models 
source( 'mov_win_vr.R' )


# 3. Check models: diagnostics, Bayes-checks, y~x ----------------------

source( 'mod_diagnostics.R' )
source( 'check_bayes.R' )
source( 'y_vs_ante.R' )


# 4. Model selection results -------------------------------------------

source( 'mod_sel_continuous_summ.R' )

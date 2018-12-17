source('main_phenofit.R')

# library(future)
# plan(multisession)
################################################################################
## global parameters for check_season
load('data/phenoflux115.rda')
nptperyear <- 365
tidy_fluxGPP() # tidy df_GPP

setting    <- list(mar = c(2, 3, 1, 1), mgp = c(1.2, 0.6, 0))
fig.height <- 225
par(setting)

param_step <- 0.1

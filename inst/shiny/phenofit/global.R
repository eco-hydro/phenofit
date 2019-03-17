# shiny::runApp("inst/shiny/phenofit/")
source('global_shiny.R')
source('global_phenofit.R')


# initial date_range
date_begin <- "2010-01-01"
date_end   <- "2014-12-31"

options_wFUN <- c("wTSM", "wBisquare", "wChen", "wKong")

# library(future)
# plan(multisession)
################################################################################
## global parameters for check_season
load('data/flux115_GPP.rda')
df$y <- df$GPP

nptperyear <- 365
# tidy_fluxGPP() # tidy df_GPP

setting    <- list(mar = c(2, 3, 1, 1), mgp = c(1.2, 0.6, 0))
fig.height <- 200 # pixel
lgd.height <- 20

par(setting)

param_step <- 0.1

#' select_var_VI
select_var_VI <- function(df){
	vars_select(colnames(df), -tidyselect::matches("id|ID|site|scale|date|t|year|doy|DayOfYear|QA|QC|qa|qc|w"))
}

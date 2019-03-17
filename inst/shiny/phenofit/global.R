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

# phenofit parameters
options.phenofit <- list(
    file_site              = "", 
    file_veg_text          = "",
    file_veg_rda           = "",
    nptperyear             = "",

    # growing season dividing parameters
    iters_rough            = 4, 
    max_extend_month_rough = 2, 
    calendarYear           = FALSE, 
    FUN_season             = "season_mov", 
    FUN_rough              = "wWHIT", 
    wFUN_rough             = "wKong", 
    r_max                  = 0.2, 
    r_min                  = 0.0, 
    rtrough_min            = 0.6, 

    # fine fitting parameters
    FUN_fine               = "Beck", 
    wFUN_fine              = "wKong", 
    iters_fine             = 2, 
    max_extend_month_fine  = 2, 
    nextend_fine           = 2, 
    use.rough              = FALSE,

    # phenology extraction
    meths_pheno            = ""   
)


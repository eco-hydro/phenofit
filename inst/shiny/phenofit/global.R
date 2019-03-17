# shiny::runApp("inst/shiny/phenofit/")
source('src_load_input.R')
source('global_phenofit.R')
# library(future)
# plan(multisession)


# initial date_range
date_begin <- "2010-01-01"
date_end   <- "2014-12-31"

options_wFUN <- c("wTSM", "wBisquare", "wChen", "wKong")
options_phenofit <- list(
    nptperyear = 365
) %>% setting.update()


file_rda <- 'data/flux115_GPP.rda' %>% normalizePath()
nptperyear = 365

load(file_rda)
df$y <- df$GPP
################################################################################
## global parameters for check_season

# tidy_fluxGPP() # tidy df_GPP

# setting    <- list(mar = c(2, 3, 1, 1), mgp = c(1.2, 0.6, 0))
# par(setting)

## global parameter for UI
fig.height <- 200 # pixel
lgd.height <- 20

param_step <- 0.1 # for r_max and rtrough_max

load_data <- function(){
    file_veg_rda  <- options$file_rda
    file_veg_text <- options$file_veg_text
    file_site <- options$file_site

    # if (is.null(file_))
}

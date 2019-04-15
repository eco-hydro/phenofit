# shiny::runApp("inst/shiny/phenofit/")
source('src_load_input.R')
library(shinyFiles)
library(jsonlite)

# library(future)
# plan(multisession)

check_setting <- function(){
    file_json <- "perference/phenofit_setting.json"
    # print(file.exists(json))
    if (!dir.exists(dirname(file_json)))
        dir.create(dirname(file_json), recursive = TRUE)

    if (file.exists(file_json)) {
        options <- read_json(file_json) %>% map(unlist)
    } else {
        options <- setting.get()
    }

    if (with(options, !(check_file(file_veg_text) || check_file(file_veg_rda)))) {
        file_rda <- 'data/flux115_GPP.rda' %>% normalizePath()
        options$file_veg_rda <- file_rda
        options$nptperyear   <- 365
        options$var_y        <- "GPP"
    }
    options
}


## MAIN SCRIPTS ----------------------------------------------------------------
# initial date_range
date_begin <- "2010-01-01"
date_end   <- "2014-12-31"

options_wFUN <- c("wTSM", "wBisquare", "wChen", "wKong")
# load(file_rda)
# df$y <- df$GPP

################################################################################
## global parameters for check_season
# tidy_fluxGPP() # tidy df_GPP

par_setting    <- list(mar = c(2, 3, 1, 1), mgp = c(1.2, 0.6, 0))
# par(setting)

## global parameter for UI
fig.height <- 200 # pixel
lgd.height <- 20

param_step <- 0.1 # for r_max and rtrough_max

options <- check_setting()
dataIN <- phenofit_loaddata(options)

sites <- dataIN$sites
if (!(options$site %in% sites)){
    options$site <- sites[1]
}

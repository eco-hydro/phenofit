# source('inst/shiny/check_season/global.R')
library(phenofit)
library(shiny)
library(DT)
library(data.table)
library(magrittr)

getSiteData  <- function(df, sitename){
    df[site == sitename, .(t = date, y = GPP_DT, w = 1)] #%T>% plotdata(365)
}

get_input <- function(df, st, sitename){
    sp       <- st[site == sitename]
    south    <- sp$lat < 0
    titlestr <- with(sp, sprintf("[%3d] %s, %s, lat = %.2f", ID, site, IGBP, lat))

    d  <- df[site == sitename, .(t = date, y = GPP_DT, w = 1)] #%T>% plotdata(365)
    # d_obs[site == sitename, .(t = date, y = GPP_DT, w = 1)] %>% plotdata(365)
    d_new <- add_HeadTail(d, south = south)
    INPUT <- do.call(check_input, d_new)
    # d <- d_obs[site == sitename, ]

    INPUT$south <- south
    INPUT$titlestr <- titlestr
    # list(INPUT = INPUT, plotdat = d)
    INPUT
}

check_season <- function(INPUT,
                         FUN_season = c("season", "season_3y"),
                         FUN_fit = "wWHIT",
                         wFUN = "wTSM",
                         lambda = 1000,
                         iters = 3,
                         IsPlot = F, ...) {
    # sitename <- "US-ARM" # "FR-LBr", "ZA-Kru", "US-ARM"

    FUN_season <- get(FUN_season[1])
    wFUN       <- get(wFUN)
    res  <- FUN_season(INPUT, south = INPUT$south, rFUN = get(FUN_fit),
                       wFUN = wFUN,
                     IsPlot = IsPlot,
                     lambda = lambda,
                     iters = iters,
                     minpeakdistance = 30,
                     MaxPeaksPerYear = 3,
                     MaxTroughsPerYear = 4,
                     ypeak_min = 1, ...,
                     IsOnlyPlotbad = FALSE
    )

    if (IsPlot){
        abline(h = 1, col = "red")
        title(INPUT$titlestr)
    }
    return(res)
}

plot_data <- function(d, title){
    par(setting)
    do.call(check_input, d) %>% plotdata()
    mtext(title, side = 2, line = 2, cex = 1.3, font = 2)
}

################################################################################
## global parameters for check_season
nptperyear <- 365
setting    <- list(mar = c(2, 3, 1, 1), mgp = c(1.2, 0.6, 0))
fig.height <- 225
par(setting)

param_step <- 0.1
# sites_grp  <- list(sites_sg, sites_mg) %>%
#     set_names(c("single season", "multiple season"))

# threshold_max = 0.1,
# threshold_min = 0,
# rytrough_max = 0.6,

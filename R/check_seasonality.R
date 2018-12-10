# Check vegetation seasonality
# @inheritParams season
check_seasonality <- function(INPUT, IsPlot = F, plotdat = INPUT, nf = 3, ...){
    # y <- d$EVI
    # t <- as.numeric(d$date - ymd(20000101))
    # w <- d$w
    # nptperyear = 23
    # INPUT <- check_input(t, y, w)
    #
    # pdat <- as.list(d[, .(y = EVI, t = date, w = w)])
    # pdat$ylu <- INPUT$ylu
    fit <- wHANTS(INPUT$y, INPUT$t, INPUT$w, nf = nf, ylu = INPUT$ylu,
                  nptperyear = 23, iters = 3, wFUN = wTSM, wmin = 0.1)
    if (IsPlot){
        plotdata(plotdat, 23)
        colors <- c("red", "blue", "green")
        for (i in 1:(ncol(fit) - 1)){
            lines(plotdat$t, fit[[i+1]], col = colors[i], lwd = 2)
        }
    }
    stat  <- GOF(INPUT$y, dplyr::last(fit), INPUT$w, include.cv = T)
    stat # quickly return
}

# #' REMOVE non-seasonality points according to HANTS smooth result
# rmNonSeasonality <- function(dt, IsPlot = T, file = 'SI.pdf'){
#     if (IsPlot){
#         CairoPDF(file, width = 10, height = 2*6)
#         par(mfrow = c(6, 1), mar = c(2, 3, 2, 1), mgp = c(1.5, 0.6, 0))
#     }

#     sites <- unique(dt$site)
#     n     <- length(sites)
#     res   <- list()
#     stats <- list()

#     for (i in 1:n){
#         runningId(i)
#         sitename  <- sites[i]
#         d         <- dt[site == sitename]
#         IGBP_code <- d$IGBPcode[1]#d$IGBP[1]
#         IGBP_name <- IGBPnames[IGBP_code]
#         lat       <- d$lat[1]

#         # vc <- v_curve(INPUT$y, w = INPUT$w, llas = seq(-2, 2, by = 0.01), d = 2, show = F); lambda <- vc$lambda
#         # lambda <- init_lambda(INPUT$y)
#         # if ( sum(d$w >= 0.5, na.rm = T) > nrow(d) * 0.3){
#             INPUT <- check_input(d$t, d$y, d$w, maxgap = nptperyear / 4, alpha = 0.02)

#             INPUT_SI <- INPUT
#             INPUT_SI$t %<>% {as.numeric(difftime(., ymd(20000101)))}

#             stat       <- check_seasonality(INPUT_SI, IsPlot = FALSE)
#             stat_str   <- stat[c("R2", "NSE", "cv")] %>% {paste(names(.), round(., 3), sep = "=", collapse = ", ")}
#             stats[[i]] <- c(site = sitename, IGBPcode = IGBP_code, stat)
#             str_title  <- sprintf("[%02d] %s, IGBP = %s, %s, lat = %.2f", i, sitename, IGBP_name, stat_str, lat)

#             NSE <- stat[['NSE']]
#             cv  <- stat[['cv']]

#             if (NSE < 0 | (cv < 0.1 & NSE < 0.1)) {
#                 # plot data for check_SI
#                 pdat <- as.list(d[, .(y, t, w)]); pdat$ylu <- INPUT$ylu
#                 temp <- check_seasonality(INPUT_SI, IsPlot, pdat)
#                 if(IsPlot) title(str_title)
#             }
#         # } else{
#         #     message(str_title)
#         # }
#     }
#     if(IsPlot){
#         dev.off()
#         file.show(file)
#     }

#     info <- do.call(rbind, stats) %>% data.table()
#     info[, IGBPname := IGBPnames[IGBPcode]]
#     info <- info[order(IGBPcode, site)]
#     info # quickly return
# }


file <- sprintf("optim_lambda, nptperyear=%d_season_bad_4_si.pdf", nptperyear)

IsPlot <- F

if (IsPlot){
    CairoPDF(file, width = 10, height = 2*6)
    par(mfrow = c(6, 1), mar = c(2, 3, 2, 1), mgp = c(1.5, 0.6, 0))
}

res   <- list()
stats <- list()
for (i in seq_along(sites)){
    runningId(i)
    sitename <- sites[i]#; grp = 1
    d    <- dt[site == sitename]
    # d    <- d[1:(23*3), ]
    lat  <- d$coords_x2[1] #d$lat[1]
    IGBP_code <- d$IGBPcode[1]#d$IGBP[1]
    IGBP_name <- IGBPnames[IGBP_code]

    if ( sum(d$w >= 0.5, na.rm = T) > nrow(d) * 0.3){
        INPUT <- check_input(d$t, d$y, d$w, trim = T, maxgap = nptperyear / 4, alpha = 0.02)

        pdat <- as.list(d[, .(y, t, w)]); pdat$ylu <- INPUT$ylu
        INPUT_SI <- INPUT
        INPUT_SI$t %<>% {as.numeric(difftime(., ymd(20000101)))}

        stat <- check_SI(INPUT_SI, IsPlot, pdat)
        if(IsPlot) title(str_title)

        stats[[i]] <- c(site = sitename, IGBPcode = IGBP_code, IGBPname = IGBP_name,
                        stat)

        # vc <- v_curve(INPUT$y, w = INPUT$w, llas = seq(-2, 2, by = 0.01), d = 2, show = F); lambda <- vc$lambda
        lambda <- init_lambda(INPUT$y)
        str_title <- sprintf("[%02d] %s, IGBP = %s, log10(lambda) = %.3f, lat = %.2f", i, sitename, IGBP_name, log10(lambda), lat)
        # params <- list(INPUT, lambda = lambda, nptperyear, iters = 3, wFUN = wFUN, IsPlot = F,
        #                rymin_less = 0.7, ymax_min = ymax_min,
        #                max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5,
        #                south = lat < 0, plotdat = d)
        #
        # brks <- do.call(season, params) #, ...
        # if (is.null(brks)) next()
        # stat       <- GOF(d$y, brks$whit$iter3) #, INPUT$w
        # stats[[i]] <- c(site = sitename, stat)
        # if (stat[['NASH']] < 0.4){
        #     params$IsPlot <- T
        #     brks <- do.call(season, params)
        #     title(str_title)
        # }
        # res[[i]]  <- listk(site = sitename, lat, IGBP_name, lambda = lambda)#, vc$
    } else{
        message(str_title)
    }
}
dev.off()
# res %<>% set_names(sites)

if(IsPlot) file.show(file)
info <- do.call(rbind, stats) %>% data.table()
info$site     %<>% as.numeric()
info$IGBPcode %<>% as.numeric()

info <- info[order(IGBPcode, site)]
# library(phenofit)
# library(data.table)
# library(lubridate)
# dt <- dt_MOD13A1
# sitename <- dt$site[1]
# d <- dt[site == sitename,]

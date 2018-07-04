#' Check vegetation seasonality
#' @export
check_seasonality <- function(INPUT, IsPlot = F, pdat = INPUT, nf = 3, ...){
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
        plotdata(pdat, 23)
        colors <- c("red", "blue", "green")
        for (i in 1:(ncol(fit) - 1)){
            lines(pdat$t, fit[[i+1]], col = colors[i], lwd = 2)
        }
    }
    stat  <- GOF(INPUT$y, dplyr::last(fit), INPUT$w, include.cv = T)
    stat # quickly return
}

#' REMOVE non-seasonality points according to HANTS smooth result
rmNonSeasonality <- function(dt, IsPlot = T, file = 'SI.pdf'){
    if (IsPlot){
        CairoPDF(file, width = 10, height = 2*6)
        par(mfrow = c(6, 1), mar = c(2, 3, 2, 1), mgp = c(1.5, 0.6, 0))
    }

    sites <- unique(dt$site)
    n     <- length(sites)
    res   <- list()
    stats <- list()

    for (i in 1:n){
        runningId(i)
        sitename  <- sites[i]
        d         <- dt[site == sitename]
        IGBP_code <- d$IGBPcode[1]#d$IGBP[1]
        IGBP_name <- IGBPnames[IGBP_code]
        lat       <- d$lat[1]

        # vc <- v_curve(INPUT$y, w = INPUT$w, llas = seq(-2, 2, by = 0.01), d = 2, show = F); lambda <- vc$lambda
        # lambda <- init_lambda(INPUT$y)
        # if ( sum(d$w >= 0.5, na.rm = T) > nrow(d) * 0.3){
            INPUT <- check_input(d$t, d$y, d$w, trim = T, maxgap = nptperyear / 4, alpha = 0.02)

            INPUT_SI <- INPUT
            INPUT_SI$t %<>% {as.numeric(difftime(., ymd(20000101)))}

            stat       <- check_seasonality(INPUT_SI, IsPlot = FALSE)
            stat_str   <- stat[c("R2", "NSE", "cv")] %>% {paste(names(.), round(., 3), sep = "=", collapse = ", ")}
            stats[[i]] <- c(site = sitename, IGBPcode = IGBP_code, stat)
            str_title  <- sprintf("[%02d] %s, IGBP = %s, %s, lat = %.2f", i, sitename, IGBP_name, stat_str, lat)

            NSE <- stat[['NSE']]
            cv  <- stat[['cv']]

            if (NSE < 0 | (cv < 0.1 & NSE < 0.1)) {
                # plot data for check_SI
                pdat <- as.list(d[, .(y, t, w)]); pdat$ylu <- INPUT$ylu
                temp <- check_seasonality(INPUT_SI, IsPlot, pdat)
                if(IsPlot) title(str_title)
            }
        # } else{
        #     message(str_title)
        # }
    }
    if(IsPlot){
        dev.off()
        file.show(file)
    }

    info <- do.call(rbind, stats) %>% data.table()
    info[, IGBPname := IGBPnames[IGBPcode]]
    info <- info[order(IGBPcode, site)]
    info # quickly return
}

#' 
#' Add one year data in the head and tail
#' 
#' 
#' @param d A data.table, should have \code{t} (compositing date) column 
#' (\code{Date} variable).
#' 
#' @return data.table
add_HeadTail <- function(d, nptperyear = 23){
    # I_beg <- floor(yday(ymd(20000218))/16) # MOD13A1 20000218 is the 4th 16-day
    ntime    <- nrow(d)
    step     <- ceiling(365/nptperyear)

    nmissing_head <- floor( yday(d$t[1])/step )
    nmissing_tail <- nptperyear - floor( yday(last(d$t))/step ) - 1

    date_beg <- first(d$t)
    date_end <- last(d$t)
    
    year_beg <- year(date_beg)
    year_end <- year(date_end)

    md_beg <- date_beg %>% {month(.)*100 + day(.)} # begin date month-day
    md_end <- date_end %>% {month(.)*100 + day(.)} # end   date month-day
    
    # in case of leap year
    if (md_beg == 0229) md_beg = 0228
    if (md_beg == 0229) md_end = 0301
     
    # if tiny missing, than this year is include to extract phenology
    nyear_add_head <- ifelse (nmissing_head < nptperyear/2, 2, 1)  
    nyear_add_tail <- ifelse (nmissing_tail < nptperyear/2, 2, 1)

    # head
    d_head <- d[t >= ymd( (year_beg+1)*1e4 + 0101) & 
                t <  ymd( (year_beg+nyear_add_head)*1e4 + md_beg)]

    d_head$t %<>% `-`(dyears(nyear_add_head)) 

    # tail 
    d_tail <- d[t >  ymd( (year_end-nyear_add_tail)*1e4 + md_end) & 
                t <= ymd( (year_end-1)*1e4 + 1231)]
    d_tail$t %<>% `+`(dyears(nyear_add_tail)) 

    res <- rbind(d_head, d, d_tail)
    res # return
}


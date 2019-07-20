#' backval
#' 
#' Calculate backgroud values for vegetation index.
#' 
#' Night temperature Tn >= Tmin (default 5 degree) defined as raw growing season.  
#' Background value is determined from two neighboring vegetation in raw growing 
#' season by assuming that the background and vegetation abundance could remain 
#' the same during a consecutive two yearperiod. 
#' Details can be seen in Zhang et al., (2015).
#' 
#' @inheritParams check_input
#' @param minT min temperature for growing season.
#' 
#' @return back If back value is NA, it is impossible to extract phenology here.
#' 
#' @keywords internal
#' @note This function only works in every growing season.
#' 
#' @references
#' 1. Zhang, X., 2015. Reconstruction of a complete global time series of daily 
#'      vegetation index trajectory from long-term AVHRR data. Remote Sens. Environ. 
#'      156, 457–472. https://doi.org/10.1016/j.rse.2014.10.012. \cr
#' 2. Zhang, Y., Xiao, X., Jin, C., Dong, J., Zhou, S., Wagle, P., Joiner, 
#'      J., Guanter, L., Zhang, Y., Zhang, G., Qin, Y., Wang, J., Moore, B., 2016. 
#'      Consistency between sun-induced chlorophyll fluorescence and gross primary 
#'      production of vegetation in North America. Remote Sens. Environ. 
#'      183, 154–169. https://doi.org/10.1016/j.rse.2016.05.015.
#' 
#' @export
backval <- function(y, t, w, Tn, minT = 5, nptperyear, ...){
    # get median of 5 smallest values in y[index]
    getBack <- function(y, index){
        yi <- y[index]
        n  <- length(yi)

        i_end <- ifelse(n <= 5, n, 5)

        # if i_end = 0, yi[1:0] = NA
        return( median(sort(yi)[1:i_end]) ) #fixed 20180913
    }

    n <- length(y)
    i_lowT <- (Tn < minT)
    n_lowT <- length(which(i_lowT))

    i_good    <- w >= 1
    i_margin  <- w >= 0.5
    
    if (n_lowT == 0){ 
        # Low latitude region, no winter identified, backval = averaging the five 
        # smallest EVI2 values with good quality during a period of two years
        # (Zhang, 2015).
        back <- getBack(y, i_good)
    } else if (n - n_lowT < max(5, nptperyear/12*2)){ 
        # Warming temperature less 2 month (e.g. high latitude, boreal region)
        # If nighttimetemperature is always below 10°C, the growing season is 
        # set to be June to August (Zhang et al., 2016)
        month    <- month(t)
        i_summer <- month %in% c(6, 7, 8)
        
        # In the so cold region, good values also should be very limited.
        I <- i_summer & i_good
        if (sum(I) < 5) I <- i_summer & i_margin
        # if (sum(I) < 5) I <- i_summer
        back <- getBack(y, I)
    } else {
        i_goodT   <- i_good & (!i_lowT) # need to search the nearest value
        # if no enough good values, consider margin values
        if (sum(i_goodT) < 5) {
            i_marginT <- i_margin & (!i_lowT) # need to search the nearest value
            i_goodT   <- i_marginT 
        }

        # enough good values to select
        back <- getBack(y, i_goodT)
    }
    
    return(back)
    # # update data
    # y[y < back] <- back
    # w[y < back] <- 1
}


# nogrowthPolygon <- function(d, Tmin){
#     x <- d$Tn#[100:400]
#     I_ends   <- findBrks(x > Tmin, zero = "-", nups = 1)
#     I_begins <- findBrks(x < Tmin, zero = "-", nups = 1) - 1

#     plot(type = "b", x); grid(); abline(a = Tmin, b = 0, col = "red")
#     points(I_ends, x[I_ends], col = "red", pch = 19)
#     points(I_begins, x[I_begins], col = "blue", pch = 19)

#     nptperyear <- 23
#     polygons <- list()
#     for (i in seq_along(I_ends)){
#         i_beg <- I_ends[i]
#         i_end <- I_begins[which(I_begins > i_beg & I_begins < i_beg + nptperyear*2/3)]

#         if (!is_empty(i_end)){
#             polygons[[i]] <- c(i_beg, -Inf, i_end, -Inf, i_end, Inf,i_beg, Inf,i_beg, -Inf) %>%
#                 matrix(ncol = 2, byrow = TRUE) %>%
#                 set_colnames(c("x", "y")) %>% as.data.frame()
#         }
#     }
#     polygons %<>% set_names(seq_along(.))
#     df_polygon <- melt_list(polygons, "id")
#     df_polygon$date <- d$date[df_polygon$x]
#     return(df_polygon)
# }


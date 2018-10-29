# source('test/07_whit/main_gee_Whittaker.R')

load("data_test/whit_lambda/lambda_formula_v013.rda") # v013
coef_extra <- matrix(c(
    0.831120,  0.035160, 0     , 1.599970, -4.094027 , -0.063533 ,
    0.8209  ,  0       , 0.0041, 1.5008  , -4.0286   , -0.1017   ,
    0.831120, -0.035160, 0     , 1.599970, - 4.094027, -0.063533),
   nrow = 3, byrow = T, dimnames = list(c("v12", "v13", "gee"), NULL))
coefs <- rbind(coef_extra, coef$mean)

# 1.3 global parameters
nptperyear = 23
print  = F
IsPlot = F # for brks

nf = 4
frame = floor(nptperyear/5*2) + 1;# print(frame)

################################################################################
#' fill_missdate
#' fill missing date
fill_missdate <- function(){
    years <- 2000:2018
    doy   <- seq(1, 366, 16)
    date  <- sprintf("%4d%03d", rep(years, each = 23), doy) %>%
        parse_date_time("%Y%j") %>% date()
    if (years[1] == 2000) date <- date[-(1:3)]
    date  <- date[1:(length(date)-12)] # for 2018
    date
}

#' GOF_season3y
#' GOF of season3y object
GOF_season3y <- function(brks2){
    # browser()
    GOF_fun <- function(d){
        vars_iter <- d %>% contain("ziter")
        varnames  <- gsub("z", "", vars_iter)

        res <- map(vars_iter, function(varname){
            GOF(d$y, d[[varname]])
        }) %>%
            do.call(rbind, .) %>%
            as.data.table() %>% cbind(iter = varnames, .)
        res
    }

    all  <- brks2$whit %>% GOF_fun()
    good <- brks2$whit[witer1 >= 1] %>% GOF_fun()

    list(all = all, good = good) %>% melt_list("type")
}

#' rough_fitting
#' @param lambda Unless lambda is constant, lambda should be null.
rough_fitting <- function(sitename, df, st, .FUN = wWHIT, lambda = NULL,
    IsPlot = FALSE, print = FALSE, ...)
{
    sp    <- st[site == sitename, ] # station point
    d     <- df[site == sitename, ] # get the first site data
    south <- sp$lat < 0
    titlestr <- with(sp, sprintf('[%03d,%s] %s, ', ID, as.character(site), IGBPname))
    cat(titlestr, "\n")

    tryCatch({
        # fill missing values
        # date : image date
        # t    : compositing date
        if (is_flux) {
            d <- merge(d, data.table(date), by = "date", all = T)
        } else {
            d <- merge(d, data.table(t = date), by = "t", all = T)
        }
        ############################################################################
        dnew  <- add_HeadTail(d, south, nptperyear)
        # 1. Check input data and initial parameters for phenofit
        INPUT <- check_input(dnew$t, dnew$y, dnew$w,
                             nptperyear, south = south,
                             maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
        INPUT$y0 <- dnew$y

# browser()

        ## 20180819 fixed lambda bug, lambda will overwrite new lambda in season_3y
        # if (is.null(lambda)) lambda <- init_lambda(INPUT$y)#*2w
        brks2 <- season_3y(INPUT,
                           rFUN = .FUN,
                           lambda = lambda, nf = nf, frame = frame,
                           adj.param = adj.param, # default is true
                           plotdat = d, IsPlot = IsPlot, print = print,
                           titlestr = titlestr, IsOnlyPlotbad = F, ...)
        brks2$GOF <- GOF_season3y(brks2)
        brks2
    # }, error = function(e){
    #     message(sprintf("[e]: %s, %s", titlestr, e$message))
    }, warning = function(w){
        message(sprintf("[w]: %s, %s", titlestr, w$message))
    })
}

# init_lambda
# global param: param
init_lambda <- function(y){
    # print("running here")
    y        <- y[!is.na(y)] #rm NA values
    mean     <- mean(y)
    sd       <- sd(y)
    cv       <- sd/mean
    skewness <- skewness(y, type = 2)
    kurtosis <- kurtosis(y, type = 2)

    # lambda was transformed by log10
    # lambda   <- 0.555484 + 1.671514*mean - 3.434064*sd - 0.052609*skewness + 0.009057*kurtosis
    # lambda   <- 0.555465 + 1.501239*mean - 3.204295*sd - 0.031902*skewness # Just three year result
    # lambda <- 0.831120 + 1.599970*mean - 4.094027*sd - 0.035160*cv - 0.063533*skewness # all year

    ## update 2018-07-31
    # lambda <- 0.7835 +1.5959*mean -4.0371*sd +0.0048*cv -0.1032*skewness +0.0036*kurtosis # yearly
    lambda <- 0.8209 +1.5008*mean -4.0286*sd -0.1017*skewness -0.0041*kurtosis            # 4-year

    lambda <- param$`(Intercept)` + param$mean*mean + param$sd*sd +
        param$skewness*skewness + param$kurtosis*kurtosis + param$cv*cv
    # lambda <- 0.817783 + 1.803588*mean - 4.263469*sd - 0.038240*cv - 0.066914*skewness - 0.011289*kurtosis  #3y
    return(10^lambda)
}


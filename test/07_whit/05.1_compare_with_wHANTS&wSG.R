source('test/stable/load_pkgs.R')

#' @param lambda Unless lambda is constant, lambda should be null.
rough_fitting <- function(sitename, df, st, .FUN = whitsmw2, lambda = NULL){
    sp    <- st[site == sitename, ] # station point
    d     <- df[site == sitename, ] # get the first site data
    titlestr <- with(sp, sprintf('[%03d,%s] %s, ', ID, as.character(site), IGBPname))
    cat(titlestr, "\n")

    tryCatch({
        # fill missing values
        d <- merge(d, data.table(t = date), by = "t", all = T)
        ############################################################################
        dnew  <- add_HeadTail(d)
        # 1. Check input data and initial parameters for phenofit
        INPUT <- check_input(dnew$t, dnew$y, dnew$w, maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
        INPUT$y0 <- dnew$y
        # 2. The detailed information of those parameters can be seen in `season`.

        ## 20180819 fixed lambda bug, lambda will overwrite new lambda in season_3y
        # if (is.null(lambda)) lambda <- init_lambda(INPUT$y)#*2w
        brks2 <- season_3y(INPUT, nptperyear, south = sp$lat[1] < 0, FUN = .FUN,
                           lambda = lambda, nf = nf, frame = frame,
                           plotdat = d, IsPlot = IsPlot, print = print,
                           titlestr = titlestr, partial = F)
        brks2
    }, error = function(e){
        message(sprintf("error: %s, %s", titlestr, e$message))
        NULL
    })
}


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

###############################################################################
load("data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda")
load("data_test/lambda_formula.rda")

sites    <- unique(df$site) %>% set_names(., .)
sitename <- sites[1]

# 1.1 lambda formula coefs
coef_extra <- matrix(c( 0.831120, 0.035160, 0, 1.599970, -4.094027, -0.063533,
   0.8209, 0, 0.0041, 1.5008, -4.0286, -0.1017,
   0.831120, -0.035160, 0, 1.599970, - 4.094027, -0.063533),
   nrow = 3, byrow = T, dimnames = list(c("v12", "v13", "gee"), NULL))
coefs <- rbind(coef_extra, coef$mean)

# 1.2 fill missing values
years <- 2000:2018
doy   <- seq(1, 366, 16)
date  <- sprintf("%4d%03d", rep(years, each = 23), doy) %>% parse_date_time("%Y%j") %>% date()
if (years[1] == 2000) date <- date[-(1:3)]
date  <- date[1:(length(date)-11)] # for 2018

#
# 1.3 global parameters
nptperyear = 23
print  = F
IsPlot = F # for brks

nf = 4
frame = floor(nptperyear/5*2) + 1;# print(frame)
noise_percs = c(0.1, 0.3, 0.5)
noise_perc  = 0.3

methods  <- c("wHANTS", "sgfitw", "whitsmw2", "whitsmw2")
methods2 <- c("wHANTS", "wSG", "wWH", "wWH2")

lst <- list()
k = 2

source("R/season_3y.R")
source("R/curvefits.R")

for (k in 1:nrow(coefs)){
# for (k in 3){
    # noise_perc <- noise_percs[k]
    # df <- select_valid(df, noise_perc = noise_perc)[, 1:10]
    runningId(k, prefix = "k | " )
    pattern <- rownames(coefs)[k]
    param   <- as.list(coefs[k, ])
    ############################################################################
    ############################################################################
    for (i in 3){ # only wWH2 this time
        method <- methods[i] #"sgfitw", "whitsmw2" and "wHANTS".
        FUN    <- get(method, envir = as.environment("package:phenofit"))

        if (IsPlot){
            file <- sprintf("%s_%s_%2d%%.pdf", prefix, methods2[i], noise_perc*100)
            CairoPDF(file, 8, 10)
            par(mar = c(1.5, 2, 2, 1), mgp = c(3, 0.6, 0), mfrow = c(5, 1), ann = F)
        }

        lambda <- NULL
        if (i >= 4) lambda <- 2
        # a <- llply(sites[23:100], rough_fitting,
        #              df = df, st = st, FUN = get(method), lambda = lambda,
        #              .progress = "text") #lst[[i]]
        # rough_fitting(sitename, df, st, FUN = whitsmw2, lambda = lambda)

        # outdir <- sprintf("result/whit_lambda/valid/%s_%2d%%", methods2[i], noise_perc*100)
        # outdir <- sprintf("/flush1/kon055/result/whit_lambda/%s_0", methods2[i])
        outdir <- sprintf("/flush1/kon055/result/whit_lambda2/%s", pattern)
        temp <- par_sbatch(sites, rough_fitting,
                           df = df, st = st, .FUN = get(method), lambda = lambda,
                           return.res = F, Save = T, outdir = outdir)

        if (IsPlot) dev.off()
        # brks2 <- season_3y(INPUT, nptperyear, south = sp$lat[1] < 0, FUN =FUN,
        #                    nf = nf, frame = frame,
        #                    IsPlot = IsPlot, print = print, partial = F)
    }
}

# lst %<>% set_names(methods2)
# save(df, lst, file = outfile)

# i = 1
# if (i == 1){
#     infile <- file_flux
#     st     <- fread(file_st_flux)
#     outfile <- sprintf("data_test/phenoflux166_rough_val_%02d%%.rda", noise_perc*100)
#     prefix  <- "phenoflux"
# } else{
#     infile <- file_cam
#     st     <- fread(file_st_cam)
#     outfile <- sprintf("data_test/phenocam133_rough_val_%02d%%.rda", noise_perc*100)
#     prefix  <- "phenocam"
# }

# sitename  <- sites[21]
# lst <- llply(sites, rough_fitting,
#              df = df, st = st, FUN = get(method),
#              .progress = "text")
# rough_fitting(sitename, df, st, FUN = get(method))
# df <- fread(infile) # , strip.white = T
# setkeyv(df, c("site", "t"))
# no noise in this version

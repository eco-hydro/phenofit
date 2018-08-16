source('test/stable/load_pkgs.R')

nptperyear = 23
print  = F
IsPlot = F # for brks

nf = 4
frame = floor(nptperyear/5*2) + 1;# print(frame)
noise_perc = 0.3

###############################################################################

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

# fill missing values
years <- 2000:2018
doy   <- seq(1, 366, 16)
date  <- sprintf("%4d%03d", rep(years, each = 23), doy) %>% parse_date_time("%Y%j") %>% date()
if (years[1] == 2000) date <- date[-(1:3)]
date  <- date[1:(length(date)-11)] # for 2018

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
        if (is.null(lambda)) lambda <- init_lambda(INPUT$y)#*2w

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

# sitename  <- sites[21]

# lst <- llply(sites, rough_fitting,
#              df = df, st = st, FUN = get(method),
#              .progress = "text")
# rough_fitting(sitename, df, st, FUN = get(method))

methods  <- c("wHANTS", "sgfitw", "whitsmw2", "whitsmw2")
methods2 <- c("wHANTS", "wSG", "wWH", "wWH2")

lst <- list()

noise_percs = c(0.1, 0.3, 0.5)
k = 2

for (k in 3){
    runningId(k)

    ############################################################################
    noise_perc <- noise_percs[k]
    # df <- fread(infile) # , strip.white = T
    load("data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda")
    # setkeyv(df, c("site", "t"))

    df <- select_valid(df, noise_perc = noise_perc)[, 1:10]
    sites <- unique(df$site) %>% set_names(., .)
    ############################################################################
    
    for (i in 1:4){
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

        outdir <- sprintf("result/whit_lambda/valid/%s_%2d%%", methods2[i], noise_perc*100)
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

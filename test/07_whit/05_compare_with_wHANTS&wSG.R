source('test/stable/load_pkgs.R')

nptperyear = 23
print  = F
IsPlot = F # for brks

nf = 4
frame = floor(nptperyear/5*2) + 1; print(frame)
noise_perc = 0.3

###############################################################################

i = 1
if (i == 1){
    infile <- file_flux
    st     <- fread(file_st_flux)
    outfile <- sprintf("data_test/phenoflux166_rough_val_%02d%%.rda", noise_perc*100)
    prefix  <- "phenoflux"
} else{
    infile <- file_cam
    st     <- fread(file_st_cam)
    outfile <- sprintf("data_test/phenocam133_rough_val_%02d%%.rda", noise_perc*100)
    prefix  <- "phenocam"
}

###############################################################################
select_validI <- function(Id, perc = 0.2) {
    n <- length(Id)
    set.seed(100)
    I <- sample(1:n, n*perc)
    unique(Id[I])
}
# only select good points (with the percentage of noise_perc)
select_valid <- function(df, noise_perc = 0.3, group = F){
    df     <- unique(df) # sometimes data is duplicated at begin and end.
    if (class(df$SummaryQA[1]) != "factor") {
        df[, `:=`( SummaryQA = factor(SummaryQA, qc_levels))]
    }
    if (class(df$t[1]) != "Date") df[, `:=`( t = ymd(t))]

    ## 1.1 get range for every site and divide grp
    alpha = 0.05
    d_range <- df[SummaryQA == "good", .(ymin = quantile(y, alpha/2, na.rm = T),
                                         ymax = quantile(y, 1-alpha/2, na.rm = T)), .(site)] %>%
        plyr::mutate(A = ymax - ymin,
                     brk_min = ymin + 0.2*A,
                     brk_max = ymin + 0.6*A) %>% .[, c(1, 4:6)]
    if (!group) d_range <- d_range[, .(site, A)]

    df <- merge(df, d_range, all.x = T) # rm points has no good points
    df[, `:=`(Id = 1:.N, y0 = y, w0 = w, I_valid = 0)]

    if (group){
        df[, grp:=1]
        df[y >= brk_max, grp := 2]
        df[y <= brk_min, grp := 0]
    }

    # 1. select cross-validation points ---------------------------------------

    ## 1.2 select validation points
    # & grp == 1
    I <- df[SummaryQA == "good", select_validI(Id, noise_perc), .(site)]$V1

    set.seed(I[1])
    desc_perc <- runif(length(I), .05, .5)

    ## 1.3 adjust y and w
    # randomly reduced by 5%–50% with of their amplitude; w set to wmin

    # df_val <- df[I, ]
    df[I, `:=`(w=0, y=y-A*desc_perc, I_valid = 1)]
    return(df)
}

rough_fitting <- function(sitename, df, st, FUN = whitsmw2, lambda = NULL){
    sp    <- st[site == sitename, ] # station point
    d     <- df[site == sitename, ] # get the first site data
    titlestr <- with(sp, sprintf('[%03d,%s] %s, ', ID, site, IGBPname))

    dnew  <- add_HeadTail(d)
    # 1. Check input data and initial parameters for phenofit
    INPUT <- check_input(dnew$t, dnew$y, dnew$w, maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
    INPUT$y0 <- dnew$y
    # 2. The detailed information of those parameters can be seen in `season`.
    if (is.null(lambda)) lambda <- init_lambda(INPUT$y)#*2w

    brks2 <- season_3y(INPUT, nptperyear, south = sp$lat[1] < 0, FUN =FUN,
                       lambda = lambda, nf = nf, frame = frame,
                       plotdat = d, IsPlot = IsPlot, print = print,
                       titlestr = titlestr, partial = F)
    brks2
}

# df <- fread(infile) # , strip.white = T
load("data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda")
setkeyv(df, c("site", "t"))

df <- select_valid(df, noise_perc = noise_perc)[, 1:10]
sites <- unique(df$site) %>% set_names(., .)
sitename  <- sites[21]

# lst <- llply(sites, rough_fitting,
#              df = df, st = st, FUN = get(method),
#              .progress = "text")
# rough_fitting(sitename, df, st, FUN = get(method))

methods  <- c("wHANTS", "sgfitw", "whitsmw2", "whitsmw2")
methods2 <- c("wHANTS", "wSG", "wWH", "wWH2")

lst <- list()

noise_percs = c(0.1, 0.3, 0.5)
k = 2

for (k in 1:3){
    noise_perc <- noise_percs[k]
    for (i in 1:4){
        method <- methods[i] #"sgfitw", "whitsmw2" and "wHANTS".
        FUN    <- get(method)

        if (IsPlot){
            file <- sprintf("%s_%s_%2d%%.pdf", prefix, methods2[i], noise_perc*100)
            CairoPDF(file, 8, 10)
            par(mar = c(1.5, 2, 2, 1), mgp = c(3, 0.6, 0), mfrow = c(5, 1), ann = F)
        }

        lambda <- NULL
        if (i >= 4) lambda <- 2
        # lst[[i]] <- llply(sites, rough_fitting,
        #              df = df, st = st, FUN = get(method), lambda = lambda,
        #              .progress = "text")

        outdir <- sprintf("result/whit_lambda/valid/%s_%2d%%", methods2[i], noise_perc*100)
        temp <- par_sbatch(sites, rough_fitting,
                           df = df, st = st, FUN = get(method), lambda = lambda,
                           return.res = F, Save = T, outdir = outdir)

        if (IsPlot) dev.off()
        # brks2 <- season_3y(INPUT, nptperyear, south = sp$lat[1] < 0, FUN =FUN,
        #                    nf = nf, frame = frame,
        #                    IsPlot = IsPlot, print = print, partial = F)
    }
}


# lst %<>% set_names(methods2)


# save(df, lst, file = outfile)

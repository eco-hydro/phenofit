source("test/load_pkgs.R")

################################################################################
# lambda     <- 5    # Whittaker parameter
ypeak_min    <- 0.1  # the maximum ymax shoud be greater than `ymax_min`
rtrough_max <- 0.8  # trough < ymin + A*rymin_less
nptperyear   <- 23   # How many points for a single year
wFUN         <- wTSM # Weights updating function, could be one of `wTSM`, 'wBisquare', `wChen` and `wSELF`.

print = T

main <- function(df, st, outdir){
    # prefix <- str_extract(infile, "\\w*(?=_MOD)")
    # outdir <- paste0("result/", prefix)

    # df     <- fread(infile) # , strip.white = T
    # df     <- unique(df) # sometimes data is duplicated at begin and end.
    # df[, `:=`( t = ymd(t),
    #     SummaryQA = factor(SummaryQA, qc_levels))]

    sites <- unique(df$site) %>% set_names(., .)
    # sitename  <- df$site[1]
    # d         <- df[site == sites[1], ] # get the first site data

    # sitename = sites[1]; a <- get_phenofit(sitename, df, st, prefix)[1:100]
    par_sbatch(sites, get_phenofit, df=df, st=st, IsPlot = F, #prefix_fig=prefix,
        Save=T, outdir=outdir)
}

# fill missing values
years <- 2000:2018
doy   <- seq(1, 366, 16)
date  <- sprintf("%4d%03d", rep(years, each = 23), doy) %>% parse_date_time("%Y%j") %>% date()
if (years[1] == 2000) date <- date[-(1:3)]
date  <- date[1:(length(date)-11)] # for 2018

# site = unique(df_org$site)
# temp = data.table(t = date, site = rep(site, rep(length(date), length(site))))

# # t as here is image date, other than pixel data.
# df_org = merge(df_org, temp, by = c("t", "site"), all = T) # fill missing values


noise_percs = c(0.1, 0.3, 0.5)
k = 2

load("data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda")

# too much negative values will lead to fitting failure (negative will be 
# replaced to zero).
for (k in 1){
    runningId(k)

    ############################################################################
    noise_perc <- noise_percs[k]
    # df <- fread(infile) # , strip.white = T
    df_org <- df
    # df_org <- select_valid(df, noise_perc = noise_perc)[, 1:10]

    # setkeyv(df, c("t", "site"))
    site = unique(df_org$site)
    temp = data.table(t = date, site = rep(site, rep(length(date), length(site))))

    # t as here is image date, other than pixel data.
    df_org = merge(df_org, temp, by = c("t", "site"), all = T) # fill missing values
    # outdir <- sprintf("/flush1/kon055/result/valid/%s_%2d%%", "phenofit", noise_perc*100)
    outdir <- sprintf("/flush1/kon055/result/valid/%s", "phenofit")
    main(df_org, st, outdir)
    # main(file_flux, st_flux)
}

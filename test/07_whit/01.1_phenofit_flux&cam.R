source("test/load_pkgs.R")
source("test/07_whit/dat_flux&cam_phenofit.R")
# source("R/plot_phenofit.R")

################################################################################
# lambda     <- 5    # Whittaker parameter
ypeak_min    <- 0.1  # the maximum ymax shoud be greater than `ymax_min`
rtrough_max <- 0.8  # trough < ymin + A*rymin_less
nptperyear   <- 23   # How many points for a single year
wFUN         <- wTSM # Weights updating function, could be one of `wTSM`, 'wBisquare', `wChen` and `wSELF`.

print = T


# fill missing values
years <- 2000:2018
doy   <- seq(1, 366, 16)
date  <- sprintf("%4d%03d", rep(years, each = 23), doy) %>% parse_date_time("%Y%j") %>% date()
if (years[1] == 2000) date <- date[-(1:3)]
date  <- date[1:(length(date)-11)] # for 2018

site = unique(df_org$site)
temp = data.table(t = date, site = rep(site, rep(length(date), length(site))))

# t as here is image date, other than pixel data.
df_org = merge(df_org, temp, by = c("t", "site"), all = T) # fill missing values
df <- df_org
################################################################################
outdir <- sprintf("/flush1/kon055/result/valid/flux/%s", "phenofit_0%")

# main <- function(infile, st){
    # prefix <- str_extract(infile, "\\w*(?=_MOD)")
    prefix <- "fluxcam"

    # df     <- fread(infile) # , strip.white = T
    df     <- unique(df) # sometimes data is duplicated at begin and end.
    df[, `:=`( t = ymd(t),
               SummaryQA = factor(SummaryQA, qc_levels))]

    sites <- unique(df$site) %>% set_names(., .)
    sitename  <- df$site[1]
    d         <- df[site == sites[1], ] # get the first site data

    # sitename = sites[1]; a <- get_phenofit(sitename, df, st, prefix)
    par_sbatch(sites, get_phenofit, df=df, st=st, prefix_fig=prefix,
               Save=T, outdir=outdir)
# }
#
# main(file_cam, st_cam)
# main(file_flux, st_flux)
# check phenocam 85 IGBPname

# phenofit::merge_pdf('phenofit_flux166_MOD13A1_v1.pdf', indir = 'Figure/', 'phenoflux166.*.pdf')
# phenofit::merge_pdf('phenofit_cam133_MOD13A1_v1.pdf', indir = 'Figure/', 'phenocam133.*.pdf')

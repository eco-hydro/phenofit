source('test/stable/load_pkgs.R')
source("test/07_whit/dat_flux&cam_phenofit.R")
# source("R/plot_phenofit.R")

################################################################################
# lambda     <- 5    # Whittaker parameter
ypeak_min    <- 0.1  # the maximum ymax shoud be greater than `ymax_min`
rytrough_max <- 0.8  # trough < ymin + A*rymin_less
nptperyear   <- 23   # How many points for a single year
wFUN         <- wTSM # Weights updating function, could be one of `wTSM`, 'wBisquare', `wChen` and `wSELF`.

print = T
outdir <- sprintf("%sresult/valid/flux/%s", dir_flush, "phenofit_0%")


df <- df_org

# main <- function(infile, st){
    # prefix <- str_extract(infile, "\\w*(?=_MOD)")
    prefix <- "fluxcam"
    outdir <- paste0("result/", prefix)

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

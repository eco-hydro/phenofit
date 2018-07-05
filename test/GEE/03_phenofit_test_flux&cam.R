library(phenofit)
source('test/stable/load_pkgs.R')
# source("R/plot_phenofit.R")

################################################################################
if (.Platform$OS.type == "windows"){
    st_flux <- fread("F:/Github/MATLAB/PML/data/flux_166.csv")
    st_cam  <- fread("D:/SciData/PhenoCam/PhenoCam_V1_1511/phenocam133_site.csv")

    st_flux <- st_flux[, .(ID = 1:.N, site, lat, lon = long, IGBPname = IGBP)]
    st_cam  <- st_cam[,  .(ID = 1:.N, site = sitename, lat, lon, primary_veg_type, secondary_veg_type,
               IGBPname = IGBPnames_006[landcover_igbp])]
    fwrite(st_cam , file_st_cam)
    fwrite(st_flux, file_st_flux)
} else if (.Platform$OS.type == "unix"){
    st_cam  <- fread(file_st_cam)
    st_flux <- fread(file_st_flux)
}

if (!file.exists(file_flux)){
    infile_flux  <- "D:/Document/GoogleDrive/phenofit/data/fluxsites212_MOD13A1_006_0m_buffer.csv"
    tidyMOD13INPUT_gee(infile_flux) %>% merge(st_flux[, .(site)]) %>% fwrite(file_flux)
}

if (!file.exists(file_cam)){
    infile_cam  <- "D:/Document/GoogleDrive/phenofit/data/phenocam133_MOD13A1_006_0m_buffer.csv"
    tidyMOD13INPUT_gee(infile_cam) %>% merge(st_cam[, .(site)]) %>% fwrite(file_cam)
}

################################################################################
# lambda     <- 5    # Whittaker parameter
ypeak_min    <- 0.1  # the maximum ymax shoud be greater than `ymax_min`
rytrough_max <- 0.8  # trough < ymin + A*rymin_less
nptperyear   <- 23   # How many points for a single year
wFUN         <- wTSM # Weights updating function, could be one of `wTSM`, 'wBisquare', `wChen` and `wSELF`.


print = T

main <- function(infile, st){
    prefix <- str_extract(infile, "\\w*(?=_MOD)")
    outdir <- paste0("result/", prefix)

    df     <- fread(infile, strip.white = F)
    df     <- unique(df) # sometimes data is duplicated at begin and end.
    df[, `:=`( t = ymd(t),
        SummaryQA = factor(SummaryQA, qc_levels))]

    sites <- unique(df$site) %>% set_names(., .)
    sitename  <- df$site[1]
    d         <- df[site == sites[1], ] # get the first site data

    # sitename = sites[1]; a <- get_phenofit(sitename, df, st, prefix)
    par_sbatch(sites, get_phenofit, df=df, st=st, prefix_fig=prefix,
        save=T, outdir=outdir)
}

main(file_cam, st_cam)
main(file_flux, st_flux)
# check phenocam 85 IGBPname

source("test/load_pkgs.R")
source("test/07_whit/main_phenofit_test.R")

GOF_fineFitting <- function(fit){
    d <- getFittings2(fit)

    # For fine curve fitting, y alreay in `d`. But the `d` has been modified by
    # `check_input`.
    
    by <- .(site, iters, meth) %>% names()
    by <- intersect(by, colnames(d))
    # by <- ""

    all  <- d[      , as.list(GOF_extra2(y, value)), by]
    good <- d[w == 1, as.list(GOF_extra2(y, value)), by]
    info <- listk(all, good) %>% melt_list("type")

    perc_good <- d[, .(perc_good = sum(w == 1)/.N), by]
    info <- merge(perc_good, info)
    info
}

################################################################################

indir <- "result/gee_whittaker/valid (20181024)/phenofit_0%/" %>% paste0(dir_flush, .)
# indir <- "result/gee_whittaker/valid (20180910)/phenofit_0%/" %>% paste0(dir_flush, .)

files <- dir(indir, "*.RDS", full.names = T)
file  <- files[1]
print(files)

outdir <- paste0(dir_flush, "result/gee_whittaker")
prefix <- "phenofit_st1e3 (20180910)"

res <- par_sbatch(files, function(file){
    lst <- readRDS(file) %>% rm_empty()
    res <- llply(lst, GOF_fineFitting, .progress = "text")
}, return.res = F, Save = T, outdir = outdir, prefix = prefix)

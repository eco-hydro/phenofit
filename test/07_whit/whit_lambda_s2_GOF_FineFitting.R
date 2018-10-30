source('test/stable/load_pkgs.R')
source("test/07_whit/main_phenofit_test.R")


indir <- "V:/result/gee_whittaker/valid (20181024)/phenofit_0%/"
files <- dir(indir, "*.RDS", full.names = T)

GOF_fineFitting <- function(fit){
    d <- getFittings2(fit)
    d[, as.list(GOF_extra2(y, value)), .(iters, meth)]
}

res <- llply(files, function(file){
    lst <- readRDS(file) %>% rm_empty()
    llply(lst, GOF_fineFitting)
})
